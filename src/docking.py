"""
Module for running AutoDock Vina molecular docking.
Handles execution via subprocess, output parsing, and batch processing.
"""
import subprocess
import time
import pathlib
from pathlib import Path
import pandas as pd
import re
import logging
import json
from typing import Dict, Any, List, Optional, Tuple, Union
from tqdm import tqdm

logger = logging.getLogger(__name__)

def check_vina_installation(vina_executable: str = 'vina') -> str:
    """
    Verifies Vina is accessible. Runs `vina --version` via subprocess.
    
    Args:
        vina_executable: Path or command for the AutoDock Vina executable.
        
    Returns:
        Version string of the installed Vina.
        
    Raises:
        RuntimeError: If Vina cannot be executed or fails.
    """
    try:
        result = subprocess.run([vina_executable, '--version'], capture_output=True, text=True, check=True)
        version = result.stdout.strip()
        logger.info(f"AutoDock Vina found: {version}")
        return version
    except FileNotFoundError:
        raise RuntimeError(
            f"Vina executable '{vina_executable}' not found in PATH or at specified location.\n"
            "Please ensure AutoDock Vina is installed and accessible."
        )
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Error running Vina: {e.stderr}")

def parse_vina_output(output_text: str) -> List[Dict[str, float]]:
    """
    Parses Vina stdout text output to extract binding affinity scores and RMSD values.
    
    Args:
        output_text: Standard output from Vina.
        
    Returns:
        List of dictionaries containing mode, affinity, rmsd_lb, rmsd_ub.
    """
    results = []
    # Match lines like '   1         -8.1      0.000      0.000'
    pattern = re.compile(r'\s*(\d+)\s+([-\.\d]+)\s+([-\.\d]+)\s+([-\.\d]+)')
    
    in_results_section = False
    for line in output_text.split('\n'):
        if "mode |   affinity | dist from best mode" in line:
            in_results_section = True
            continue
        if in_results_section and "---" in line:
            continue
        if in_results_section:
            match = pattern.match(line)
            if match:
                mode = int(match.group(1))
                affinity = float(match.group(2))
                rmsd_lb = float(match.group(3))
                rmsd_ub = float(match.group(4))
                results.append({
                    'mode': mode,
                    'affinity': affinity,
                    'rmsd_lb': rmsd_lb,
                    'rmsd_ub': rmsd_ub
                })
            elif line.strip() == "":
                break # End of table
    return results

def parse_pdbqt_results(pdbqt_path: Union[str, Path]) -> List[Dict[str, Any]]:
    """
    Parses a multi-model PDBQT output file to extract individual poses and their scores from REMARK lines.
    
    Args:
        pdbqt_path: Path to the output PDBQT file.
        
    Returns:
        List of dictionaries with pose information.
    """
    path = Path(pdbqt_path)
    if not path.exists():
        return []

    poses = []
    current_pose = None
    with open(path, 'r') as f:
        for line in f:
            if line.startswith("MODEL"):
                current_pose = {'model': int(line.split()[1]), 'affinity': None, 'rmsd_lb': None, 'rmsd_ub': None, 'lines': []}
            elif line.startswith("REMARK VINA RESULT:"):
                parts = line.split()
                if len(parts) >= 6:
                    current_pose['affinity'] = float(parts[3])
                    current_pose['rmsd_lb'] = float(parts[4])
                    current_pose['rmsd_ub'] = float(parts[5])
            
            if current_pose is not None:
                current_pose['lines'].append(line)
            
            if line.startswith("ENDMDL"):
                if current_pose:
                    poses.append(current_pose)
                current_pose = None
                
    return poses

def dock_single(receptor_pdbqt: Union[str, Path], 
                ligand_pdbqt: Union[str, Path], 
                center: Tuple[float, float, float], 
                box_size: Tuple[float, float, float], 
                output_path: Union[str, Path], 
                exhaustiveness: int = 32, 
                num_modes: int = 9, 
                energy_range: float = 3.0, 
                seed: int = 42, 
                cpu: int = 4, 
                vina_executable: str = 'vina') -> Dict[str, Any]:
    """
    Docks a single ligand to the receptor using AutoDock Vina.
    
    Args:
        receptor_pdbqt: Path to receptor PDBQT.
        ligand_pdbqt: Path to ligand PDBQT.
        center: (x, y, z) coordinates for the center of the grid box.
        box_size: (x, y, z) dimensions of the grid box.
        output_path: Path for the output PDBQT.
        exhaustiveness: Exhaustiveness of the global search.
        num_modes: Maximum number of binding modes to generate.
        energy_range: Maximum energy difference between the best binding mode and the worst one displayed.
        seed: Random seed.
        cpu: Number of CPUs to use.
        vina_executable: Path to Vina executable.
        
    Returns:
        Dictionary with docking results.
    """
    receptor_pdbqt = Path(receptor_pdbqt)
    ligand_pdbqt = Path(ligand_pdbqt)
    output_path = Path(output_path)
    
    result_dict = {
        'scores': [],
        'best_score': None,
        'output_path': str(output_path),
        'runtime_seconds': 0.0,
        'status': 'Failed',
        'error_message': None
    }
    
    if not receptor_pdbqt.exists():
        result_dict['error_message'] = f"Receptor file not found: {receptor_pdbqt}"
        return result_dict
    if not ligand_pdbqt.exists():
        result_dict['error_message'] = f"Ligand file not found: {ligand_pdbqt}"
        return result_dict
        
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    cmd = [
        vina_executable,
        '--receptor', str(receptor_pdbqt),
        '--ligand', str(ligand_pdbqt),
        '--center_x', str(center[0]),
        '--center_y', str(center[1]),
        '--center_z', str(center[2]),
        '--size_x', str(box_size[0]),
        '--size_y', str(box_size[1]),
        '--size_z', str(box_size[2]),
        '--exhaustiveness', str(exhaustiveness),
        '--num_modes', str(num_modes),
        '--energy_range', str(energy_range),
        '--seed', str(seed),
        '--cpu', str(cpu),
        '--out', str(output_path)
    ]
    
    start_time = time.time()
    try:
        process = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        runtime = time.time() - start_time
        result_dict['runtime_seconds'] = runtime
        
        if process.returncode != 0:
            result_dict['error_message'] = f"Vina exited with code {process.returncode}: {process.stderr}"
            return result_dict
            
        parsed_results = parse_vina_output(process.stdout)
        if parsed_results:
            result_dict['scores'] = [res['affinity'] for res in parsed_results]
            result_dict['best_score'] = parsed_results[0]['affinity']
            result_dict['status'] = 'Success'
        else:
            result_dict['error_message'] = "Could not parse Vina output."
            
    except subprocess.TimeoutExpired:
        result_dict['runtime_seconds'] = time.time() - start_time
        result_dict['error_message'] = "Vina execution timed out after 300 seconds."
    except Exception as e:
        result_dict['runtime_seconds'] = time.time() - start_time
        result_dict['error_message'] = str(e)
        
    return result_dict

def dock_library(receptor_pdbqt: Union[str, Path], 
                 ligand_dir: Union[str, Path], 
                 output_dir: Union[str, Path], 
                 center: Tuple[float, float, float], 
                 box_size: Tuple[float, float, float], 
                 compounds_df: pd.DataFrame, 
                 exhaustiveness: int = 32, 
                 num_modes: int = 9, 
                 energy_range: float = 3.0, 
                 seed: int = 42, 
                 cpu: int = 4, 
                 vina_executable: str = 'vina') -> pd.DataFrame:
    """
    Docks a library of prepared ligands.
    
    Args:
        receptor_pdbqt: Path to receptor PDBQT.
        ligand_dir: Directory containing ligand PDBQT files.
        output_dir: Directory for output PDBQT files.
        center: (x, y, z) coordinates.
        box_size: (x, y, z) dimensions.
        compounds_df: DataFrame with 'compound_id' and 'compound_name'.
        exhaustiveness: Exhaustiveness of search.
        num_modes: Maximum modes.
        energy_range: Energy range.
        seed: Random seed.
        cpu: Number of CPUs.
        vina_executable: Path to Vina executable.
        
    Returns:
        DataFrame with docking results.
    """
    ligand_dir = Path(ligand_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    results = []
    
    for _, row in tqdm(compounds_df.iterrows(), total=len(compounds_df), desc="Docking library"):
        comp_id = row['compound_id']
        comp_name = row.get('name', comp_id)
        
        safe_name = str(comp_name).replace(" ", "_").replace("/", "_").replace("\\", "_")
        if not safe_name: safe_name = str(comp_id)
        
        ligand_file = ligand_dir / f"{safe_name}.pdbqt"
        output_file = output_dir / f"{comp_id}_out.pdbqt"
        
        res = dock_single(
            receptor_pdbqt=receptor_pdbqt,
            ligand_pdbqt=ligand_file,
            center=center,
            box_size=box_size,
            output_path=output_file,
            exhaustiveness=exhaustiveness,
            num_modes=num_modes,
            energy_range=energy_range,
            seed=seed,
            cpu=cpu,
            vina_executable=vina_executable
        )
        
        results.append({
            'compound_id': comp_id,
            'compound_name': comp_name,
            'best_score': res['best_score'],
            'all_scores': res['scores'],
            'runtime_seconds': res['runtime_seconds'],
            'status': res['status'],
            'error_message': res['error_message'],
            'output_pdbqt_path': res['output_path']
        })
        
    return pd.DataFrame(results)

def save_docking_results(results_df: pd.DataFrame, output_path: Union[str, Path]) -> None:
    """
    Saves docking results to CSV, expanding score columns for top poses.
    
    Args:
        results_df: DataFrame containing docking results.
        output_path: Path to save the CSV.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    df = results_df.copy()
    
    # Expand all_scores list into separate columns
    max_poses = df['all_scores'].apply(lambda x: len(x) if isinstance(x, list) else 0).max()
    
    if pd.isna(max_poses) or max_poses == 0:
        logger.warning("No valid scores to expand.")
        df.to_csv(output_path, index=False)
        return

    for i in range(1, int(max_poses) + 1):
        df[f'score_pose_{i}'] = df['all_scores'].apply(
            lambda x: x[i-1] if isinstance(x, list) and len(x) >= i else None
        )
        
    # Drop the original list column to make it flat and clean
    if 'all_scores' in df.columns:
        df = df.drop(columns=['all_scores'])
        
    df.to_csv(output_path, index=False)
    logger.info(f"Docking results saved to {output_path}")

def generate_vina_config(receptor_path: Union[str, Path], 
                         ligand_path: Union[str, Path], 
                         output_path: Union[str, Path], 
                         center: Tuple[float, float, float], 
                         box_size: Tuple[float, float, float], 
                         exhaustiveness: int, 
                         num_modes: int, 
                         energy_range: float, 
                         seed: int, 
                         cpu: int, 
                         config_path: Union[str, Path]) -> Path:
    """
    Writes a Vina config file (.txt) as an alternative to command-line arguments.
    
    Args:
        receptor_path: Path to receptor PDBQT.
        ligand_path: Path to ligand PDBQT.
        output_path: Path to output PDBQT.
        center: (x, y, z) coordinates.
        box_size: (x, y, z) dimensions.
        exhaustiveness: Exhaustiveness of search.
        num_modes: Maximum modes.
        energy_range: Energy range.
        seed: Random seed.
        cpu: Number of CPUs.
        config_path: Path to save the config file.
        
    Returns:
        Path to the generated config file.
    """
    config_path = Path(config_path)
    config_path.parent.mkdir(parents=True, exist_ok=True)
    
    config_content = [
        f"receptor = {receptor_path}",
        f"ligand = {ligand_path}",
        "",
        f"center_x = {center[0]}",
        f"center_y = {center[1]}",
        f"center_z = {center[2]}",
        "",
        f"size_x = {box_size[0]}",
        f"size_y = {box_size[1]}",
        f"size_z = {box_size[2]}",
        "",
        f"out = {output_path}",
        "",
        f"exhaustiveness = {exhaustiveness}",
        f"num_modes = {num_modes}",
        f"energy_range = {energy_range}",
        f"seed = {seed}",
        f"cpu = {cpu}"
    ]
    
    with open(config_path, 'w') as f:
        f.write("\n".join(config_content) + "\n")
        
    return config_path
