"""
Ligand preprocessing module for computational drug discovery pipeline.
Handles SMILES conversion, 3D structure generation, and PDBQT preparation.
"""

import json
import logging
from pathlib import Path
from typing import Dict, List, Optional, Union, Any
import multiprocessing
from functools import partial

import pandas as pd
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import MolStandardize

try:
    from meeko import MoleculePreparation
except ImportError:
    from meeko import PDBQTMolecule

logger = logging.getLogger(__name__)

def smiles_to_mol(smiles: str, name: Optional[str] = None) -> Optional[Chem.Mol]:
    """
    Converts a SMILES string to an RDKit Mol object.

    Args:
        smiles: SMILES string of the molecule.
        name: Optional name for the molecule.

    Returns:
        RDKit Mol object or None if invalid.
    """
    if not smiles or not isinstance(smiles, str):
        logger.warning(f"Invalid SMILES provided for {name}: {smiles}")
        return None

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.warning(f"Failed to parse SMILES for {name}: {smiles}")
            return None
        
        if name:
            mol.SetProp("_Name", name)
            
        return mol
    except Exception as e:
        logger.warning(f"Error parsing SMILES for {name}: {str(e)}")
        return None

def standardize_mol(mol: Chem.Mol) -> Optional[Chem.Mol]:
    if mol is None: return None
    try:
        from rdkit.Chem.MolStandardize import rdMolStandardize
        mol = Chem.RemoveHs(mol)
        mol = rdMolStandardize.MetalDisconnector().Disconnect(mol)
        mol = rdMolStandardize.LargestFragmentChooser().choose(mol)
        mol = rdMolStandardize.Uncharger().uncharge(mol)
        return Chem.AddHs(mol)
    except Exception as e:
        return Chem.AddHs(mol)

def generate_3d_coordinates(mol: Chem.Mol, n_conformers: int = 1, seed: int = 42) -> Optional[Chem.Mol]:
    """
    Generates 3D coordinates using RDKit's AllChem.EmbedMolecule with ETKDGv3().
    Adds hydrogens first. Optimizes with MMFF94 force field.

    Args:
        mol: RDKit Mol object.
        n_conformers: Number of conformers to generate.
        seed: Random seed for reproducibility.

    Returns:
        RDKit Mol object with 3D coordinates or None on failure.
    """
    if mol is None:
        return None

    try:
        # Add hydrogens
        mol = Chem.AddHs(mol)

        # Generate 3D coordinates
        params = AllChem.ETKDGv3()
        params.randomSeed = seed
        
        embed_res = AllChem.EmbedMolecule(mol, params)
        if embed_res == -1:
            name = mol.GetProp("_Name") if mol.HasProp("_Name") else "Unknown"
            logger.warning(f"Failed to embed molecule {name}")
            return None

        # Optimize with MMFF94
        try:
            if AllChem.MMFFHasAllMoleculeParams(mol):
                AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
            else:
                AllChem.UFFOptimizeMolecule(mol, maxIters=500)
        except Exception as e:
            name = mol.GetProp("_Name") if mol.HasProp("_Name") else "Unknown"
            logger.warning(f"Force field optimization failed for {name}: {str(e)}")

        return mol
    except Exception as e:
        name = mol.GetProp("_Name") if mol.HasProp("_Name") else "Unknown"
        logger.warning(f"3D coordinate generation failed for {name}: {str(e)}")
        return None

def mol_to_pdbqt(mol: Chem.Mol, output_path: Union[str, Path], name: Optional[str] = None) -> Optional[str]:
    """
    Converts RDKit Mol to PDBQT format using Meeko.

    Args:
        mol: RDKit Mol object (with 3D coordinates).
        output_path: Path to save the PDBQT file.
        name: Optional molecule name.

    Returns:
        Path to the saved PDBQT file as a string, or None on failure.
    """
    if mol is None:
        return None

    try:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Meeko preparation
        preparator = MoleculePreparation()
        preparator.prepare(mol)
        
        pdbqt_string = preparator.write_pdbqt_string()
        
        with open(output_path, 'w') as f:
            f.write(pdbqt_string)
            
        return str(output_path)
    except Exception as e:
        mol_name = name or (mol.GetProp("_Name") if mol.HasProp("_Name") else "Unknown")
        logger.warning(f"PDBQT conversion failed for {mol_name}: {str(e)}")
        return None

def prepare_ligand_from_smiles(smiles: str, name: str, output_dir: Union[str, Path]) -> Dict[str, Any]:
    """
    Full pipeline: SMILES -> validate -> standardize -> 3D -> PDBQT.

    Args:
        smiles: SMILES string.
        name: Name of the molecule.
        output_dir: Directory to save the output PDBQT.

    Returns:
        Dict with status, path, and errors.
    """
    result = {
        'status': 'failed',
        'pdbqt_path': None,
        'error_message': None
    }
    
    output_dir = Path(output_dir)
    output_path = output_dir / f"{name}.pdbqt"
    
    # 1. Parse
    mol = smiles_to_mol(smiles, name)
    if mol is None:
        result['error_message'] = "Invalid SMILES"
        return result
        
    # 2. Standardize
    mol = standardize_mol(mol)
    if mol is None:
        result['error_message'] = "Standardization failed"
        return result
        
    # 3. 3D Coordinates
    mol_3d = generate_3d_coordinates(mol)
    if mol_3d is None:
        result['error_message'] = "3D coordinate generation failed"
        return result
        
    # 4. PDBQT
    pdbqt_path = mol_to_pdbqt(mol_3d, output_path, name)
    if pdbqt_path is None:
        result['error_message'] = "PDBQT conversion failed"
        return result
        
    result['status'] = 'success'
    result['pdbqt_path'] = pdbqt_path
    return result

def _prepare_worker(row, output_dir: Path):
    """Worker function for multiprocessing"""
    compound_id = row['compound_id']
    name = row['name']
    smiles = row['smiles']
    
    safe_name = str(name).replace(" ", "_").replace("/", "_").replace("\\", "_")
    if not safe_name:
        safe_name = str(compound_id)
        
    res = prepare_ligand_from_smiles(smiles, safe_name, output_dir)
    
    return {
        'compound_id': compound_id,
        'name': name,
        'smiles': smiles,
        'pdbqt_path': res['pdbqt_path'],
        'status': res['status'],
        'error_message': res['error_message']
    }

def prepare_ligand_library(compounds_df: pd.DataFrame, output_dir: Union[str, Path], n_workers: int = 1) -> pd.DataFrame:
    """
    Processes a DataFrame of compounds.

    Args:
        compounds_df: DataFrame with columns 'name', 'smiles', 'compound_id'.
        output_dir: Output directory for PDBQT files.
        n_workers: Number of parallel workers.

    Returns:
        DataFrame with results.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    required_cols = {'name', 'smiles', 'compound_id'}
    if not required_cols.issubset(compounds_df.columns):
        missing = required_cols - set(compounds_df.columns)
        raise ValueError(f"Missing required columns in DataFrame: {missing}")
        
    rows = [row for _, row in compounds_df.iterrows()]
    worker_fn = partial(_prepare_worker, output_dir=output_dir)
    
    results = []
    
    if n_workers > 1:
        with multiprocessing.Pool(processes=n_workers) as pool:
            results = list(tqdm(
                pool.imap(worker_fn, rows),
                total=len(rows),
                desc="Preparing ligands"
            ))
    else:
        for row in tqdm(rows, desc="Preparing ligands"):
            results.append(worker_fn(row))
            
    return pd.DataFrame(results)

def extract_reference_ligand_smiles(pdb_path: Union[str, Path], ligand_code: str) -> Optional[str]:
    """
    Extracts the SMILES of a co-crystallized ligand from a PDB file.

    Args:
        pdb_path: Path to the PDB file.
        ligand_code: Three-letter code of the ligand.

    Returns:
        SMILES string or None if extraction fails.
    """
    try:
        mol = Chem.MolFromPDBFile(str(pdb_path), sanitize=False)
        if mol is None:
            logger.warning(f"Could not read PDB file {pdb_path}")
            return None
            
        ligand_atoms = []
        for atom in mol.GetAtoms():
            info = atom.GetPDBResidueInfo()
            if info and info.GetResidueName().strip() == ligand_code:
                ligand_atoms.append(atom.GetIdx())
                
        if not ligand_atoms:
            logger.warning(f"Ligand {ligand_code} not found in PDB file {pdb_path}")
            return None
            
        rwmol = Chem.RWMol(mol)
        atoms_to_remove = [i for i in range(rwmol.GetNumAtoms()) if i not in ligand_atoms]
        for i in sorted(atoms_to_remove, reverse=True):
            rwmol.RemoveAtom(i)
            
        ligand_mol = rwmol.GetMol()
        Chem.SanitizeMol(ligand_mol)
        
        return Chem.MolToSmiles(ligand_mol)
    except Exception as e:
        logger.warning(f"Error extracting reference ligand SMILES: {str(e)}")
        return None

def prepare_reference_ligand(pdb_path: Union[str, Path], ligand_code: str, chain_id: str, output_dir: Union[str, Path]) -> Dict[str, Any]:
    """
    Extracts co-crystallized ligand from PDB, generates PDBQT for redocking.
    Uses the 3D coordinates directly from the crystal structure.

    Args:
        pdb_path: Path to PDB file.
        ligand_code: Three-letter code of the ligand.
        chain_id: Chain ID where the ligand is bound.
        output_dir: Directory to save the output PDBQT.

    Returns:
        Dict with status, path, and errors.
    """
    result = {
        'status': 'failed',
        'pdbqt_path': None,
        'error_message': None
    }
    
    try:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        mol = Chem.MolFromPDBFile(str(pdb_path), sanitize=False)
        if mol is None:
            result['error_message'] = f"Could not read PDB file {pdb_path}"
            return result
            
        ligand_atoms = []
        for atom in mol.GetAtoms():
            info = atom.GetPDBResidueInfo()
            if info and info.GetResidueName().strip() == ligand_code and info.GetChainId() == chain_id:
                ligand_atoms.append(atom.GetIdx())
                
        if not ligand_atoms:
            result['error_message'] = f"Ligand {ligand_code} on chain {chain_id} not found"
            return result
            
        rwmol = Chem.RWMol(mol)
        atoms_to_remove = [i for i in range(rwmol.GetNumAtoms()) if i not in ligand_atoms]
        for i in sorted(atoms_to_remove, reverse=True):
            rwmol.RemoveAtom(i)
            
        ligand_mol = rwmol.GetMol()
        
        try:
            Chem.SanitizeMol(ligand_mol, catchErrors=True)
            ligand_mol = Chem.AddHs(ligand_mol, addCoords=True)
        except Exception as e:
            logger.warning(f"Failed to add hydrogens cleanly, proceeding with what we have: {e}")
            pass
            
        output_path = output_dir / f"ref_{ligand_code}_{chain_id}.pdbqt"
        
        pdbqt_path = mol_to_pdbqt(ligand_mol, output_path, f"ref_{ligand_code}")
        if pdbqt_path is None:
            result['error_message'] = "PDBQT conversion failed for reference ligand"
            return result
            
        result['status'] = 'success'
        result['pdbqt_path'] = pdbqt_path
        
    except Exception as e:
        result['error_message'] = str(e)
        
    return result

def generate_preprocessing_report(results_df: pd.DataFrame, output_path: Union[str, Path]) -> None:
    """
    Creates a JSON preprocessing report.

    Args:
        results_df: Results DataFrame from prepare_ligand_library.
        output_path: Path to save JSON report.
    """
    output_path = Path(output_path)
    
    total = len(results_df)
    successful = len(results_df[results_df['status'] == 'success'])
    failed = total - successful
    
    failure_reasons = {}
    for error in results_df[results_df['status'] == 'failed']['error_message']:
        if error:
            failure_reasons[error] = failure_reasons.get(error, 0) + 1
            
    report = {
        'summary': {
            'total_compounds': total,
            'successful': successful,
            'failed': failed,
            'success_rate': round((successful / total * 100) if total > 0 else 0, 2)
        },
        'failure_breakdown': failure_reasons
    }
    
    with open(output_path, 'w') as f:
        json.dump(report, f, indent=4)
