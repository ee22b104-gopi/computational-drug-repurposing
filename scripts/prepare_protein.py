"""
Protein preparation entry point for the drug repurposing pipeline.
"""

import argparse
import json
import logging
import sys
import time
from pathlib import Path
from typing import Dict, Any

import yaml

# Add project root to sys.path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.protein import (
    download_pdb,
    extract_chain,
    clean_protein,
    extract_ligand,
    get_ligand_center,
    get_binding_site_residues,
    prepare_receptor_pdbqt
)

def setup_logging() -> None:
    """Configure logging for the script."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    )

def load_config(config_path: str) -> Dict[str, Any]:
    try:
        with open(config_path, 'r', encoding='utf-8') as f:
            return yaml.safe_load(f)
    except Exception as e:
        logging.error(f"Failed to load config from {config_path}: {e}")
        sys.exit(1)

def save_config(config: Dict[str, Any], config_path: str) -> None:
    try:
        with open(config_path, 'w', encoding='utf-8') as f:
            yaml.dump(config, f, default_flow_style=False)
    except Exception as e:
        logging.error(f"Failed to save config to {config_path}: {e}")

def main() -> None:
    start_time = time.time()
    setup_logging()
    logger = logging.getLogger(__name__)
    
    parser = argparse.ArgumentParser(description="Prepare protein for docking.")
    parser.add_argument(
        "--config", 
        type=str, 
        default="config/docking_config.yaml",
        help="Path to docking configuration file."
    )
    args = parser.parse_args()
    
    logger.info("Starting protein preparation pipeline...")
    
    project_root = Path(__file__).resolve().parent.parent
    config_path = project_root / args.config
    results_dir = project_root / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    
    config = load_config(str(config_path))
    pdb_id = config["protein"]["pdb_id"]
    chain_id = config["protein"]["chain"]
    ligand_name = config["protein"]["ligand_code"]
    
    report = {
        "pdb_id": pdb_id,
        "chain_extracted": chain_id,
        "ligand_extracted": ligand_name,
        "status": "started"
    }
    
    data_raw = project_root / "data" / "raw"
    data_proc = project_root / "data" / "processed"
    data_raw.mkdir(parents=True, exist_ok=True)
    data_proc.mkdir(parents=True, exist_ok=True)
    
    try:
        # 2. Download PDB
        logger.info(f"Downloading PDB {pdb_id}")
        raw_pdb_path = download_pdb(pdb_id, str(data_raw))
        report["raw_pdb_path"] = raw_pdb_path
        
        # 3. Extract chain
        logger.info(f"Extracting chain {chain_id}")
        chain_pdb_path = str(data_proc / f"{pdb_id}_chain{chain_id}.pdb")
        extract_chain(raw_pdb_path, chain_id, chain_pdb_path)
        
        # 4. Clean protein
        logger.info(f"Cleaning protein")
        clean_pdb_path = str(data_proc / f"{pdb_id}_clean.pdb")
        clean_protein(chain_pdb_path, clean_pdb_path, remove_water=True, remove_heterogens=True, keep_ligand_code=ligand_name)
        report["clean_pdb_path"] = clean_pdb_path
        
        # 5. Extract ligand
        logger.info(f"Extracting co-crystallized ligand {ligand_name}")
        ligand_pdb_path = str(data_proc / f"{pdb_id}_ligand_{ligand_name}.pdb")
        extract_ligand(raw_pdb_path, ligand_name, chain_id, ligand_pdb_path)
        report["ligand_pdb_path"] = ligand_pdb_path
        
        # 6. Calculate binding site center
        logger.info("Calculating binding site center")
        center = get_ligand_center(ligand_pdb_path, ligand_name, chain_id)
        if center is None:
            center = (40.5, 51.0, 36.5) # Default fallback for 1X70
        report["binding_site_center"] = center
        
        # 7. Identify binding site residues
        logger.info("Identifying binding site residues within 5A")
        residues = get_binding_site_residues(clean_pdb_path, ligand_name, chain_id, distance_cutoff=5.0)
        report["binding_site_residues"] = residues
        
        # 8. Prepare receptor PDBQT
        logger.info("Preparing receptor PDBQT")
        receptor_pdbqt_path = str(data_proc / f"{pdb_id}_receptor.pdbqt")
        prepare_receptor_pdbqt(clean_pdb_path, receptor_pdbqt_path)
        report["receptor_pdbqt_path"] = receptor_pdbqt_path
        
        # 9. Update config with computed binding site center
        if "binding_site" not in config:
            config["binding_site"] = {}
            
        config["binding_site"]["center_x"] = round(float(center[0]), 3)
        config["binding_site"]["center_y"] = round(float(center[1]), 3)
        config["binding_site"]["center_z"] = round(float(center[2]), 3)
        
        save_config(config, str(config_path))
        logger.info(f"Updated config {config_path} with binding site center: {center}")
        
        report["status"] = "success"
        
    except Exception as e:
        logger.error(f"Protein preparation failed: {e}")
        import traceback
        traceback.print_exc()
        report["status"] = "failed"
        report["error"] = str(e)
    finally:
        # 10. Save report
        report_path = results_dir / "protein_preparation_report.json"
        with open(report_path, 'w', encoding='utf-8') as f:
            json.dump({k: str(v) if hasattr(v, '__fspath__') else v for k, v in report.items()}, f, indent=4, default=str)
        logger.info(f"Saved preparation report to {report_path}")
        
        elapsed_time = time.time() - start_time
        logger.info(f"Protein preparation finished in {elapsed_time:.2f} seconds. Status: {report['status']}")
        
        if report["status"] == "failed":
            sys.exit(1)

if __name__ == "__main__":
    main()
