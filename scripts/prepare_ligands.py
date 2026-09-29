"""
Ligand preparation entry point for the drug repurposing pipeline.
"""

import argparse
import json
import logging
import sys
import time
from pathlib import Path
from typing import Dict, Any

import pandas as pd
import yaml

# Add project root to sys.path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.ligand import prepare_ligand_library, prepare_reference_ligand
from src.descriptors import calculate_library_descriptors

def setup_logging() -> None:
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

def main() -> None:
    start_time = time.time()
    setup_logging()
    logger = logging.getLogger(__name__)
    
    parser = argparse.ArgumentParser(description="Prepare ligands for docking.")
    parser.add_argument(
        "--config", 
        type=str, 
        default="config/docking_config.yaml",
        help="Path to docking configuration file."
    )
    args = parser.parse_args()
    
    logger.info("Starting ligand preparation pipeline...")
    
    project_root = Path(__file__).resolve().parent.parent
    config_path = project_root / args.config
    data_dir = project_root / "data"
    results_dir = project_root / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    
    compounds_csv = data_dir / "compounds.csv"
    ligand_pdb_path = data_dir / "processed" / "1X70_ligand_715.pdb"
    
    config = load_config(str(config_path))
    
    report = {
        "status": "started",
        "successful_count": 0,
        "failed_count": 0,
        "failures": {}
    }
    
    try:
        # 3. Prepare ligand library
        logger.info(f"Preparing ligand library from {compounds_csv}")
        compounds_df = pd.read_csv(compounds_csv)
        library_results = prepare_ligand_library(
            compounds_df,
            str(data_dir / "processed" / "ligands")
        )
        
        report["successful_count"] = library_results.get("successful_count", 0)
        report["failed_count"] = library_results.get("failed_count", 0)
        report["failures"] = library_results.get("failures", {})
        
        # 4. Prepare reference ligand
        if ligand_pdb_path.exists():
            logger.info("Preparing reference ligand for redocking")
            ref_pdbqt = prepare_reference_ligand(str(ligand_pdb_path), "715", "A", str(data_dir / "processed"))
            report["reference_ligand_pdbqt"] = ref_pdbqt
        else:
            logger.warning(f"Reference ligand PDB not found at {ligand_pdb_path}, skipping reference preparation.")
        
        # 5 & 6. Calculate and save descriptors
        logger.info("Calculating molecular descriptors")
        descriptors_df = calculate_library_descriptors(compounds_df)
        descriptors_out = results_dir / "compound_descriptors.csv"
        descriptors_df.to_csv(descriptors_out, index=False)
        logger.info(f"Saved descriptors to {descriptors_out}")
        
        report["status"] = "success"
        
    except Exception as e:
        logger.error(f"Ligand preparation failed: {e}")
        import traceback
        traceback.print_exc()
        report["status"] = "failed"
        report["error"] = str(e)
    finally:
        # 7. Save report
        report_path = results_dir / "ligand_preparation_report.json"
        with open(report_path, 'w', encoding='utf-8') as f:
            json.dump({k: str(v) if hasattr(v, '__fspath__') else v for k, v in report.items()}, f, indent=4, default=str)
        logger.info(f"Saved preparation report to {report_path}")
        
        # 8. Print summary
        elapsed_time = time.time() - start_time
        logger.info(f"Ligand preparation finished in {elapsed_time:.2f} seconds.")
        logger.info(f"Summary: {report['successful_count']} successful, {report['failed_count']} failed.")
        if report["failed_count"] > 0:
            logger.info(f"Failure reasons: {report['failures']}")

if __name__ == "__main__":
    main()
