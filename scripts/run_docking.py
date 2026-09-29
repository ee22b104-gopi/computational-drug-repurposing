"""
Docking entry point for the drug repurposing pipeline.
"""

import argparse
import logging
import sys
import time
from pathlib import Path
from typing import Dict, Any
import csv
import pandas as pd

import yaml

# Add project root to sys.path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.docking import check_vina_installation, dock_library
from src.validation import redocking_validation

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
    
    parser = argparse.ArgumentParser(description="Run AutoDock Vina pipeline.")
    parser.add_argument("--config", type=str, default="config/docking_config.yaml", help="Path to config")
    args = parser.parse_args()
    
    logger.info("Starting docking pipeline...")
    
    project_root = Path(__file__).resolve().parent.parent
    config_path = project_root / args.config
    data_dir = project_root / "data"
    results_dir = project_root / "results"
    
    config = load_config(str(config_path))
    
    vina_executable = config.get("docking", {}).get("vina_executable", "bin/vina.exe")
    receptor_pdbqt = data_dir / "processed" / "1X70_receptor.pdbqt"
    ligands_dir = data_dir / "processed" / "ligands"
    compounds_csv = data_dir / "compounds.csv"
    
    try:
        logger.info("Checking Vina installation")
        if not check_vina_installation(vina_executable=vina_executable):
            raise RuntimeError("AutoDock Vina executable not found.")
            
        logger.info("Step 1: Running full library docking")
        compounds_df = pd.read_csv(compounds_csv)
        center = (
            config["binding_site"]["center_x"],
            config["binding_site"]["center_y"],
            config["binding_site"]["center_z"]
        )
        box_size = (
            config["binding_site"]["size_x"],
            config["binding_site"]["size_y"],
            config["binding_site"]["size_z"]
        )
        
        results_df = dock_library(
            receptor_pdbqt=str(receptor_pdbqt),
            ligand_dir=str(ligands_dir),
            output_dir=str(results_dir / "docking"),
            center=center,
            box_size=box_size,
            compounds_df=compounds_df,
            exhaustiveness=config["docking"]["exhaustiveness"],
            num_modes=config["docking"]["num_modes"],
            energy_range=config["docking"]["energy_range"],
            seed=config["docking"]["seed"],
            cpu=config["docking"]["cpu"],
            vina_executable=vina_executable
        )
        
        results_csv = results_dir / "docking_results.csv"
        results_df.to_csv(results_csv, index=False)
        logger.info(f"Docking results saved to {results_csv}")
        
    except Exception as e:
        logger.error(f"Docking pipeline failed: {e}")
        import traceback
        traceback.print_exc()
        
    finally:
        elapsed = time.time() - start_time
        logger.info(f"Docking pipeline finished in {elapsed:.2f} seconds.")

if __name__ == "__main__":
    main()
