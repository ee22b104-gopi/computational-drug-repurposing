#!/usr/bin/env python3
"""
Master pipeline script — runs all stages sequentially.

Usage:
    python run_pipeline.py [--config config/docking_config.yaml] [--skip-download]

This script orchestrates the complete computational drug repurposing workflow:
    Stage 1-3:  python scripts/prepare_protein.py
    Stage 4-5:  python data/download_data.py + python scripts/prepare_ligands.py
    Stage 7-9:  python scripts/run_docking.py
    Stage 10-15: python scripts/analyze_results.py
    Stage 13:   python scripts/generate_figures.py
"""

import argparse
import logging
import subprocess
import sys
import time
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger("pipeline")


def run_stage(script_path: str, description: str, config_path: str) -> bool:
    """Run a pipeline stage as a subprocess.

    Args:
        script_path: Path to the Python script.
        description: Human-readable stage description.
        config_path: Path to the config file.

    Returns:
        True if successful, False otherwise.
    """
    logger.info("=" * 70)
    logger.info(f"  {description}")
    logger.info("=" * 70)

    start = time.time()
    cmd = [sys.executable, str(script_path)]
    if "--config" not in str(script_path):
        cmd.extend(["--config", str(config_path)])

    try:
        result = subprocess.run(
            cmd,
            cwd=str(PROJECT_ROOT),
            capture_output=False,
            text=True,
            timeout=7200,  # 2 hour timeout per stage
        )
        elapsed = time.time() - start

        if result.returncode == 0:
            logger.info(f"  ✓ Completed in {elapsed:.1f}s")
            return True
        else:
            logger.error(f"  ✗ Failed (exit code {result.returncode}) after {elapsed:.1f}s")
            return False

    except subprocess.TimeoutExpired:
        logger.error(f"  ✗ Timed out after 7200s")
        return False
    except Exception as e:
        logger.error(f"  ✗ Error: {e}")
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Run the complete computational drug repurposing pipeline"
    )
    parser.add_argument(
        "--config",
        default=str(PROJECT_ROOT / "config" / "docking_config.yaml"),
        help="Path to configuration file",
    )
    parser.add_argument(
        "--skip-download",
        action="store_true",
        help="Skip data download (use existing data)",
    )
    parser.add_argument(
        "--skip-docking",
        action="store_true",
        help="Skip docking (use existing results)",
    )
    args = parser.parse_args()

    pipeline_start = time.time()
    config_path = Path(args.config)

    logger.info("╔" + "═" * 68 + "╗")
    logger.info("║  Structure-Based Virtual Screening for Drug Repurposing Pipeline  ║")
    logger.info("║  Target: DPP-4 (Human Dipeptidyl Peptidase-4)                    ║")
    logger.info("╚" + "═" * 68 + "╝")
    logger.info("")

    stages = []

    # Stage 1: Download data
    if not args.skip_download:
        stages.append((
            PROJECT_ROOT / "data" / "download_data.py",
            "STAGE 1-4: Downloading PDB structures and assembling compound library",
        ))

    # Stage 2: Prepare protein
    stages.append((
        PROJECT_ROOT / "scripts" / "prepare_protein.py",
        "STAGE 2-3,6: Protein structure preparation and binding site identification",
    ))

    # Stage 3: Prepare ligands
    stages.append((
        PROJECT_ROOT / "scripts" / "prepare_ligands.py",
        "STAGE 5: Ligand preprocessing (SMILES → 3D → PDBQT)",
    ))

    # Stage 4: Run docking
    if not args.skip_docking:
        stages.append((
            PROJECT_ROOT / "scripts" / "run_docking.py",
            "STAGE 7-9: Docking pipeline (validation + virtual screening)",
        ))

    # Stage 5: Analyze results
    stages.append((
        PROJECT_ROOT / "scripts" / "analyze_results.py",
        "STAGE 10-15: Chemical filtering, interaction analysis, candidate selection",
    ))

    # Stage 6: Generate figures
    stages.append((
        PROJECT_ROOT / "scripts" / "generate_figures.py",
        "STAGE 13: Generating publication-quality figures",
    ))

    # Execute stages
    results = {}
    for script_path, description in stages:
        if not script_path.exists():
            logger.warning(f"Script not found: {script_path}")
            results[description] = "SKIPPED"
            continue

        success = run_stage(script_path, description, str(config_path))
        results[description] = "SUCCESS" if success else "FAILED"

        if not success:
            logger.warning("Stage failed — continuing with remaining stages")

    # Summary
    total_time = time.time() - pipeline_start
    logger.info("")
    logger.info("╔" + "═" * 68 + "╗")
    logger.info("║  Pipeline Summary                                                ║")
    logger.info("╚" + "═" * 68 + "╝")
    for description, status in results.items():
        icon = "✓" if status == "SUCCESS" else "✗" if status == "FAILED" else "⊘"
        logger.info(f"  {icon} {description}: {status}")
    logger.info("")
    logger.info(f"  Total pipeline time: {total_time / 60:.1f} minutes")
    logger.info("")

    # Check for failures
    failures = [d for d, s in results.items() if s == "FAILED"]
    if failures:
        logger.warning(f"  {len(failures)} stage(s) failed. Check logs above.")
        sys.exit(1)
    else:
        logger.info("  All stages completed successfully!")
        logger.info("")
        logger.info("  Results are in: results/")
        logger.info("  Figures are in: results/figures/")
        logger.info("  Run 'python scripts/generate_figures.py' to regenerate figures")


if __name__ == "__main__":
    main()
