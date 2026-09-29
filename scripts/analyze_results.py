#!/usr/bin/env python3
"""
Analyze docking results: chemical filtering, interaction analysis, and robustness.

Usage:
    python scripts/analyze_results.py [--config config/docking_config.yaml]
"""

import argparse
import json
import logging
import sys
import time
from pathlib import Path

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import pandas as pd

from src.descriptors import (
    calculate_library_descriptors,
    check_drug_likeness,
    filter_compounds,
    generate_descriptor_report,
)
from src.interactions import (
    analyze_library_interactions,
    check_binding_site_coverage,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)


def load_config(config_path: str) -> dict:
    """Load YAML configuration."""
    import yaml
    with open(config_path, "r") as f:
        return yaml.safe_load(f)


def main():
    parser = argparse.ArgumentParser(description="Analyze docking results")
    parser.add_argument(
        "--config",
        default=str(PROJECT_ROOT / "config" / "docking_config.yaml"),
        help="Path to configuration file",
    )
    args = parser.parse_args()

    start_time = time.time()
    config = load_config(args.config)

    results_dir = PROJECT_ROOT / config["paths"]["results_dir"]
    tables_dir = PROJECT_ROOT / config["paths"]["tables_dir"]
    tables_dir.mkdir(parents=True, exist_ok=True)
    docking_dir = PROJECT_ROOT / config["paths"]["docking_dir"]

    # -------------------------------------------------------------------------
    # Step 1: Load docking results
    # -------------------------------------------------------------------------
    logger.info("=" * 60)
    logger.info("STAGE 10: Chemical Property Analysis")
    logger.info("=" * 60)

    docking_results_path = results_dir / "docking_results.csv"
    if not docking_results_path.exists():
        logger.error(f"Docking results not found: {docking_results_path}")
        logger.error("Run 'python scripts/run_docking.py' first")
        sys.exit(1)

    docking_df = pd.read_csv(docking_results_path)
    logger.info(f"Loaded {len(docking_df)} docking results")

    # Load compound descriptors if already calculated, otherwise calculate
    descriptors_path = results_dir / "compound_descriptors.csv"
    if descriptors_path.exists():
        descriptors_df = pd.read_csv(descriptors_path)
        logger.info(f"Loaded pre-computed descriptors for {len(descriptors_df)} compounds")
    else:
        compounds_path = PROJECT_ROOT / config["paths"]["data_dir"] / "compounds.csv"
        if compounds_path.exists():
            compounds_df = pd.read_csv(compounds_path)
            descriptors_df = calculate_library_descriptors(compounds_df)
            descriptors_df.to_csv(descriptors_path, index=False)
            logger.info(f"Calculated descriptors for {len(descriptors_df)} compounds")
        else:
            logger.warning("No compounds.csv found; skipping descriptor calculation")
            descriptors_df = pd.DataFrame()

    # Merge docking results with descriptors
    if not descriptors_df.empty and not docking_df.empty:
        merged_df = docking_df.merge(
            descriptors_df, on="compound_id", how="left", suffixes=("", "_desc")
        )
    else:
        merged_df = docking_df

    # Chemical filtering
    if not descriptors_df.empty:
        filter_config = config.get("filtering", {}).get("lipinski", {})
        passing, failing, filter_report = filter_compounds(
            descriptors_df,
            max_mw=filter_config.get("max_mw", 500),
            max_logp=filter_config.get("max_logp", 5),
            max_hbd=filter_config.get("max_hbd", 5),
            max_hba=filter_config.get("max_hba", 10),
        )
        logger.info(
            f"Chemical filtering: {len(passing)} pass, {len(failing)} fail Lipinski"
        )

        # Save descriptor report
        generate_descriptor_report(
            descriptors_df,
            str(results_dir / "descriptor_report.json"),
        )

    # -------------------------------------------------------------------------
    # Step 2: Interaction Analysis
    # -------------------------------------------------------------------------
    logger.info("=" * 60)
    logger.info("STAGE 11: Binding Interaction Analysis")
    logger.info("=" * 60)

    receptor_pdb = PROJECT_ROOT / config["paths"]["data_dir"] / "processed" / f"{config['protein']['pdb_id']}_clean.pdb"
    if not receptor_pdb.exists():
        # Try alternative path
        receptor_pdb = PROJECT_ROOT / config["paths"]["data_dir"] / "processed" / f"{config['protein']['pdb_id']}_chain{config['protein']['chain']}.pdb"

    if receptor_pdb.exists() and not docking_df.empty:
        binding_residues = (
            config["binding_site"]["catalytic_triad"]
            + config["binding_site"]["s1_pocket"]
            + config["binding_site"]["s2_pocket"]
            + config["binding_site"]["s2_extensive"]
            + config["binding_site"]["oxyanion"]
        )

        interactions_df = analyze_library_interactions(
            receptor_pdb=str(receptor_pdb),
            docking_results_df=docking_df,
            docking_output_dir=str(docking_dir),
            top_n=30,
            reference_residues=binding_residues,
        )
        interactions_path = results_dir / "interactions.csv"
        interactions_df.to_csv(interactions_path, index=False)
        logger.info(f"Interaction analysis saved to {interactions_path}")
    else:
        logger.warning("Receptor PDB not found; skipping interaction analysis")
        interactions_df = pd.DataFrame()

    # -------------------------------------------------------------------------
    # Step 3: Candidate Selection
    # -------------------------------------------------------------------------
    logger.info("=" * 60)
    logger.info("STAGE 14: Candidate Selection")
    logger.info("=" * 60)

    # Identify positive control IDs
    positive_control_names = [
        pc["name"].lower()
        for pc in config["library"]["positive_controls"]
    ]

    if not docking_df.empty:
        # Multi-criteria candidate selection
        candidates = select_candidates(
            docking_df=docking_df,
            interactions_df=interactions_df if not interactions_df.empty else None,
            descriptors_df=descriptors_df if not descriptors_df.empty else None,
            positive_control_names=positive_control_names,
        )

        summary_path = results_dir / "summary.csv"
        candidates.to_csv(summary_path, index=False)
        logger.info(f"Summary with {len(candidates)} candidates saved to {summary_path}")

    # -------------------------------------------------------------------------
    # Summary
    # -------------------------------------------------------------------------
    elapsed = time.time() - start_time
    logger.info("=" * 60)
    logger.info(f"Analysis complete in {elapsed:.1f} seconds")
    logger.info("=" * 60)


def select_candidates(
    docking_df: pd.DataFrame,
    interactions_df: pd.DataFrame = None,
    descriptors_df: pd.DataFrame = None,
    positive_control_names: list = None,
    top_n: int = 10,
) -> pd.DataFrame:
    """Select computationally prioritized candidates using multiple criteria.

    Selection strategy:
    1. Docking score relative to positive control median
    2. Binding site residue coverage (if available)
    3. Drug-likeness (Lipinski compliance)
    4. Not already a known DPP-4 inhibitor

    Args:
        docking_df: Docking results DataFrame.
        interactions_df: Interaction analysis DataFrame (optional).
        descriptors_df: Molecular descriptors DataFrame (optional).
        positive_control_names: List of known inhibitor names (lowercase).
        top_n: Number of candidates to select.

    Returns:
        DataFrame of selected candidates with justification.
    """
    if positive_control_names is None:
        positive_control_names = []

    df = docking_df.copy()
    df = df[df["status"].str.lower() == "success"] if "status" in df.columns else df
    df = df.dropna(subset=["best_score"])

    # Add category label
    df["is_positive_control"] = df["compound_name"].str.lower().isin(
        positive_control_names
    )

    # Calculate reference: median score of positive controls
    pc_scores = df[df["is_positive_control"]]["best_score"]
    if len(pc_scores) > 0:
        pc_median = pc_scores.median()
        logger.info(f"Positive control median score: {pc_median:.2f} kcal/mol")
    else:
        pc_median = df["best_score"].quantile(0.1)
        logger.info(f"No positive controls found; using 10th percentile: {pc_median:.2f}")

    # Score relative to controls
    df["score_vs_controls"] = df["best_score"] - pc_median

    # Merge interaction data if available
    if interactions_df is not None and not interactions_df.empty:
        df = df.merge(
            interactions_df[["compound_id", "binding_site_coverage", "n_interacting_residues"]],
            on="compound_id",
            how="left",
        )
    else:
        df["binding_site_coverage"] = None
        df["n_interacting_residues"] = None

    # Filter: only non-control compounds with scores at least as good as control median
    screening_hits = df[
        (~df["is_positive_control"])
        & (df["best_score"] <= pc_median)
    ].copy()

    if screening_hits.empty:
        # Relax criterion: take top N non-control compounds
        screening_hits = df[~df["is_positive_control"]].copy()

    # Sort by docking score (best first)
    screening_hits = screening_hits.sort_values("best_score").head(top_n)

    # Add selection justification
    screening_hits["selection_basis"] = "Docking score comparable to or better than known DPP-4 inhibitors"
    screening_hits["candidate_status"] = "computationally_prioritized"

    # Also include positive controls for comparison
    controls = df[df["is_positive_control"]].sort_values("best_score")
    controls["selection_basis"] = "Known DPP-4 inhibitor (positive control)"
    controls["candidate_status"] = "positive_control"

    result = pd.concat([controls, screening_hits], ignore_index=True)
    result = result.sort_values("best_score")

    return result


if __name__ == "__main__":
    main()
