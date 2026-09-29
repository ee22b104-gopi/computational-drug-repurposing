#!/usr/bin/env python3
"""
Generate all publication-quality figures.

Usage:
    python scripts/generate_figures.py [--config config/docking_config.yaml]
"""

import argparse
import logging
import sys
import time
from pathlib import Path

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import pandas as pd

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
    parser = argparse.ArgumentParser(description="Generate publication figures")
    parser.add_argument(
        "--config",
        default=str(PROJECT_ROOT / "config" / "docking_config.yaml"),
        help="Path to configuration file",
    )
    args = parser.parse_args()

    start_time = time.time()
    config = load_config(args.config)

    results_dir = PROJECT_ROOT / config["paths"]["results_dir"]
    figures_dir = PROJECT_ROOT / config["paths"]["figures_dir"]
    figures_dir.mkdir(parents=True, exist_ok=True)

    # Import visualization module
    from src.visualization import (
        plot_workflow_diagram,
        plot_binding_pocket,
        plot_score_distribution,
        plot_controls_vs_library,
        plot_top_candidate_scores,
        plot_interaction_heatmap,
        plot_chemical_properties,
        plot_robustness_scores,
        plot_rank_correlation,
    )

    # Identify positive control IDs
    positive_control_ids = [
        pc["name"].lower().replace(" ", "_")
        for pc in config["library"]["positive_controls"]
    ]

    # -------------------------------------------------------------------------
    # Figure 1: Workflow diagram
    # -------------------------------------------------------------------------
    logger.info("Generating Figure 1: Workflow diagram")
    try:
        plot_workflow_diagram(str(figures_dir / "fig1_workflow.png"))
    except Exception as e:
        logger.error(f"Figure 1 failed: {e}")

    # -------------------------------------------------------------------------
    # Figure 2: Binding pocket (schematic)
    # -------------------------------------------------------------------------
    logger.info("Generating Figure 2: Binding pocket")
    try:
        plot_binding_pocket(
            receptor_pdb=None,  # Schematic only
            ligand_pdb=None,
            center=None,
            box_size=None,
            residue_list=None,
            output_path=str(figures_dir / "fig2_binding_pocket.png"),
        )
    except Exception as e:
        logger.error(f"Figure 2 failed: {e}")

    # -------------------------------------------------------------------------
    # Load data for remaining figures
    # -------------------------------------------------------------------------
    docking_path = results_dir / "docking_results.csv"
    descriptors_path = results_dir / "compound_descriptors.csv"
    interactions_path = results_dir / "interactions.csv"
    robustness_path = results_dir / "ranking_robustness.csv"

    docking_df = pd.read_csv(docking_path) if docking_path.exists() else None
    descriptors_df = pd.read_csv(descriptors_path) if descriptors_path.exists() else None
    interactions_df = pd.read_csv(interactions_path) if interactions_path.exists() else None
    robustness_df = pd.read_csv(robustness_path) if robustness_path.exists() else None

    # -------------------------------------------------------------------------
    # Figure 3: Score distribution
    # -------------------------------------------------------------------------
    if docking_df is not None:
        logger.info("Generating Figure 3: Score distribution")
        try:
            plot_score_distribution(
                docking_df, str(figures_dir / "fig3_score_distribution.png")
            )
        except Exception as e:
            logger.error(f"Figure 3 failed: {e}")

        # Figure 4: Controls vs library
        logger.info("Generating Figure 4: Controls vs library")
        try:
            plot_controls_vs_library(
                docking_df,
                positive_control_ids,
                str(figures_dir / "fig4_controls_vs_library.png"),
            )
        except Exception as e:
            logger.error(f"Figure 4 failed: {e}")

        # Figure 5: Top candidate scores
        logger.info("Generating Figure 5: Top candidate scores")
        try:
            plot_top_candidate_scores(
                docking_df,
                top_n=20,
                output_path=str(figures_dir / "fig5_top_candidates.png"),
            )
        except Exception as e:
            logger.error(f"Figure 5 failed: {e}")

    # Figure 6: Interaction heatmap
    if interactions_df is not None:
        logger.info("Generating Figure 6: Interaction heatmap")
        try:
            plot_interaction_heatmap(
                interactions_df, str(figures_dir / "fig6_interactions.png")
            )
        except Exception as e:
            logger.error(f"Figure 6 failed: {e}")

    # Figure 7: Chemical properties
    if descriptors_df is not None:
        logger.info("Generating Figure 7: Chemical properties")
        try:
            plot_chemical_properties(
                descriptors_df, str(figures_dir / "fig7_chemical_properties.png")
            )
        except Exception as e:
            logger.error(f"Figure 7 failed: {e}")

    # Figure 8: Robustness scores
    if robustness_df is not None:
        logger.info("Generating Figure 8: Robustness scores")
        try:
            plot_robustness_scores(
                robustness_df, str(figures_dir / "fig8_robustness.png")
            )
        except Exception as e:
            logger.error(f"Figure 8 failed: {e}")

        # Figure 9: Rank correlation
        logger.info("Generating Figure 9: Rank correlation")
        try:
            rank_corr_path = results_dir / "rank_correlations.csv"
            if rank_corr_path.exists():
                rank_corr_df = pd.read_csv(rank_corr_path)
                plot_rank_correlation(
                    rank_corr_df, str(figures_dir / "fig9_rank_correlation.png")
                )
        except Exception as e:
            logger.error(f"Figure 9 failed: {e}")

    elapsed = time.time() - start_time
    logger.info(f"Figure generation complete in {elapsed:.1f} seconds")
    logger.info(f"Figures saved to {figures_dir}")


if __name__ == "__main__":
    main()
