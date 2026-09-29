"""
Robustness analysis module.

Assesses the sensitivity of docking rankings to reasonable changes
in docking parameters: box dimensions, exhaustiveness, and random seed.
"""

import logging
import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any

import numpy as np
import pandas as pd
from scipy import stats

logger = logging.getLogger(__name__)


def run_robustness_docking(
    receptor_pdbqt: str,
    ligand_pdbqt_dir: str,
    output_base_dir: str,
    compound_ids: List[str],
    center: Tuple[float, float, float],
    default_box_size: float = 22.0,
    box_sizes: List[float] = None,
    exhaustiveness_values: List[int] = None,
    random_seeds: List[int] = None,
    vina_executable: str = "vina",
    cpu: int = 4,
) -> pd.DataFrame:
    """Run docking under multiple parameter configurations for robustness analysis.

    Args:
        receptor_pdbqt: Path to receptor PDBQT file.
        ligand_pdbqt_dir: Directory containing prepared ligand PDBQT files.
        output_base_dir: Base directory for robustness output.
        compound_ids: List of compound IDs to test.
        center: Docking box center (x, y, z).
        default_box_size: Default box size (used when varying other params).
        box_sizes: Box sizes to test (Å per side).
        exhaustiveness_values: Exhaustiveness values to test.
        random_seeds: Random seeds to test.
        vina_executable: Path to Vina executable.
        cpu: Number of CPU threads.

    Returns:
        DataFrame with scores for each compound under each condition.
    """
    from src.docking import dock_single

    if box_sizes is None:
        box_sizes = [18.0, 22.0, 26.0, 30.0]
    if exhaustiveness_values is None:
        exhaustiveness_values = [8, 16, 32, 64]
    if random_seeds is None:
        random_seeds = [42, 137, 256, 314, 500]

    output_base = Path(output_base_dir)
    output_base.mkdir(parents=True, exist_ok=True)
    all_results = []

    # Test 1: Box size variation (fixed exhaustiveness=32, seed=42)
    logger.info("=== Robustness Test 1: Box Size Variation ===")
    for box_size in box_sizes:
        condition = f"box_{box_size}"
        condition_dir = output_base / condition
        condition_dir.mkdir(parents=True, exist_ok=True)

        for cid in compound_ids:
            ligand_path = Path(ligand_pdbqt_dir) / f"{cid}.pdbqt"
            if not ligand_path.exists():
                all_results.append({
                    "compound_id": cid,
                    "condition_type": "box_size",
                    "condition_value": box_size,
                    "condition_label": condition,
                    "best_score": None,
                    "status": "ligand_not_found",
                })
                continue

            output_path = condition_dir / f"{cid}_out.pdbqt"
            result = dock_single(
                receptor_pdbqt=receptor_pdbqt,
                ligand_pdbqt=str(ligand_path),
                center=center,
                box_size=(box_size, box_size, box_size),
                output_path=str(output_path),
                exhaustiveness=32,
                seed=42,
                cpu=cpu,
                vina_executable=vina_executable,
            )
            all_results.append({
                "compound_id": cid,
                "condition_type": "box_size",
                "condition_value": box_size,
                "condition_label": condition,
                "best_score": result.get("best_score"),
                "status": result.get("status", "unknown"),
            })

    # Test 2: Exhaustiveness variation (fixed box=22, seed=42)
    logger.info("=== Robustness Test 2: Exhaustiveness Variation ===")
    for exh in exhaustiveness_values:
        condition = f"exh_{exh}"
        condition_dir = output_base / condition
        condition_dir.mkdir(parents=True, exist_ok=True)

        for cid in compound_ids:
            ligand_path = Path(ligand_pdbqt_dir) / f"{cid}.pdbqt"
            if not ligand_path.exists():
                all_results.append({
                    "compound_id": cid,
                    "condition_type": "exhaustiveness",
                    "condition_value": exh,
                    "condition_label": condition,
                    "best_score": None,
                    "status": "ligand_not_found",
                })
                continue

            output_path = condition_dir / f"{cid}_out.pdbqt"
            result = dock_single(
                receptor_pdbqt=receptor_pdbqt,
                ligand_pdbqt=str(ligand_path),
                center=center,
                box_size=(default_box_size, default_box_size, default_box_size),
                output_path=str(output_path),
                exhaustiveness=exh,
                seed=42,
                cpu=cpu,
                vina_executable=vina_executable,
            )
            all_results.append({
                "compound_id": cid,
                "condition_type": "exhaustiveness",
                "condition_value": exh,
                "condition_label": condition,
                "best_score": result.get("best_score"),
                "status": result.get("status", "unknown"),
            })

    # Test 3: Random seed variation (fixed box=22, exh=32)
    logger.info("=== Robustness Test 3: Random Seed Variation ===")
    for seed in random_seeds:
        condition = f"seed_{seed}"
        condition_dir = output_base / condition
        condition_dir.mkdir(parents=True, exist_ok=True)

        for cid in compound_ids:
            ligand_path = Path(ligand_pdbqt_dir) / f"{cid}.pdbqt"
            if not ligand_path.exists():
                all_results.append({
                    "compound_id": cid,
                    "condition_type": "random_seed",
                    "condition_value": seed,
                    "condition_label": condition,
                    "best_score": None,
                    "status": "ligand_not_found",
                })
                continue

            output_path = condition_dir / f"{cid}_out.pdbqt"
            result = dock_single(
                receptor_pdbqt=receptor_pdbqt,
                ligand_pdbqt=str(ligand_path),
                center=center,
                box_size=(default_box_size, default_box_size, default_box_size),
                output_path=str(output_path),
                exhaustiveness=32,
                seed=seed,
                cpu=cpu,
                vina_executable=vina_executable,
            )
            all_results.append({
                "compound_id": cid,
                "condition_type": "random_seed",
                "condition_value": seed,
                "condition_label": condition,
                "best_score": result.get("best_score"),
                "status": result.get("status", "unknown"),
            })

    return pd.DataFrame(all_results)


def calculate_rank_correlation(
    robustness_df: pd.DataFrame,
    reference_condition: str = "box_22.0",
) -> pd.DataFrame:
    """Calculate Spearman rank correlation between conditions.

    Args:
        robustness_df: Output from run_robustness_docking.
        reference_condition: Condition label to use as reference.

    Returns:
        DataFrame with pairwise rank correlations.
    """
    # Pivot to get scores per condition
    pivot = robustness_df.pivot_table(
        index="compound_id",
        columns="condition_label",
        values="best_score",
        aggfunc="first",
    )

    # Drop compounds with any missing scores
    pivot = pivot.dropna()

    if pivot.empty:
        logger.warning("No compounds with complete data across all conditions")
        return pd.DataFrame()

    conditions = pivot.columns.tolist()
    correlations = []

    for i, cond1 in enumerate(conditions):
        for j, cond2 in enumerate(conditions):
            if i <= j:
                rho, pval = stats.spearmanr(pivot[cond1], pivot[cond2])
                correlations.append({
                    "condition_1": cond1,
                    "condition_2": cond2,
                    "spearman_rho": round(rho, 4),
                    "p_value": round(pval, 6),
                    "n_compounds": len(pivot),
                })

    return pd.DataFrame(correlations)


def analyze_rank_stability(
    robustness_df: pd.DataFrame,
    top_n: int = 10,
) -> Dict[str, Any]:
    """Analyze whether top-ranked compounds remain stable across conditions.

    Args:
        robustness_df: Output from run_robustness_docking.
        top_n: Number of top compounds to track.

    Returns:
        Dictionary with stability metrics.
    """
    conditions = robustness_df["condition_label"].unique()
    top_sets = {}

    for condition in conditions:
        cond_df = robustness_df[
            (robustness_df["condition_label"] == condition)
            & (robustness_df["best_score"].notna())
        ].sort_values("best_score")

        top_ids = set(cond_df.head(top_n)["compound_id"].tolist())
        top_sets[condition] = top_ids

    # Calculate overlap between all pairs
    n_conditions = len(conditions)
    overlaps = []
    for i in range(n_conditions):
        for j in range(i + 1, n_conditions):
            c1, c2 = conditions[i], conditions[j]
            overlap = len(top_sets[c1].intersection(top_sets[c2]))
            jaccard = (
                overlap / len(top_sets[c1].union(top_sets[c2]))
                if top_sets[c1] or top_sets[c2]
                else 0.0
            )
            overlaps.append({
                "condition_1": c1,
                "condition_2": c2,
                "overlap_count": overlap,
                "jaccard_index": round(jaccard, 3),
            })

    # Find consistently top-ranked compounds
    all_top = set()
    for ids in top_sets.values():
        all_top.update(ids)

    consistency = {}
    for cid in all_top:
        n_appearances = sum(1 for ids in top_sets.values() if cid in ids)
        consistency[cid] = {
            "appearances_in_top_n": n_appearances,
            "fraction_of_conditions": round(n_appearances / n_conditions, 3),
        }

    # Stable candidates: appear in top_n in >= 60% of conditions
    stable_threshold = 0.6
    stable_candidates = [
        cid for cid, info in consistency.items()
        if info["fraction_of_conditions"] >= stable_threshold
    ]

    return {
        "top_n": top_n,
        "n_conditions": n_conditions,
        "condition_labels": list(conditions),
        "pairwise_overlaps": overlaps,
        "mean_jaccard": round(
            np.mean([o["jaccard_index"] for o in overlaps]), 3
        ) if overlaps else 0.0,
        "compound_consistency": consistency,
        "stable_candidates": sorted(stable_candidates),
        "n_stable": len(stable_candidates),
        "stability_threshold": stable_threshold,
        "assessment": (
            "Rankings are stable"
            if len(stable_candidates) >= top_n * 0.5
            else "Rankings show moderate instability — interpret top compounds cautiously"
        ),
    }


def score_variation_analysis(
    robustness_df: pd.DataFrame,
) -> pd.DataFrame:
    """Calculate score variation statistics for each compound across conditions.

    Args:
        robustness_df: Output from run_robustness_docking.

    Returns:
        DataFrame with variation statistics per compound.
    """
    valid = robustness_df[robustness_df["best_score"].notna()]

    stats_list = []
    for cid, group in valid.groupby("compound_id"):
        scores = group["best_score"].values
        stats_list.append({
            "compound_id": cid,
            "mean_score": round(np.mean(scores), 3),
            "std_score": round(np.std(scores), 3),
            "min_score": round(np.min(scores), 3),
            "max_score": round(np.max(scores), 3),
            "range_score": round(np.max(scores) - np.min(scores), 3),
            "cv_score": round(np.std(scores) / abs(np.mean(scores)), 4)
            if np.mean(scores) != 0 else 0.0,
            "n_conditions": len(scores),
        })

    return pd.DataFrame(stats_list).sort_values("mean_score")


def generate_robustness_report(
    robustness_df: pd.DataFrame,
    rank_correlations: pd.DataFrame,
    stability_analysis: Dict[str, Any],
    score_variation: pd.DataFrame,
    output_path: str,
) -> None:
    """Generate comprehensive robustness analysis report.

    Args:
        robustness_df: Raw robustness docking results.
        rank_correlations: Pairwise rank correlation results.
        stability_analysis: Output from analyze_rank_stability.
        score_variation: Output from score_variation_analysis.
        output_path: Path to save the JSON report.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    report = {
        "summary": {
            "n_compounds_tested": robustness_df["compound_id"].nunique(),
            "n_conditions_tested": robustness_df["condition_label"].nunique(),
            "total_docking_runs": len(robustness_df),
            "mean_rank_correlation": round(
                rank_correlations["spearman_rho"].mean(), 3
            ) if not rank_correlations.empty else None,
            "n_stable_candidates": stability_analysis.get("n_stable", 0),
            "overall_assessment": stability_analysis.get("assessment", "Unknown"),
        },
        "rank_correlations": rank_correlations.to_dict(orient="records")
        if not rank_correlations.empty else [],
        "stability_analysis": stability_analysis,
        "score_variation_summary": {
            "mean_score_std": round(score_variation["std_score"].mean(), 3)
            if not score_variation.empty else None,
            "mean_score_range": round(score_variation["range_score"].mean(), 3)
            if not score_variation.empty else None,
        },
        "interpretation": (
            "The robustness analysis tests whether docking rankings are sensitive "
            "to reasonable changes in docking box dimensions, search exhaustiveness, "
            "and random seed. Stable rankings (Spearman ρ ≥ 0.85 across conditions) "
            "provide confidence that top candidates are not artifacts of specific "
            "parameter choices. Unstable rankings suggest that the docking results "
            "should be interpreted with additional caution."
        ),
    }

    with open(output_path, "w") as f:
        json.dump(report, f, indent=2, default=str)

    logger.info(f"Robustness report saved to {output_path}")
