"""
Validation module for docking protocol assessment.

Performs redocking of co-crystallized ligands, RMSD calculation,
enrichment analysis, and positive-control benchmarking.
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def calculate_rmsd(
    coords_ref: np.ndarray,
    coords_pred: np.ndarray,
) -> float:
    """Calculate heavy-atom RMSD between reference and predicted coordinates.

    Args:
        coords_ref: Reference coordinates array of shape (N, 3).
        coords_pred: Predicted coordinates array of shape (N, 3).

    Returns:
        RMSD value in Angstroms.

    Raises:
        ValueError: If coordinate arrays have different shapes.
    """
    coords_ref = np.asarray(coords_ref, dtype=np.float64)
    coords_pred = np.asarray(coords_pred, dtype=np.float64)

    if coords_ref.shape != coords_pred.shape:
        raise ValueError(
            f"Coordinate shape mismatch: {coords_ref.shape} vs {coords_pred.shape}"
        )

    diff = coords_ref - coords_pred
    return float(np.sqrt(np.mean(np.sum(diff ** 2, axis=1))))


def extract_heavy_atom_coords(pdb_block: str) -> np.ndarray:
    """Extract heavy atom coordinates from a PDB/PDBQT format text block.

    Args:
        pdb_block: String containing PDB/PDBQT format atom records.

    Returns:
        Numpy array of shape (N, 3) with heavy atom coordinates.
    """
    coords = []
    for line in pdb_block.split("\n"):
        if not line.startswith(("ATOM", "HETATM")):
            continue
        if len(line) < 54:
            continue

        # Check element (column 76-78) or atom name
        element = line[76:78].strip() if len(line) >= 78 else ""
        atom_name = line[12:16].strip()

        # Skip hydrogens
        if element in ("H", "HD") or (not element and atom_name.startswith("H")):
            continue

        try:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            coords.append([x, y, z])
        except ValueError:
            continue

    return np.array(coords)


def extract_coords_from_pdbqt_file(
    pdbqt_path: str,
    model_index: int = 0,
) -> np.ndarray:
    """Extract heavy atom coordinates from a PDBQT file (specific model).

    Args:
        pdbqt_path: Path to PDBQT file (may contain multiple models).
        model_index: Zero-based index of the model to extract.

    Returns:
        Numpy array of shape (N, 3) with heavy atom coordinates.
    """
    path = Path(pdbqt_path)
    if not path.exists():
        raise FileNotFoundError(f"PDBQT file not found: {path}")

    models = []
    current_model = []
    in_model = False

    with open(path, "r") as f:
        for line in f:
            if line.startswith("MODEL"):
                in_model = True
                current_model = []
            elif line.startswith("ENDMDL"):
                models.append("".join(current_model))
                in_model = False
            elif in_model or not models:
                current_model.append(line)

    # If no MODEL/ENDMDL markers, treat entire file as one model
    if not models and current_model:
        models.append("".join(current_model))

    if model_index >= len(models):
        raise IndexError(
            f"Model index {model_index} out of range (file has {len(models)} models)"
        )

    return extract_heavy_atom_coords(models[model_index])


def redocking_validation(
    receptor_pdbqt: str,
    reference_ligand_pdb: str,
    docked_ligand_pdbqt: str,
    compound_name: str = "reference_ligand",
) -> Dict[str, Any]:
    """Perform redocking validation by comparing docked pose to crystal pose.

    Args:
        receptor_pdbqt: Path to receptor PDBQT (for context).
        reference_ligand_pdb: Path to the crystal ligand PDB/PDBQT file.
        docked_ligand_pdbqt: Path to Vina output PDBQT with docked poses.
        compound_name: Name of the ligand being validated.

    Returns:
        Dictionary with RMSD values, assessment, and interpretation.
    """
    logger.info(f"Running redocking validation for {compound_name}")

    # Extract reference coordinates
    try:
        with open(reference_ligand_pdb, "r") as f:
            ref_block = f.read()
        ref_coords = extract_heavy_atom_coords(ref_block)
        logger.info(f"Reference ligand: {len(ref_coords)} heavy atoms")
    except Exception as e:
        return {
            "compound_name": compound_name,
            "status": "error",
            "error": f"Failed to read reference ligand: {e}",
        }

    if len(ref_coords) == 0:
        return {
            "compound_name": compound_name,
            "status": "error",
            "error": "No heavy atoms found in reference ligand",
        }

    # Extract docked pose coordinates (all models)
    results = {
        "compound_name": compound_name,
        "n_reference_atoms": len(ref_coords),
        "pose_rmsds": [],
        "best_rmsd": None,
        "best_pose_index": None,
        "status": "success",
    }

    try:
        path = Path(docked_ligand_pdbqt)
        models = []
        current_model = []

        with open(path, "r") as f:
            for line in f:
                if line.startswith("MODEL"):
                    current_model = []
                elif line.startswith("ENDMDL"):
                    models.append("".join(current_model))
                else:
                    current_model.append(line)

        if not models and current_model:
            models.append("".join(current_model))

        for i, model_block in enumerate(models):
            pred_coords = extract_heavy_atom_coords(model_block)
            if len(pred_coords) == len(ref_coords):
                rmsd = calculate_rmsd(ref_coords, pred_coords)
                results["pose_rmsds"].append({"pose": i + 1, "rmsd": round(rmsd, 3)})
            else:
                logger.warning(
                    f"Atom count mismatch for pose {i + 1}: "
                    f"ref={len(ref_coords)}, pred={len(pred_coords)}"
                )
                results["pose_rmsds"].append({
                    "pose": i + 1,
                    "rmsd": None,
                    "note": f"Atom count mismatch ({len(pred_coords)} vs {len(ref_coords)})",
                })

    except Exception as e:
        results["status"] = "error"
        results["error"] = str(e)
        return results

    # Find best RMSD
    valid_rmsds = [p for p in results["pose_rmsds"] if p["rmsd"] is not None]
    if valid_rmsds:
        best = min(valid_rmsds, key=lambda p: p["rmsd"])
        results["best_rmsd"] = best["rmsd"]
        results["best_pose_index"] = best["pose"]

        # Assessment
        if best["rmsd"] <= 1.0:
            results["assessment"] = "Excellent — pose closely reproduces crystal structure"
        elif best["rmsd"] <= 2.0:
            results["assessment"] = "Good — pose reasonably reproduces crystal binding mode"
        elif best["rmsd"] <= 3.0:
            results["assessment"] = "Marginal — partial reproduction of binding mode"
        else:
            results["assessment"] = "Poor — docked pose deviates substantially from crystal"

        results["interpretation"] = (
            "Redocking RMSD measures the geometric deviation between the computationally "
            "predicted binding pose and the experimentally observed crystal structure pose. "
            "An RMSD ≤ 2.0 Å is conventionally considered successful in molecular docking "
            "validation. However, successful redocking demonstrates only that the scoring "
            "function can identify the correct binding geometry for THIS specific ligand "
            "in THIS specific receptor conformation. It does not guarantee accurate ranking "
            "of other compounds or prediction of binding affinity."
        )
    else:
        results["status"] = "error"
        results["error"] = "No valid RMSD calculations (atom count mismatches)"

    return results


def calculate_enrichment(
    docking_results_df: pd.DataFrame,
    known_active_ids: List[str],
    score_col: str = "best_score",
    id_col: str = "compound_id",
    percentages: List[float] = None,
) -> Dict[str, Any]:
    """Calculate enrichment factors for known active compounds.

    The enrichment factor (EF) at X% measures how many times more likely
    we are to find known actives in the top X% of the ranked list compared
    to random selection.

    EF(X%) = (actives_found_in_top_X% / N_top_X%) / (total_actives / total_compounds)

    Args:
        docking_results_df: DataFrame with docking results.
        known_active_ids: List of compound IDs that are known actives.
        score_col: Column name for the docking score (lower = better).
        id_col: Column name for compound identifiers.
        percentages: List of percentages at which to calculate EF.

    Returns:
        Dictionary with enrichment metrics.
    """
    if percentages is None:
        percentages = [1.0, 2.0, 5.0, 10.0, 20.0]

    df = docking_results_df.copy()
    df = df.dropna(subset=[score_col])
    df = df.sort_values(score_col, ascending=True)  # Lower score = better
    df = df.reset_index(drop=True)

    total = len(df)
    total_actives = len([x for x in df[id_col] if x in known_active_ids])

    if total_actives == 0:
        logger.warning("No known actives found in docking results")
        return {"error": "No known actives in results", "enrichment_factors": {}}

    expected_rate = total_actives / total

    enrichment_factors = {}
    for pct in percentages:
        n_top = max(1, int(total * pct / 100.0))
        top_ids = set(df.head(n_top)[id_col])
        actives_found = len(top_ids.intersection(known_active_ids))
        observed_rate = actives_found / n_top if n_top > 0 else 0.0
        ef = observed_rate / expected_rate if expected_rate > 0 else 0.0

        enrichment_factors[f"EF_{pct}%"] = {
            "n_compounds_screened": n_top,
            "actives_found": actives_found,
            "enrichment_factor": round(ef, 2),
            "expected_actives_random": round(n_top * expected_rate, 2),
        }

    # Ranks of known actives
    active_ranks = []
    for idx, row in df.iterrows():
        if row[id_col] in known_active_ids:
            active_ranks.append({
                "compound_id": row[id_col],
                "rank": idx + 1,
                "percentile": round((idx + 1) / total * 100, 1),
                "score": row[score_col],
            })

    return {
        "total_compounds": total,
        "total_actives": total_actives,
        "expected_hit_rate": round(expected_rate * 100, 2),
        "enrichment_factors": enrichment_factors,
        "active_ranks": active_ranks,
        "median_active_rank": (
            round(np.median([a["rank"] for a in active_ranks]), 1)
            if active_ranks else None
        ),
        "median_active_percentile": (
            round(np.median([a["percentile"] for a in active_ranks]), 1)
            if active_ranks else None
        ),
    }


def generate_validation_report(
    redocking_result: Dict[str, Any],
    enrichment_result: Dict[str, Any],
    output_path: str,
) -> None:
    """Generate a comprehensive validation report as JSON.

    Args:
        redocking_result: Output from redocking_validation.
        enrichment_result: Output from calculate_enrichment.
        output_path: Path to save the JSON report.
    """
    import json

    report = {
        "redocking_validation": redocking_result,
        "enrichment_analysis": enrichment_result,
        "caveats": [
            "Successful redocking demonstrates the scoring function can identify "
            "the correct geometry for one known ligand. It does not validate "
            "the accuracy of binding affinity predictions.",
            "Enrichment factors depend on the composition and size of the screening "
            "library and should not be compared across differently composed libraries.",
            "Docking scores are approximate and should not be interpreted as "
            "experimental binding free energies.",
        ],
    }

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w") as f:
        json.dump(report, f, indent=2, default=str)

    logger.info(f"Validation report saved to {output_path}")
