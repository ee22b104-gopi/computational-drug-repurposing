"""
Protein-ligand interaction analysis module.

Uses ProLIF and/or manual distance-based analysis to identify
hydrogen bonds, hydrophobic contacts, aromatic interactions,
salt bridges, and interacting residues from docked poses.
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def parse_pdbqt_poses(pdbqt_path: str) -> List[str]:
    """Parse a multi-model PDBQT file into individual pose blocks.

    Args:
        pdbqt_path: Path to multi-model PDBQT output from Vina.

    Returns:
        List of PDB-format string blocks, one per pose.
    """
    pdbqt_path = Path(pdbqt_path)
    if not pdbqt_path.exists():
        raise FileNotFoundError(f"PDBQT file not found: {pdbqt_path}")

    poses = []
    current_pose = []
    with open(pdbqt_path, "r") as f:
        for line in f:
            if line.startswith("MODEL"):
                current_pose = []
            elif line.startswith("ENDMDL"):
                poses.append("".join(current_pose))
            elif line.startswith(("ATOM", "HETATM")):
                current_pose.append(line)

    if not poses and current_pose:
        poses.append("".join(current_pose))

    logger.info(f"Parsed {len(poses)} poses from {pdbqt_path.name}")
    return poses


def calculate_contacts(
    receptor_pdb: str,
    ligand_pose_block: str,
    distance_cutoff: float = 4.0,
) -> List[Dict[str, Any]]:
    """Calculate atom-level contacts between receptor and a ligand pose.

    Uses simple Euclidean distance between heavy atoms.

    Args:
        receptor_pdb: Path to receptor PDB file.
        ligand_pose_block: PDB-format string of the ligand pose.
        distance_cutoff: Maximum distance (Å) for a contact.

    Returns:
        List of contact dictionaries with residue info and distances.
    """
    try:
        from Bio.PDB import PDBParser

        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("receptor", receptor_pdb)
        model = structure[0]
    except Exception as e:
        logger.error(f"Failed to parse receptor PDB: {e}")
        return []

    # Parse ligand atoms from the pose block
    ligand_atoms = []
    for line in ligand_pose_block.split("\n"):
        if line.startswith(("ATOM", "HETATM")) and len(line) >= 54:
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                atom_name = line[12:16].strip()
                element = line[76:78].strip() if len(line) >= 78 else atom_name[0]
                if element not in ("H", "HD"):
                    ligand_atoms.append({
                        "name": atom_name,
                        "coord": np.array([x, y, z]),
                        "element": element,
                    })
            except (ValueError, IndexError):
                continue

    if not ligand_atoms:
        logger.warning("No ligand heavy atoms parsed from pose block")
        return []

    contacts = []
    for chain in model:
        for residue in chain:
            res_name = residue.get_resname()
            res_id = residue.get_id()[1]
            chain_id = chain.get_id()

            # Skip water and non-standard residues
            if res_name in ("HOH", "WAT"):
                continue

            for atom in residue:
                if atom.element == "H":
                    continue
                rec_coord = atom.get_vector().get_array()

                for lig_atom in ligand_atoms:
                    dist = np.linalg.norm(rec_coord - lig_atom["coord"])
                    if dist <= distance_cutoff:
                        contacts.append({
                            "residue_name": res_name,
                            "residue_number": res_id,
                            "chain": chain_id,
                            "receptor_atom": atom.get_name(),
                            "receptor_element": atom.element,
                            "ligand_atom": lig_atom["name"],
                            "ligand_element": lig_atom["element"],
                            "distance": round(dist, 2),
                        })

    logger.info(f"Found {len(contacts)} contacts within {distance_cutoff} Å")
    return contacts


def classify_interaction(
    receptor_atom: str,
    receptor_element: str,
    ligand_atom: str,
    ligand_element: str,
    distance: float,
) -> str:
    """Classify a contact into an interaction type.

    This is a simplified heuristic classification. For rigorous analysis,
    use ProLIF or PLIP.

    Args:
        receptor_atom: Receptor atom name.
        receptor_element: Receptor element symbol.
        ligand_atom: Ligand atom name.
        ligand_element: Ligand element symbol.
        distance: Distance in Angstroms.

    Returns:
        Interaction type string.
    """
    h_bond_donors = {"N", "O", "S"}
    h_bond_acceptors = {"N", "O", "S", "F"}
    hydrophobic = {"C"}
    charged_positive = {"N"}  # Simplified
    charged_negative = {"O"}  # Simplified (carboxylate O)

    # Hydrogen bond (simplified: donor-acceptor within 3.5 Å)
    if distance <= 3.5:
        if (receptor_element in h_bond_donors and ligand_element in h_bond_acceptors):
            return "hydrogen_bond"
        if (receptor_element in h_bond_acceptors and ligand_element in h_bond_donors):
            return "hydrogen_bond"

    # Salt bridge (charged atoms within 4.0 Å)
    if distance <= 4.0:
        if (receptor_element in charged_positive and ligand_element in charged_negative):
            return "salt_bridge"
        if (receptor_element in charged_negative and ligand_element in charged_positive):
            return "salt_bridge"

    # Hydrophobic contact (C-C within 4.0 Å)
    if distance <= 4.0:
        if receptor_element in hydrophobic and ligand_element in hydrophobic:
            return "hydrophobic"

    return "van_der_waals"


def analyze_pose_interactions(
    receptor_pdb: str,
    ligand_pose_block: str,
    distance_cutoff: float = 4.5,
) -> Dict[str, Any]:
    """Analyze all interactions for a single ligand pose.

    Args:
        receptor_pdb: Path to receptor PDB file.
        ligand_pose_block: PDB-format string of the ligand pose.
        distance_cutoff: Maximum distance for contacts.

    Returns:
        Dictionary with interaction summary.
    """
    contacts = calculate_contacts(receptor_pdb, ligand_pose_block, distance_cutoff)

    if not contacts:
        return {
            "n_contacts": 0,
            "interacting_residues": [],
            "interaction_types": {},
            "contacts": [],
        }

    # Classify each contact
    for contact in contacts:
        contact["interaction_type"] = classify_interaction(
            contact["receptor_atom"],
            contact["receptor_element"],
            contact["ligand_atom"],
            contact["ligand_element"],
            contact["distance"],
        )

    # Summarize
    residues = list(set(
        (c["residue_name"], c["residue_number"], c["chain"])
        for c in contacts
    ))
    residues.sort(key=lambda r: r[1])

    interaction_counts = {}
    for c in contacts:
        itype = c["interaction_type"]
        interaction_counts[itype] = interaction_counts.get(itype, 0) + 1

    return {
        "n_contacts": len(contacts),
        "interacting_residues": [
            {"name": r[0], "number": r[1], "chain": r[2]} for r in residues
        ],
        "interaction_types": interaction_counts,
        "contacts": contacts,
    }


def generate_interaction_fingerprint(
    interaction_data: Dict[str, Any],
    reference_residues: Optional[List[int]] = None,
) -> Dict[int, List[str]]:
    """Generate a residue-level interaction fingerprint.

    Args:
        interaction_data: Output from analyze_pose_interactions.
        reference_residues: List of residue numbers to include in fingerprint.
            If None, uses all interacting residues.

    Returns:
        Dict mapping residue number to list of interaction types observed.
    """
    fingerprint = {}
    for contact in interaction_data.get("contacts", []):
        res_num = contact["residue_number"]
        itype = contact["interaction_type"]

        if reference_residues and res_num not in reference_residues:
            continue

        if res_num not in fingerprint:
            fingerprint[res_num] = set()
        fingerprint[res_num].add(itype)

    return {k: sorted(list(v)) for k, v in sorted(fingerprint.items())}


def check_binding_site_coverage(
    interaction_data: Dict[str, Any],
    catalytic_triad: List[int] = None,
    s1_pocket: List[int] = None,
    s2_pocket: List[int] = None,
    s2_extensive: List[int] = None,
) -> Dict[str, Any]:
    """Check whether a pose interacts with known DPP-4 binding site residues.

    Args:
        interaction_data: Output from analyze_pose_interactions.
        catalytic_triad: Residue numbers for catalytic triad (default: [630, 708, 740]).
        s1_pocket: Residue numbers for S1 pocket.
        s2_pocket: Residue numbers for S2 pocket.
        s2_extensive: Residue numbers for S2-extensive pocket.

    Returns:
        Dictionary with coverage metrics for each subsite.
    """
    if catalytic_triad is None:
        catalytic_triad = [630, 708, 740]
    if s1_pocket is None:
        s1_pocket = [631, 656, 659, 662, 666, 711]
    if s2_pocket is None:
        s2_pocket = [125, 205, 206]
    if s2_extensive is None:
        s2_extensive = [207, 209, 357, 358]

    interacting_nums = set(
        r["number"] for r in interaction_data.get("interacting_residues", [])
    )

    def subsite_coverage(subsite_residues: List[int]) -> Dict[str, Any]:
        contacted = interacting_nums.intersection(subsite_residues)
        return {
            "total_residues": len(subsite_residues),
            "contacted_residues": sorted(list(contacted)),
            "n_contacted": len(contacted),
            "coverage_fraction": (
                len(contacted) / len(subsite_residues) if subsite_residues else 0.0
            ),
        }

    return {
        "catalytic_triad": subsite_coverage(catalytic_triad),
        "s1_pocket": subsite_coverage(s1_pocket),
        "s2_pocket": subsite_coverage(s2_pocket),
        "s2_extensive": subsite_coverage(s2_extensive),
        "total_known_residues_contacted": len(
            interacting_nums.intersection(
                set(catalytic_triad + s1_pocket + s2_pocket + s2_extensive)
            )
        ),
    }


def analyze_library_interactions(
    receptor_pdb: str,
    docking_results_df: pd.DataFrame,
    docking_output_dir: str,
    top_n: int = 30,
    reference_residues: Optional[List[int]] = None,
) -> pd.DataFrame:
    """Analyze interactions for top-ranked compounds from docking.

    Args:
        receptor_pdb: Path to receptor PDB file.
        docking_results_df: DataFrame with docking results (must have 'compound_id',
            'compound_name', 'best_score' columns).
        docking_output_dir: Directory containing docked PDBQT files.
        top_n: Number of top-ranked compounds to analyze.
        reference_residues: DPP-4 binding site residue numbers.

    Returns:
        DataFrame with interaction analysis for each compound.
    """
    if reference_residues is None:
        reference_residues = [125, 205, 206, 207, 209, 357, 358,
                              547, 630, 631, 656, 659, 662, 666,
                              708, 711, 740]

    sorted_df = docking_results_df.sort_values("best_score").head(top_n)
    results = []

    for _, row in sorted_df.iterrows():
        compound_id = row["compound_id"]
        pdbqt_path = Path(docking_output_dir) / f"{compound_id}_out.pdbqt"

        if not pdbqt_path.exists():
            logger.warning(f"Output PDBQT not found for {compound_id}")
            results.append({
                "compound_id": compound_id,
                "compound_name": row.get("compound_name", ""),
                "best_score": row["best_score"],
                "n_contacts": 0,
                "n_interacting_residues": 0,
                "interaction_types": "{}",
                "interacting_residue_numbers": "",
                "binding_site_coverage": 0.0,
                "status": "pdbqt_not_found",
            })
            continue

        try:
            poses = parse_pdbqt_poses(str(pdbqt_path))
            if not poses:
                raise ValueError("No poses found in PDBQT file")

            # Analyze best (first) pose
            interactions = analyze_pose_interactions(
                receptor_pdb, poses[0]
            )
            coverage = check_binding_site_coverage(interactions)
            fingerprint = generate_interaction_fingerprint(
                interactions, reference_residues
            )

            results.append({
                "compound_id": compound_id,
                "compound_name": row.get("compound_name", ""),
                "best_score": row["best_score"],
                "n_contacts": interactions["n_contacts"],
                "n_interacting_residues": len(interactions["interacting_residues"]),
                "interaction_types": str(interactions["interaction_types"]),
                "interacting_residue_numbers": ",".join(
                    str(r["number"]) for r in interactions["interacting_residues"]
                ),
                "catalytic_triad_contacts": coverage["catalytic_triad"]["n_contacted"],
                "s1_contacts": coverage["s1_pocket"]["n_contacted"],
                "s2_contacts": coverage["s2_pocket"]["n_contacted"],
                "s2_ext_contacts": coverage["s2_extensive"]["n_contacted"],
                "binding_site_coverage": coverage["total_known_residues_contacted"],
                "fingerprint": str(fingerprint),
                "status": "success",
            })
        except Exception as e:
            logger.error(f"Interaction analysis failed for {compound_id}: {e}")
            results.append({
                "compound_id": compound_id,
                "compound_name": row.get("compound_name", ""),
                "best_score": row["best_score"],
                "n_contacts": 0,
                "n_interacting_residues": 0,
                "interaction_types": "{}",
                "interacting_residue_numbers": "",
                "binding_site_coverage": 0.0,
                "status": f"error: {str(e)}",
            })

    return pd.DataFrame(results)


def try_prolif_analysis(
    receptor_pdb: str,
    ligand_sdf: str,
) -> Optional[pd.DataFrame]:
    """Attempt interaction fingerprint analysis using ProLIF if available.

    Args:
        receptor_pdb: Path to receptor PDB file.
        ligand_sdf: Path to ligand SDF file.

    Returns:
        ProLIF interaction DataFrame or None if ProLIF is not available.
    """
    try:
        import prolif as plf
        import MDAnalysis as mda

        protein = mda.Universe(receptor_pdb)
        protein = plf.Molecule.from_mda(protein)

        ligand_universe = mda.Universe(ligand_sdf)
        ligand = plf.Molecule.from_mda(ligand_universe)

        fp = plf.Fingerprint()
        fp.run_from_iterable([ligand], protein)

        return fp.to_dataframe()

    except ImportError:
        logger.info("ProLIF not available; using distance-based analysis")
        return None
    except Exception as e:
        logger.warning(f"ProLIF analysis failed: {e}")
        return None
