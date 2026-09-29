"""
Protein structure operations module for computational drug repurposing.
Handles downloading, cleaning, and preparing proteins for docking.
"""

import json
import logging
import subprocess
import requests
import numpy as np
from pathlib import Path
from typing import List, Tuple, Optional, Union
import Bio.PDB
from Bio.PDB import PDBParser, PDBIO, Select, NeighborSearch
from Bio.PDB.Structure import Structure
from Bio.PDB.Residue import Residue

logger = logging.getLogger(__name__)

def download_pdb(pdb_id: str, output_dir: Union[str, Path]) -> Path:
    """
    Downloads PDB file from RCSB PDB.
    
    Args:
        pdb_id: 4-character PDB ID
        output_dir: Directory to save the PDB file
        
    Returns:
        Path to the downloaded PDB file
    """
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    
    pdb_id = pdb_id.upper()
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    output_path = out_dir / f"{pdb_id}.pdb"
    
    try:
        logger.info(f"Downloading PDB {pdb_id} from {url}")
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        
        with open(output_path, 'wb') as f:
            f.write(response.content)
            
        logger.info(f"Successfully downloaded {pdb_id} to {output_path}")
        return output_path
    except requests.RequestException as e:
        logger.error(f"Failed to download PDB {pdb_id}: {e}")
        raise

class ChainSelect(Select):
    """Selector class for Bio.PDB to extract specific chains."""
    def __init__(self, chain_id: str):
        self.chain_id = chain_id
        
    def accept_chain(self, chain):
        return 1 if chain.get_id() == self.chain_id else 0

def extract_chain(pdb_path: Union[str, Path], chain_id: str, output_path: Union[str, Path]) -> Path:
    """
    Uses BioPython to extract a single chain.
    
    Args:
        pdb_path: Path to input PDB file
        chain_id: Chain identifier to extract
        output_path: Path for output PDB file
        
    Returns:
        Path to output PDB file
    """
    pdb_path = Path(pdb_path)
    output_path = Path(output_path)
    
    try:
        logger.info(f"Extracting chain {chain_id} from {pdb_path}")
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("protein", str(pdb_path))
        
        io = PDBIO()
        io.set_structure(structure)
        io.save(str(output_path), ChainSelect(chain_id))
        
        logger.info(f"Saved extracted chain to {output_path}")
        return output_path
    except Exception as e:
        logger.error(f"Failed to extract chain {chain_id} from {pdb_path}: {e}")
        raise

class CleanProteinSelect(Select):
    """Selector class for Bio.PDB to clean protein structures."""
    def __init__(self, remove_water: bool = True, remove_heterogens: bool = True, keep_ligand_code: Optional[str] = None):
        self.remove_water = remove_water
        self.remove_heterogens = remove_heterogens
        self.keep_ligand_code = keep_ligand_code
        
    def accept_residue(self, residue):
        res_id = residue.get_id()
        hetfield, seqid, icode = res_id
        
        # Remove water
        if self.remove_water and hetfield == "W":
            return 0
            
        # Handle heterogens
        if hetfield != " ":
            # If it's the ligand we want to keep
            if self.keep_ligand_code and residue.get_resname() == self.keep_ligand_code:
                return 1
            # If we want to remove other heterogens
            if self.remove_heterogens:
                return 0
                
        return 1
        
    def accept_atom(self, atom):
        # Keep only alternate conformation A or atoms without alternate conformation
        altloc = atom.get_altloc()
        if altloc != ' ' and altloc != 'A':
            return 0
            
        # Remove hydrogens for later re-addition
        if atom.element == 'H':
            return 0
            
        return 1

def clean_protein(pdb_path: Union[str, Path], output_path: Union[str, Path], 
                  remove_water: bool = True, remove_heterogens: bool = True, 
                  keep_ligand_code: Optional[str] = None) -> Path:
    """
    Cleans PDB: removes waters, HETATM records (except specified ligand), 
    alternate conformations, and hydrogens.
    """
    pdb_path = Path(pdb_path)
    output_path = Path(output_path)
    
    try:
        logger.info(f"Cleaning protein {pdb_path}")
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("protein", str(pdb_path))
        
        # Strip alternate conformations from names
        for atom in structure.get_atoms():
            if atom.get_altloc() == 'A':
                atom.set_altloc(' ')
        
        io = PDBIO()
        io.set_structure(structure)
        selector = CleanProteinSelect(remove_water, remove_heterogens, keep_ligand_code)
        io.save(str(output_path), selector)
        
        logger.info(f"Saved cleaned protein to {output_path}")
        return output_path
    except Exception as e:
        logger.error(f"Error cleaning protein {pdb_path}: {e}")
        raise

def add_hydrogens(pdb_path: Union[str, Path], output_path: Union[str, Path], ph: float = 7.4) -> Path:
    """
    Adds hydrogens to protein. Placeholder since Meeko/ADFR handles this during PDBQT conversion.
    """
    logger.info(f"add_hydrogens called for {pdb_path}. Relying on PDBQT conversion tools to handle protonation.")
    # For now, just copy the file over or do a simple save.
    output_path = Path(output_path)
    if pdb_path != output_path:
        output_path.write_text(Path(pdb_path).read_text())
    return output_path

def prepare_receptor_pdbqt(pdb_path: Union[str, Path], output_path: Union[str, Path]) -> Path:
    """
    Converts clean PDB to PDBQT format using Meeko or fallback.
    """
    pdb_path = Path(pdb_path)
    output_path = Path(output_path)
    
    logger.info(f"Preparing receptor PDBQT for {pdb_path}")
    
    try:
        # Try meeko programmatic approach first
        from meeko import PDBQTMolecule
        from meeko.preparation import ReceptorPreparation
        
        prep = ReceptorPreparation()
        # This is a hypothetical API, real Meeko receptor prep might differ or be handled mostly by ADFRSuite
        # Meeko provides PDBQTReceptor in recent versions for simpler receptor prep
        try:
            from meeko import PDBQTReceptor
            receptor = PDBQTReceptor(str(pdb_path))
            receptor.write_pdbqt(str(output_path))
            logger.info("Successfully prepared PDBQT using meeko.PDBQTReceptor")
            return output_path
        except ImportError:
            pass
            
    except ImportError:
        logger.warning("Meeko not available for programmatic receptor preparation.")
        
    # Fallback to subprocess
    logger.info("Trying subprocess fallback for mk_prepare_receptor.py")
    import sys
    import os
    script_paths = [
        "mk_prepare_receptor.exe",
        "mk_prepare_receptor.py",
        r"C:\Users\gaadi\AppData\Local\Packages\PythonSoftwareFoundation.Python.3.13_qbz5n2kfra8p0\LocalCache\local-packages\Python313\Scripts\mk_prepare_receptor.exe",
    ]
    for sp in script_paths:
        try:
            cmd = [sp, "--read_pdb", str(pdb_path), "-p", str(output_path)]
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            logger.info(f"Successfully prepared PDBQT using {sp}")
            return output_path
        except (subprocess.CalledProcessError, FileNotFoundError):
            continue
            
    logger.error("All receptor preparation methods failed.")
    raise RuntimeError("Failed to prepare receptor PDBQT.")

class LigandSelect(Select):
    """Selector class to extract specific ligand."""
    def __init__(self, ligand_code: str, chain_id: str):
        self.ligand_code = ligand_code
        self.chain_id = chain_id
        
    def accept_residue(self, residue):
        parent_chain = residue.get_parent().get_id()
        if parent_chain != self.chain_id:
            return 0
        if residue.get_resname() == self.ligand_code:
            return 1
        return 0

def extract_ligand(pdb_path: Union[str, Path], ligand_code: str, chain_id: str, output_path: Union[str, Path]) -> Path:
    """
    Extracts the co-crystallized ligand from PDB into a separate PDB file.
    """
    pdb_path = Path(pdb_path)
    output_path = Path(output_path)
    
    try:
        logger.info(f"Extracting ligand {ligand_code} from chain {chain_id} in {pdb_path}")
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("protein", str(pdb_path))
        
        io = PDBIO()
        io.set_structure(structure)
        io.save(str(output_path), LigandSelect(ligand_code, chain_id))
        
        logger.info(f"Saved ligand to {output_path}")
        return output_path
    except Exception as e:
        logger.error(f"Error extracting ligand: {e}")
        raise

def get_ligand_center(pdb_path: Union[str, Path], ligand_code: str, chain_id: str) -> Tuple[float, float, float]:
    """
    Calculates the geometric center (centroid) of the co-crystallized ligand heavy atoms.
    """
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("protein", str(pdb_path))
        
        coords = []
        for model in structure:
            for chain in model:
                if chain.get_id() == chain_id:
                    for residue in chain:
                        if residue.get_resname() == ligand_code:
                            for atom in residue:
                                if atom.element != 'H':
                                    coords.append(atom.coord)
                                    
        if not coords:
            raise ValueError(f"Ligand {ligand_code} not found in chain {chain_id}")
            
        center = np.mean(coords, axis=0)
        logger.info(f"Calculated center for {ligand_code}: {center}")
        return tuple(float(x) for x in center)
    except Exception as e:
        logger.error(f"Error calculating ligand center: {e}")
        raise

def get_binding_site_residues(pdb_path: Union[str, Path], ligand_code: str, chain_id: str, distance_cutoff: float = 5.0) -> List[Tuple[str, int]]:
    """
    Identifies protein residues within distance_cutoff Angstroms of any ligand atom.
    """
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("protein", str(pdb_path))
        
        ligand_atoms = []
        protein_atoms = []
        
        for model in structure:
            for chain in model:
                if chain.get_id() == chain_id:
                    for residue in chain:
                        if residue.get_resname() == ligand_code:
                            ligand_atoms.extend(residue.get_atoms())
                        elif residue.get_id()[0] == ' ': # Standard amino acid
                            protein_atoms.extend(residue.get_atoms())
                            
        if not ligand_atoms:
            raise ValueError(f"Ligand {ligand_code} not found in chain {chain_id}")
            
        ns = NeighborSearch(protein_atoms)
        nearby_residues = set()
        
        for atom in ligand_atoms:
            nearby_atoms = ns.search(atom.coord, distance_cutoff)
            for nearby_atom in nearby_atoms:
                res = nearby_atom.get_parent()
                nearby_residues.add((res.get_resname(), res.get_id()[1]))
                
        result = sorted(list(nearby_residues), key=lambda x: x[1])
        logger.info(f"Found {len(result)} binding site residues within {distance_cutoff}A")
        return result
    except Exception as e:
        logger.error(f"Error finding binding site residues: {e}")
        raise

def generate_protein_report(pdb_path: Union[str, Path], chain_id: str, ligand_code: str, output_path: Union[str, Path]) -> None:
    """
    Generates a JSON report with protein info.
    """
    try:
        pdb_path = Path(pdb_path)
        
        center = get_ligand_center(pdb_path, ligand_code, chain_id)
        nearby_res = get_binding_site_residues(pdb_path, ligand_code, chain_id)
        
        report = {
            "pdb_file": pdb_path.name,
            "chain_id": chain_id,
            "ligand_code": ligand_code,
            "binding_site_center": {
                "x": center[0],
                "y": center[1],
                "z": center[2]
            },
            "binding_site_residues": [f"{res[0]}{res[1]}" for res in nearby_res],
            "binding_site_residue_count": len(nearby_res)
        }
        
        output_path = Path(output_path)
        with open(output_path, 'w') as f:
            json.dump(report, f, indent=4)
            
        logger.info(f"Generated protein report at {output_path}")
    except Exception as e:
        logger.error(f"Error generating protein report: {e}")
        raise
