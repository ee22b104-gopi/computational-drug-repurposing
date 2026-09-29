"""
descriptors.py

This module calculates molecular descriptors using RDKit for computational drug discovery pipelines.
"""

import json
import logging
from typing import Dict, List, Optional, Tuple, Any, Union

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors
from rdkit import DataStructs
from tqdm import tqdm

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def calculate_descriptors(smiles: str, name: Optional[str] = None) -> Dict[str, Any]:
    """
    Calculates key descriptors for a single molecule from SMILES.
    
    Args:
        smiles (str): SMILES string of the molecule.
        name (Optional[str]): Optional name of the molecule.
        
    Returns:
        Dict[str, Any]: Dictionary containing molecular descriptors. 
                        Returns None for all descriptor values if SMILES is invalid.
    """
    mol = Chem.MolFromSmiles(smiles)
    
    if mol is None:
        logger.warning(f"Invalid SMILES string: {smiles}")
        return {
            'molecular_weight': None,
            'logp': None,
            'hbd': None,
            'hba': None,
            'rotatable_bonds': None,
            'tpsa': None,
            'num_rings': None,
            'num_aromatic_rings': None,
            'num_heavy_atoms': None,
            'molecular_formula': None,
            'fraction_csp3': None
        }
        
    return {
        'molecular_weight': Descriptors.ExactMolWt(mol),
        'logp': Descriptors.MolLogP(mol),
        'hbd': Descriptors.NumHDonors(mol),
        'hba': Descriptors.NumHAcceptors(mol),
        'rotatable_bonds': Descriptors.NumRotatableBonds(mol),
        'tpsa': Descriptors.TPSA(mol),
        'num_rings': Descriptors.RingCount(mol),
        'num_aromatic_rings': Descriptors.NumAromaticRings(mol),
        'num_heavy_atoms': Descriptors.HeavyAtomCount(mol),
        'molecular_formula': rdMolDescriptors.CalcMolFormula(mol),
        'fraction_csp3': Descriptors.FractionCSP3(mol)
    }

def calculate_library_descriptors(compounds_df: pd.DataFrame, 
                                  smiles_col: str = 'smiles', 
                                  name_col: str = 'name', 
                                  id_col: str = 'compound_id') -> pd.DataFrame:
    """
    Calculates descriptors for all compounds in a DataFrame.
    
    Args:
        compounds_df (pd.DataFrame): DataFrame containing compound data.
        smiles_col (str): Column name containing SMILES strings.
        name_col (str): Column name containing molecule names.
        id_col (str): Column name containing compound IDs.
        
    Returns:
        pd.DataFrame: New DataFrame with descriptor columns added.
    """
    logger.info(f"Calculating descriptors for {len(compounds_df)} compounds...")
    
    results = []
    for _, row in tqdm(compounds_df.iterrows(), total=len(compounds_df), desc="Processing Compounds"):
        smiles = row.get(smiles_col)
        name = row.get(name_col) if name_col in row else None
        
        if pd.isna(smiles) or not isinstance(smiles, str):
            continue
            
        desc = calculate_descriptors(smiles, name=name)
        
        # Add original row data
        row_dict = row.to_dict()
        row_dict.update(desc)
        results.append(row_dict)
        
    return pd.DataFrame(results)

def check_lipinski(descriptors_dict: Dict[str, Any]) -> Dict[str, Any]:
    """
    Checks Lipinski's Rule of Five for a molecule based on its descriptors.
    
    Args:
        descriptors_dict (Dict[str, Any]): Dictionary of molecular descriptors.
        
    Returns:
        Dict[str, Any]: Dictionary with pass/fail status, number of violations, and violation descriptions.
    """
    violations = []
    
    if descriptors_dict.get('molecular_weight') is not None and descriptors_dict['molecular_weight'] > 500:
        violations.append('Molecular weight > 500')
    if descriptors_dict.get('logp') is not None and descriptors_dict['logp'] > 5:
        violations.append('LogP > 5')
    if descriptors_dict.get('hbd') is not None and descriptors_dict['hbd'] > 5:
        violations.append('H-bond donors > 5')
    if descriptors_dict.get('hba') is not None and descriptors_dict['hba'] > 10:
        violations.append('H-bond acceptors > 10')
        
    return {
        'passes': len(violations) <= 1,
        'n_violations': len(violations),
        'violations': violations
    }

def check_drug_likeness(descriptors_dict: Dict[str, Any]) -> Dict[str, Any]:
    """
    Extended drug-likeness check combining Lipinski and Veber rules.
    Veber rules: rotatable bonds <= 10, TPSA <= 140.
    
    Args:
        descriptors_dict (Dict[str, Any]): Dictionary of molecular descriptors.
        
    Returns:
        Dict[str, Any]: Dictionary with pass/fail status and violation descriptions.
    """
    lipinski_result = check_lipinski(descriptors_dict)
    violations = lipinski_result['violations'].copy()
    
    if descriptors_dict.get('rotatable_bonds') is not None and descriptors_dict['rotatable_bonds'] > 10:
        violations.append('Rotatable bonds > 10')
    if descriptors_dict.get('tpsa') is not None and descriptors_dict['tpsa'] > 140:
        violations.append('TPSA > 140')
        
    return {
        'passes': len(violations) == 0, # strict check for drug-likeness
        'n_violations': len(violations),
        'violations': violations
    }

def filter_compounds(descriptors_df: pd.DataFrame, 
                     max_mw: float = 500.0, 
                     max_logp: float = 5.0, 
                     max_hbd: int = 5, 
                     max_hba: int = 10, 
                     max_rotatable: int = 10, 
                     max_tpsa: float = 140.0) -> Tuple[pd.DataFrame, pd.DataFrame, Dict[str, Any]]:
    """
    Filters DataFrame based on descriptor thresholds.
    
    Args:
        descriptors_df (pd.DataFrame): DataFrame containing molecular descriptors.
        max_mw (float): Maximum molecular weight.
        max_logp (float): Maximum LogP.
        max_hbd (int): Maximum H-bond donors.
        max_hba (int): Maximum H-bond acceptors.
        max_rotatable (int): Maximum rotatable bonds.
        max_tpsa (float): Maximum TPSA.
        
    Returns:
        Tuple[pd.DataFrame, pd.DataFrame, Dict[str, Any]]: Passing DataFrame, failing DataFrame, and report dictionary.
    """
    logger.info("Filtering compounds based on descriptors...")
    
    # Create mask for passing compounds
    mask = pd.Series([True] * len(descriptors_df), index=descriptors_df.index)
    
    if 'molecular_weight' in descriptors_df.columns:
        mask &= (descriptors_df['molecular_weight'] <= max_mw)
    if 'logp' in descriptors_df.columns:
        mask &= (descriptors_df['logp'] <= max_logp)
    if 'hbd' in descriptors_df.columns:
        mask &= (descriptors_df['hbd'] <= max_hbd)
    if 'hba' in descriptors_df.columns:
        mask &= (descriptors_df['hba'] <= max_hba)
    if 'rotatable_bonds' in descriptors_df.columns:
        mask &= (descriptors_df['rotatable_bonds'] <= max_rotatable)
    if 'tpsa' in descriptors_df.columns:
        mask &= (descriptors_df['tpsa'] <= max_tpsa)
        
    passing_df = descriptors_df[mask].copy()
    failing_df = descriptors_df[~mask].copy()
    
    report = {
        'total_compounds': len(descriptors_df),
        'passing_compounds': len(passing_df),
        'failing_compounds': len(failing_df),
        'pass_rate': len(passing_df) / len(descriptors_df) if len(descriptors_df) > 0 else 0.0
    }
    
    return passing_df, failing_df, report

def generate_descriptor_report(descriptors_df: pd.DataFrame, output_path: str) -> None:
    """
    Generates a JSON summary report of the descriptors.
    
    Args:
        descriptors_df (pd.DataFrame): DataFrame containing molecular descriptors.
        output_path (str): Path to save the JSON report.
    """
    numeric_cols = ['molecular_weight', 'logp', 'hbd', 'hba', 'rotatable_bonds', 'tpsa', 
                    'num_rings', 'num_aromatic_rings', 'num_heavy_atoms', 'fraction_csp3']
    
    available_cols = [col for col in numeric_cols if col in descriptors_df.columns]
    
    report = {}
    
    for col in available_cols:
        col_data = descriptors_df[col].dropna()
        if len(col_data) > 0:
            report[col] = {
                'mean': float(np.mean(col_data)),
                'median': float(np.median(col_data)),
                'std': float(np.std(col_data)),
                'min': float(np.min(col_data)),
                'max': float(np.max(col_data))
            }
            
    # Calculate Lipinski pass rate
    if all(c in descriptors_df.columns for c in ['molecular_weight', 'logp', 'hbd', 'hba']):
        pass_count = 0
        total_count = 0
        for _, row in descriptors_df.iterrows():
            if pd.isna(row.get('molecular_weight')):
                continue
            total_count += 1
            if check_lipinski(row.to_dict())['passes']:
                pass_count += 1
                
        report['lipinski_pass_rate'] = pass_count / total_count if total_count > 0 else 0.0

    try:
        with open(output_path, 'w') as f:
            json.dump(report, f, indent=4)
        logger.info(f"Descriptor report saved to {output_path}")
    except Exception as e:
        logger.error(f"Failed to save descriptor report to {output_path}: {e}")


def calculate_similarity_matrix(smiles_list: List[str], 
                                fingerprint_type: str = 'morgan', 
                                radius: int = 2, 
                                n_bits: int = 2048) -> np.ndarray:
    """
    Calculates pairwise Tanimoto similarity using Morgan fingerprints.
    
    Args:
        smiles_list (List[str]): List of SMILES strings.
        fingerprint_type (str): Type of fingerprint to use (currently only 'morgan' is supported).
        radius (int): Morgan fingerprint radius.
        n_bits (int): Number of bits for the fingerprint.
        
    Returns:
        np.ndarray: N x N pairwise Tanimoto similarity matrix.
    """
    mols = [Chem.MolFromSmiles(s) for s in smiles_list]
    valid_mols = []
    valid_indices = []
    
    for i, m in enumerate(mols):
        if m is not None:
            valid_mols.append(m)
            valid_indices.append(i)
        else:
            logger.warning(f"Invalid SMILES at index {i}: {smiles_list[i]}")
            
    n_mols = len(smiles_list)
    sim_matrix = np.zeros((n_mols, n_mols))
    
    if fingerprint_type.lower() == 'morgan':
        fps = [AllChem.GetMorganFingerprintAsBitVect(m, radius, nBits=n_bits) for m in valid_mols]
    else:
        logger.error(f"Unsupported fingerprint type: {fingerprint_type}")
        return sim_matrix
        
    # Calculate pairwise similarities for valid molecules
    for i in range(len(valid_mols)):
        for j in range(i, len(valid_mols)):
            sim = DataStructs.TanimotoSimilarity(fps[i], fps[j])
            orig_i = valid_indices[i]
            orig_j = valid_indices[j]
            sim_matrix[orig_i, orig_j] = sim
            sim_matrix[orig_j, orig_i] = sim
            
    return sim_matrix
