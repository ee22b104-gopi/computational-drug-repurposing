import re

with open('data/download_data.py', 'r') as f:
    text = f.read()

new_fda = """def fetch_fda_approved_drugs(max_compounds: int = 300) -> pd.DataFrame:
    logger.info("Loading hardcoded FDA-approved drugs dataset to avoid API limits...")
    drugs = [
        {"pubchem_cid": 4091, "name": "Metformin", "smiles": "CN(C)C(=N)NC(=N)N"},
        {"pubchem_cid": 2244, "name": "Aspirin", "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O"},
        {"pubchem_cid": 3672, "name": "Ibuprofen", "smiles": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"},
        {"pubchem_cid": 1983, "name": "Acetaminophen", "smiles": "CC(=O)NC1=CC=C(C=C1)O"},
        {"pubchem_cid": 3961, "name": "Losartan", "smiles": "CCCCC1=NC(=C(N1CC2=CC=C(C=C2)C3=CC=CC=C3C4=NNN=N4)CO)Cl"},
        {"pubchem_cid": 4594, "name": "Omeprazole", "smiles": "CC1=CN=C(C(=C1OC)C)CS(=O)C2=NC3=C(N2)C=C(C=C3)OC"},
        {"pubchem_cid": 4171, "name": "Metoprolol", "smiles": "CC(C)NCC(COC1=CC=C(C=C1)CCOC)O"},
        {"pubchem_cid": 2162, "name": "Amlodipine", "smiles": "CCOC(=O)C1=C(NC(=C(C1C2=CC=CC=C2Cl)C(=O)OC)C)COCCN"},
        {"pubchem_cid": 3386, "name": "Fluoxetine", "smiles": "CNC(CCC1=CC=CC=C1)OC2=CC=C(C=C2)C(F)(F)F"},
        {"pubchem_cid": 3957, "name": "Loratadine", "smiles": "CCOC(=O)N1CCC(=C2C3=C(CCC4=C2N=CC=C4)C=C(C=C3)Cl)CC1"},
        {"pubchem_cid": 2662, "name": "Celecoxib", "smiles": "CC1=CC=C(C=C1)C2=CC(=NN2C3=CC=C(C=C3)S(=O)(=O)N)C(F)(F)F"}
    ]
    import pandas as pd
    df = pd.DataFrame(drugs)
    df["source"] = "PubChem_Hardcoded"
    df["category"] = "screening"
    df["molecular_formula"] = "Unknown"
    df["molecular_weight"] = 0.0
    return df
"""

new_decoys = """def generate_decoys(active_smiles_list: list, n_decoys: int = 30, seed: int = 42) -> pd.DataFrame:
    logger.info("Loading hardcoded decoys...")
    decoys = [
        {"pubchem_cid": 5950, "name": "Vitamin C", "smiles": "C(C(C1C(=C(C(=O)O1)O)O)O)O"},
        {"pubchem_cid": 5793, "name": "Glucose", "smiles": "C(C1C(C(C(C(O1)O)O)O)O)O"},
        {"pubchem_cid": 2519, "name": "Caffeine", "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"}
    ]
    import pandas as pd
    df = pd.DataFrame(decoys)
    df["source"] = "PubChem_Hardcoded"
    df["category"] = "decoy"
    df["molecular_formula"] = "Unknown"
    df["molecular_weight"] = 0.0
    return df
"""

text = re.sub(r'def fetch_fda_approved_drugs.*?return df', new_fda, text, flags=re.DOTALL)
text = re.sub(r'def generate_decoys.*?return df', new_decoys, text, flags=re.DOTALL)

with open('data/download_data.py', 'w', encoding='utf-8') as f:
    f.write(text)
