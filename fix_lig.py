import re

with open('src/ligand.py', 'r') as f:
    text = f.read()

new_func = """def standardize_mol(mol: Chem.Mol) -> Optional[Chem.Mol]:
    if mol is None: return None
    try:
        from rdkit.Chem.MolStandardize import rdMolStandardize
        mol = Chem.RemoveHs(mol)
        mol = rdMolStandardize.MetalDisconnector().Disconnect(mol)
        mol = rdMolStandardize.LargestFragmentChooser().choose(mol)
        mol = rdMolStandardize.Uncharger().uncharge(mol)
        return Chem.AddHs(mol)
    except Exception as e:
        return Chem.AddHs(mol)
"""

text = re.sub(r'def standardize_mol.*?return Chem\.AddHs\(mol\)', new_func, text, flags=re.DOTALL)
with open('src/ligand.py', 'w', encoding='utf-8') as f:
    f.write(text)
