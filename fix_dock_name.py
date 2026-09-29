import re

with open('src/docking.py', 'r') as f:
    text = f.read()

new_ligand_file = """        safe_name = str(comp_name).replace(" ", "_").replace("/", "_").replace("\\\\", "_")
        if not safe_name: safe_name = str(comp_id)
        ligand_file = ligand_dir / f"{safe_name}.pdbqt\"""

text = re.sub(r'ligand_file = ligand_dir / f"\{comp_id\}\.pdbqt"', new_ligand_file, text)

with open('src/docking.py', 'w', encoding='utf-8') as f:
    f.write(text)
