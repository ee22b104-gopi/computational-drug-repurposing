import re

with open('scripts/prepare_ligands.py', 'r') as f:
    text = f.read()

import_pd = "import pandas as pd\nimport yaml"
text = text.replace("import yaml", import_pd)

call_1 = """compounds_df = pd.read_csv(compounds_csv)
        library_results = prepare_ligand_library(
            compounds_df,
            str(data_dir / "processed" / "ligands")
        )"""
text = re.sub(r'library_results = prepare_ligand_library\(\s*str\(compounds_csv\),\s*str\(data_dir / "processed" / "ligands"\)\s*\)', call_1, text, flags=re.DOTALL)

call_2 = """descriptors_df = calculate_library_descriptors(compounds_df)"""
text = text.replace("descriptors_df = calculate_library_descriptors(str(compounds_csv))", call_2)

text = text.replace('ligand_pdb_path = data_dir / "processed" / "ligand_715.pdb"', 'ligand_pdb_path = data_dir / "processed" / "1X70_ligand_715.pdb"')

call_ref = """ref_pdbqt = prepare_reference_ligand(
                str(ligand_pdb_path),
                "715", "A",
                str(data_dir / "processed")
            )"""
text = re.sub(r'ref_pdbqt = prepare_reference_ligand\([^)]+\)', call_ref, text)

with open('scripts/prepare_ligands.py', 'w', encoding='utf-8') as f:
    f.write(text)
