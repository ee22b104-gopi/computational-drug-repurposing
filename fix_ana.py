import re

with open('scripts/analyze_results.py', 'r') as f:
    text = f.read()

# Fix status case
text = text.replace('df[df["status"] == "success"]', 'df[df["status"].str.lower() == "success"]')

# Fix receptor path
new_receptor = """receptor_pdb = PROJECT_ROOT / "data" / "processed" / f"{config['protein']['pdb_id']}_clean.pdb\"""
text = re.sub(r'receptor_pdb = PROJECT_ROOT.*?pdb_id.*?_clean\.pdb"', new_receptor, text)

with open('scripts/analyze_results.py', 'w', encoding='utf-8') as f:
    f.write(text)
