import re

with open('scripts/run_docking.py', 'r') as f:
    text = f.read()

text = text.replace('if not check_vina_installation():', 'if not check_vina_installation(vina_executable=config.get("docking", {}).get("vina_executable", "bin/vina.exe")):')

with open('scripts/run_docking.py', 'w', encoding='utf-8') as f:
    f.write(text)
