import yaml

with open('config/docking_config.yaml', 'r') as f:
    config = yaml.safe_load(f)

config['docking']['vina_executable'] = 'bin/vina.exe'

with open('config/docking_config.yaml', 'w') as f:
    yaml.dump(config, f, default_flow_style=False)
