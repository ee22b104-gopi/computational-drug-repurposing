# Computational Drug Repurposing Pipeline: Progress Summary

This document summarizes the current state of the computational drug repurposing pipeline for the DPP-4 target.

## Project Structure & Setup
The repository was initialized successfully at `c:\Users\gaadi\OneDrive\Desktop\PLACEMENTS\Computational Drug Discovery\computational-drug-repurposing`.

We have structured the pipeline into independent stages:
1. **Data Download** (`download_data.py`)
2. **Protein Preparation** (`prepare_protein.py`)
3. **Ligand Preparation** (`prepare_ligands.py`)
4. **Docking Execution** (`run_docking.py`)
5. **Results Analysis & Figure Generation** (`analyze_results.py`, `generate_figures.py`)

## What has been accomplished so far:

### 1. Data Download (`data/download_data.py`)
- **Action**: Ran the script to compile our compound screening library.
- **Result**: Due to severe rate limiting on the PubChem API servers (returning 503 Server Busy errors), we bypassed the API by hardcoding a diverse list of 25+ FDA-approved drugs (including controls like Sitagliptin, Vildagliptin, and non-target drugs like Aspirin, Losartan, Metformin). 
- **Files Created**: 
  - `data/compounds.csv`: The finalized list of 24 screening compounds, ready for processing.

### 2. Protein Preparation (`scripts/prepare_protein.py`)
- **Action**: Downloaded the PDB structure 1X70 (DPP-4 bound to sitagliptin), extracted chain A, cleaned it (removed waters/heteroatoms), extracted the co-crystallized ligand `715` to find the exact binding pocket coordinates, and prepared the receptor.
- **Files Created**:
  - `data/raw/1X70.pdb` (Original crystal structure)
  - `data/processed/1X70_chainA.pdb` (Extracted chain A)
  - `data/processed/1X70_clean.pdb` (Cleaned protein)
  - `data/processed/1X70_ligand_715.pdb` (Co-crystallized ligand reference)
  - `data/processed/1X70_receptor.pdbqt` (Final prepared receptor, ready for docking)
  - `results/protein_preparation_report.json` (Summary log)
- **Impact**: We calculated the exact center of the binding site `(41.184, 51.045, 35.629)` using the native ligand, allowing for highly targeted molecular docking.

### 3. Ligand Preparation (`scripts/prepare_ligands.py`)
- **Action**: Loaded `compounds.csv`, built 3D conformers for all molecules using RDKit, standardized their structures, assigned Gasteiger charges, and exported them to AutoDock Vina's PDBQT format. 
- **Files Created**:
  - `data/processed/ligands/*.pdbqt` (24 individual 3D prepared structures, e.g., `Aspirin.pdbqt`, `Sitagliptin.pdbqt`)
  - `results/compound_descriptors.csv` (Computed molecular descriptors for the library)
  - `results/ligand_preparation_report.json` (Summary log)

### 4. Docking Execution (`scripts/run_docking.py`)
- **Action**: The main AutoDock Vina pipeline has been successfully configured and launched.
- **Status**: **Currently running in the background.** It iterates through all 24 prepared ligands, testing their binding affinities against the DPP-4 active site we identified.
- **Files Being Created**:
  - `results/docking/*_out.pdbqt` (The 3D binding poses for each compound)
  - `results/docking_results.csv` (The final table compiling binding affinities and scores)

## Next Steps
Once the docking job finishes in the background, the pipeline will have all the raw data needed. The final commands to run will be:
```bash
# 1. Analyze the docking results and filter candidates
python scripts/analyze_results.py

# 2. Generate visualization figures for your portfolio/paper
python scripts/generate_figures.py
```
