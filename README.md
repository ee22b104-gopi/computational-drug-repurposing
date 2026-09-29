# Structure-Based Virtual Screening for Drug Repurposing

A reproducible computational pipeline for structure-based virtual screening
targeting human Dipeptidyl Peptidase-4 (DPP-4), built with Python, AutoDock Vina,
and RDKit.

## Overview

This project implements a complete drug repurposing workflow:

1. **Target Validation** — DPP-4 as a validated Type 2 Diabetes target
2. **Structure Selection** — PDB 1X70 (2.10 Å, sitagliptin-bound)
3. **Binding Site Identification** — From co-crystallized ligand coordinates
4. **Compound Library** — ~250–350 compounds (approved drugs + positive controls + decoys)
5. **Molecular Docking** — AutoDock Vina with systematic validation
6. **Chemical Analysis** — RDKit-based property calculation and Lipinski filtering
7. **Interaction Analysis** — Protein–ligand contact analysis
8. **Robustness Testing** — Parameter sensitivity and ranking stability
9. **Candidate Prioritization** — Multi-criteria selection with scientific caution

> **⚠️ Important:** Molecular docking is a computational prediction technique.
> Favorable docking scores do not prove therapeutic efficacy, safety, or selectivity.
> All identified candidates are "computationally prioritized" and require
> experimental validation.

## Quick Start

### Prerequisites

- Python 3.10+
- [AutoDock Vina](https://github.com/ccsb-scripps/AutoDock-Vina/releases) binary (`vina.exe` on Windows)

### Installation

```bash
# Clone the repository
git clone https://github.com/ee22b104-gopi/computational-drug-repurposing.git
cd computational-drug-repurposing

# Install dependencies
pip install -r requirements.txt

# Or use conda
conda env create -f environment.yml
conda activate drug-repurposing

# Download Vina (Windows)
# Place vina.exe in your PATH or project root
```

### Running the Pipeline

```bash
# Complete pipeline (recommended)
python run_pipeline.py

# Or run stages individually:
python data/download_data.py          # Download structures and compounds
python scripts/prepare_protein.py     # Prepare receptor
python scripts/prepare_ligands.py     # Prepare ligand library
python scripts/run_docking.py         # Run docking + validation
python scripts/analyze_results.py     # Chemical and interaction analysis
python scripts/generate_figures.py    # Generate publication figures
```

## Project Structure

```
computational-drug-repurposing/
├── README.md                  # This file
├── LICENSE                    # MIT License
├── CITATION.cff               # Citation metadata
├── requirements.txt           # Pip dependencies
├── environment.yml            # Conda environment
├── run_pipeline.py            # Master automation script
│
├── config/
│   └── docking_config.yaml    # All parameters (single source of truth)
│
├── data/
│   ├── README.md              # Data provenance documentation
│   ├── download_data.py       # Programmatic data retrieval
│   └── protein_metadata.json  # Structure selection metadata
│
├── src/
│   ├── protein.py             # Protein download, cleaning, PDBQT prep
│   ├── ligand.py              # Ligand standardization, 3D generation
│   ├── docking.py             # AutoDock Vina execution and parsing
│   ├── descriptors.py         # RDKit molecular descriptors
│   ├── interactions.py        # Protein–ligand interaction analysis
│   ├── validation.py          # Redocking RMSD and enrichment
│   ├── robustness.py          # Parameter sensitivity analysis
│   └── visualization.py       # Publication figure generation
│
├── scripts/
│   ├── prepare_protein.py     # Entry point: protein preparation
│   ├── prepare_ligands.py     # Entry point: ligand preparation
│   ├── run_docking.py         # Entry point: docking execution
│   ├── analyze_results.py     # Entry point: results analysis
│   └── generate_figures.py    # Entry point: figure generation
│
├── results/
│   ├── figures/               # Generated figures (Fig 1–9)
│   ├── tables/                # Summary tables
│   ├── docking/               # Raw docking output
│   ├── docking_results.csv    # All docking scores
│   ├── compound_descriptors.csv
│   ├── interactions.csv
│   ├── ranking_robustness.csv
│   └── summary.csv
│
├── paper/
│   ├── manuscript.tex         # Research manuscript
│   ├── references.bib         # Verified BibTeX references
│   └── figures/               # Publication figures
│
└── tests/
    └── test_pipeline.py       # Unit tests
```

## Target: DPP-4

| Property | Value |
|---|---|
| Protein | Dipeptidyl Peptidase-4 (DPP-4 / CD26) |
| UniProt | [P27487](https://www.uniprot.org/uniprotkb/P27487) |
| EC | 3.4.14.5 |
| PDB Structure | [1X70](https://www.rcsb.org/structure/1X70) (2.10 Å, sitagliptin) |
| Binding Site | Ser630–Asp708–His740 catalytic triad + S1/S2/S2-ext subsites |

## Methods

- **Docking Engine:** AutoDock Vina 1.2+ (Trott & Olson, 2010; Eberhardt et al., 2021)
- **Cheminformatics:** RDKit (Landrum, 2006)
- **Ligand Preparation:** Meeko (Forli Lab, Scripps)
- **Interaction Analysis:** Distance-based + ProLIF (optional)
- **Validation:** Redocking RMSD + enrichment factor analysis

## Key References

1. Kim D et al. *J Med Chem.* 2005;48(1):141-151. [DOI: 10.1021/jm0493156](https://doi.org/10.1021/jm0493156)
2. Trott O, Olson AJ. *J Comput Chem.* 2010;31(2):455-461. [DOI: 10.1002/jcc.21334](https://doi.org/10.1002/jcc.21334)
3. Nabeno M et al. *BBRC.* 2013;434(2):191-196. [DOI: 10.1016/j.bbrc.2013.03.010](https://doi.org/10.1016/j.bbrc.2013.03.010)

## License

MIT License. See [LICENSE](LICENSE).

## Disclaimer

This project is for educational and research purposes. Computational predictions
do not constitute drug discovery claims. No compound identified through this
pipeline should be considered a therapeutic candidate without extensive
experimental validation including biochemical assays, cell-based studies,
pharmacokinetic profiling, and clinical trials.
