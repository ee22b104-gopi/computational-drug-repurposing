"""
Computational Drug Repurposing Pipeline
========================================

Structure-based virtual screening for drug repurposing targeting
human Dipeptidyl Peptidase-4 (DPP-4).

Modules
-------
protein : Protein structure download, cleaning, and PDBQT preparation
ligand : Ligand standardization, 3D coordinate generation, PDBQT conversion
docking : AutoDock Vina execution and result parsing
descriptors : RDKit-based molecular descriptor calculation
interactions : Protein-ligand interaction fingerprint analysis
validation : Redocking validation and enrichment metrics
robustness : Docking parameter sensitivity analysis
visualization : Publication-quality figure generation
"""

__version__ = "1.0.0"
__author__ = "Computational Drug Discovery Project"
