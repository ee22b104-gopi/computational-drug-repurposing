"""
Basic tests for the computational drug repurposing pipeline.

Run with: python -m pytest tests/ -v
"""

import sys
from pathlib import Path

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))


class TestConfig:
    """Test configuration loading."""

    def test_config_exists(self):
        config_path = PROJECT_ROOT / "config" / "docking_config.yaml"
        assert config_path.exists(), f"Config not found: {config_path}"

    def test_config_loads(self):
        import yaml
        config_path = PROJECT_ROOT / "config" / "docking_config.yaml"
        with open(config_path) as f:
            config = yaml.safe_load(f)
        assert config is not None
        assert "protein" in config
        assert "docking" in config
        assert "library" in config

    def test_positive_controls_present(self):
        import yaml
        config_path = PROJECT_ROOT / "config" / "docking_config.yaml"
        with open(config_path) as f:
            config = yaml.safe_load(f)
        controls = config["library"]["positive_controls"]
        assert len(controls) >= 10, "Should have at least 10 positive controls"
        for ctrl in controls:
            assert "name" in ctrl
            assert "smiles" in ctrl
            assert "pubchem_cid" in ctrl

    def test_pdb_id_valid(self):
        import yaml
        config_path = PROJECT_ROOT / "config" / "docking_config.yaml"
        with open(config_path) as f:
            config = yaml.safe_load(f)
        pdb_id = config["protein"]["pdb_id"]
        assert len(pdb_id) == 4, "PDB ID should be 4 characters"
        assert pdb_id == "1X70"


class TestDescriptors:
    """Test molecular descriptor calculation."""

    def test_import(self):
        from src.descriptors import calculate_descriptors
        assert callable(calculate_descriptors)

    def test_valid_smiles(self):
        from src.descriptors import calculate_descriptors
        # Aspirin
        result = calculate_descriptors("CC(=O)Oc1ccccc1C(=O)O", name="aspirin")
        assert result is not None
        assert result["molecular_weight"] > 0
        assert "logp" in result
        assert "hbd" in result
        assert "hba" in result

    def test_invalid_smiles(self):
        from src.descriptors import calculate_descriptors
        result = calculate_descriptors("INVALID_SMILES")
        assert result is None or result.get("molecular_weight") is None

    def test_sitagliptin_descriptors(self):
        from src.descriptors import calculate_descriptors
        smiles = "C1CN2C(=NN=C2C(F)(F)F)CN1C(=O)C[C@@H](CC3=CC(=C(C=C3F)F)F)N"
        result = calculate_descriptors(smiles, name="sitagliptin")
        assert result is not None
        # Sitagliptin MW ~407
        assert 390 < result["molecular_weight"] < 420

    def test_lipinski_check(self):
        from src.descriptors import calculate_descriptors, check_lipinski
        result = calculate_descriptors("CC(=O)Oc1ccccc1C(=O)O")
        lipinski = check_lipinski(result)
        assert "passes" in lipinski
        assert "n_violations" in lipinski
        # Aspirin should pass Lipinski
        assert lipinski["passes"] is True


class TestValidation:
    """Test validation utilities."""

    def test_import(self):
        from src.validation import calculate_rmsd
        assert callable(calculate_rmsd)

    def test_rmsd_identical(self):
        import numpy as np
        from src.validation import calculate_rmsd
        coords = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
        rmsd = calculate_rmsd(coords, coords)
        assert abs(rmsd) < 1e-10

    def test_rmsd_known(self):
        import numpy as np
        from src.validation import calculate_rmsd
        ref = np.array([[0.0, 0.0, 0.0]])
        pred = np.array([[1.0, 0.0, 0.0]])
        rmsd = calculate_rmsd(ref, pred)
        assert abs(rmsd - 1.0) < 1e-10

    def test_enrichment_calculation(self):
        import pandas as pd
        from src.validation import calculate_enrichment
        df = pd.DataFrame({
            "compound_id": [f"C{i}" for i in range(100)],
            "best_score": list(range(-10, -10 + 100)),
        })
        known_actives = ["C0", "C1", "C2"]  # Best scores
        result = calculate_enrichment(df, known_actives)
        assert result["total_compounds"] == 100
        assert result["total_actives"] == 3
        # All actives should be in top 5%
        assert result["enrichment_factors"]["EF_5.0%"]["actives_found"] == 3


class TestInteractions:
    """Test interaction analysis."""

    def test_import(self):
        from src.interactions import classify_interaction
        assert callable(classify_interaction)

    def test_hydrogen_bond_classification(self):
        from src.interactions import classify_interaction
        itype = classify_interaction("N", "N", "O", "O", 2.8)
        assert itype == "hydrogen_bond"

    def test_hydrophobic_classification(self):
        from src.interactions import classify_interaction
        itype = classify_interaction("CA", "C", "C1", "C", 3.5)
        assert itype == "hydrophobic"


class TestProteinMetadata:
    """Test protein metadata file."""

    def test_metadata_exists(self):
        metadata_path = PROJECT_ROOT / "data" / "protein_metadata.json"
        assert metadata_path.exists()

    def test_metadata_valid(self):
        import json
        metadata_path = PROJECT_ROOT / "data" / "protein_metadata.json"
        with open(metadata_path) as f:
            metadata = json.load(f)
        assert metadata["primary_structure"]["pdb_id"] == "1X70"
        assert metadata["primary_structure"]["chain"] == "A"
        assert metadata["primary_structure"]["ligand"]["name"] == "Sitagliptin"
        assert metadata["primary_structure"]["uniprot_id"] == "P27487"


if __name__ == "__main__":
    # Simple test runner without pytest
    import traceback

    test_classes = [TestConfig, TestDescriptors, TestValidation,
                    TestInteractions, TestProteinMetadata]
    passed = 0
    failed = 0
    errors = []

    for test_class in test_classes:
        instance = test_class()
        for method_name in dir(instance):
            if method_name.startswith("test_"):
                try:
                    getattr(instance, method_name)()
                    print(f"  PASS {test_class.__name__}.{method_name}")
                    passed += 1
                except Exception as e:
                    print(f"  FAIL {test_class.__name__}.{method_name}: {e}")
                    failed += 1
                    errors.append((f"{test_class.__name__}.{method_name}", str(e)))

    print(f"\n{passed} passed, {failed} failed")
    if errors:
        for name, err in errors:
            print(f"  FAIL: {name} — {err}")
