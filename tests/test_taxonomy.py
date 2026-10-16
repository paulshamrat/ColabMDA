import pytest

from colabmda.core.ions import get_ion, parse_ion_spec
from colabmda.core.project import ProjectPaths
from colabmda.core.taxonomy import SUBSTRATE_PRESETS, SystemSpec


class TestIonRegistry:
    def test_known_ions(self):
        mg = get_ion("MG")
        assert mg.element == "Mg"
        assert mg.charge == 2
        assert mg.amber_unit == "Mg+"
        assert mg.openmm_name == "MG"

        mn = get_ion("mn")  # Case insensitive
        assert mn.element == "Mn"
        assert mn.charge == 2
        assert mn.amber_unit == "Mn2+"

    def test_unknown_ion(self):
        with pytest.raises(ValueError, match="Unknown ion symbol 'XX'"):
            get_ion("XX")

    def test_parse_ion_spec(self):
        assert parse_ion_spec("MG:2") == {"MG": 2}
        assert parse_ion_spec("mg:1, mn:1") == {"MG": 1, "MN": 1}
        assert parse_ion_spec("ZN:1") == {"ZN": 1}
        assert parse_ion_spec("0") == {}
        assert parse_ion_spec("none") == {}
        assert parse_ion_spec(None) == {}
        assert parse_ion_spec({"MG": 2, "MN": 0}) == {"MG": 2}


class TestSystemSpec:
    def test_apo_state(self):
        spec = SystemSpec.from_tag("apo")
        assert spec.is_apo is True
        assert spec.ligand is None
        assert spec.ions == {}
        assert spec.state_tag == "apo"
        assert spec.is_phosphorylated is False

    def test_atp_1mg_state(self):
        spec = SystemSpec.from_tag("atp_1mg")
        assert spec.is_apo is False
        assert spec.ligand == "ATP"
        assert spec.ions == {"MG": 1}
        assert spec.state_tag == "atp_1mg"

        # atp_mg alias
        spec_alias = SystemSpec.from_tag("atp_mg")
        assert spec_alias.ligand == "ATP"
        assert spec_alias.ions == {"MG": 1}

    def test_atp_2mg_state(self):
        spec = SystemSpec.from_tag("atp_2mg")
        assert spec.ligand == "ATP"
        assert spec.ions == {"MG": 2}
        assert spec.state_tag == "atp_2mg"

    def test_adp_2mg_py171_state(self):
        spec = SystemSpec.from_tag("adp_2mg_py171")
        assert spec.ligand == "ADP"
        assert spec.ions == {"MG": 2}
        assert spec.phos_residues == ["Y171"]
        assert spec.is_phosphorylated is True

    def test_adp_1mg_state(self):
        spec = SystemSpec.from_tag("adp_1mg")
        assert spec.ligand == "ADP"
        assert spec.ions == {"MG": 1}
        assert spec.state_tag == "adp_1mg"

    def test_gdp_1mg_state(self):
        spec = SystemSpec.from_tag("gdp_1mg")
        assert spec.ligand == "GDP"
        assert spec.ions == {"MG": 1}
        assert spec.state_tag == "gdp_1mg"

    def test_active_substrate_complex(self):
        spec = SystemSpec.from_tag("atp_2mg_eb2", variant="R59W")
        assert spec.variant == "R59W"
        assert spec.ligand == "ATP"
        assert spec.ions == {"MG": 2}
        assert spec.phos_residues == ["Y171"]
        assert spec.is_phosphorylated is True
        assert spec.peptide_name == "eb2"
        assert spec.peptide_seq == SUBSTRATE_PRESETS["eb2"]["seq"]
        assert spec.state_tag == "atp_2mg_eb2"

    def test_manganese_substrate_complex(self):
        spec = SystemSpec.from_tag("atp_2mn_eb2")
        assert spec.ligand == "ATP"
        assert spec.ions == {"MN": 2}
        assert spec.peptide_name == "eb2"
        assert spec.state_tag == "atp_2mn_eb2"

    def test_declarative_custom_composition(self):
        spec = SystemSpec(
            protein="cdkl5",
            variant="L70P",
            ligand="ANP",
            ions={"MG": 1, "MN": 1},
            peptide_name="eb1",
        )
        assert spec.ligand == "ANP"
        assert spec.ions == {"MG": 1, "MN": 1}
        assert spec.peptide_seq == SUBSTRATE_PRESETS["eb1"]["seq"]
        assert spec.state_tag == "anp_1mg_1mn_eb1"


class TestProjectPaths:
    def test_path_resolution(self, tmp_path):
        paths = ProjectPaths(root=tmp_path)
        spec = SystemSpec.from_tag("atp_2mg_eb2", variant="WT")

        assert paths.structure_dir == tmp_path / "structure"
        assert paths.template_dir == tmp_path / "structure" / "template"
        assert paths.scripts_dir == tmp_path / "scripts"
        assert paths.variant_dir(spec) == tmp_path / "simulation" / "atp_2mg_eb2" / "WT"
        assert paths.prep_dir(spec) == tmp_path / "simulation" / "atp_2mg_eb2" / "WT" / "prep"
        assert (
            paths.model_pdb(spec)
            == tmp_path / "simulation" / "atp_2mg_eb2" / "WT" / "prep" / "model.pdb"
        )
        assert (
            paths.prepared_prmtop(spec)
            == tmp_path / "simulation" / "atp_2mg_eb2" / "WT" / "prep" / "system.prmtop"
        )
        assert paths.eq_dir(spec) == tmp_path / "simulation" / "atp_2mg_eb2" / "WT" / "eq"
        assert (
            paths.replica_dir(spec, "r1") == tmp_path / "simulation" / "atp_2mg_eb2" / "WT" / "r1"
        )
        assert (
            paths.analysis_dir(spec) == tmp_path / "simulation" / "atp_2mg_eb2" / "WT" / "analysis"
        )

    def test_init_layout(self, tmp_path):
        paths = ProjectPaths(root=tmp_path)
        paths.init_layout()
        assert (tmp_path / "scripts").exists()
        assert (tmp_path / "structure" / "template").exists()
        assert (tmp_path / "structure" / "protein").exists()
        assert (tmp_path / "simulation").exists()
        assert (tmp_path / "01_structure" / "protein").exists()
        assert (tmp_path / "02_simulation").exists()
        assert (tmp_path / "03_analysis").exists()

    def test_numbered_path_resolution(self, tmp_path):
        numbered_root = tmp_path / "numbered_project"
        paths = ProjectPaths(root=numbered_root)
        spec = SystemSpec.from_tag("atp_2mg_eb2", variant="WT")

        # Create pure numbered hierarchy
        (numbered_root / "01_structure" / "variants" / "WT").mkdir(parents=True)
        (numbered_root / "02_simulation").mkdir(parents=True)
        (numbered_root / "03_analysis").mkdir(parents=True)

        assert paths.structure_dir == numbered_root / "01_structure"
        assert paths.simulation_base == numbered_root / "02_simulation"
        assert paths.analysis_base == numbered_root / "03_analysis"
        assert paths.variant_dir(spec) == numbered_root / "02_simulation" / "atp_2mg_eb2" / "WT"
        assert (
            paths.prep_dir(spec)
            == numbered_root / "02_simulation" / "atp_2mg_eb2" / "WT" / "01_prep"
        )
        assert (
            paths.eq_dir(spec) == numbered_root / "02_simulation" / "atp_2mg_eb2" / "WT" / "02_eq"
        )
        assert (
            paths.prepared_solvated(spec)
            == numbered_root / "02_simulation" / "atp_2mg_eb2" / "WT" / "01_prep" / "solvated.pdb"
        )
        assert (
            paths.prep_json(spec)
            == numbered_root / "02_simulation" / "atp_2mg_eb2" / "WT" / "01_prep" / "prep.json"
        )
        assert (
            paths.compare_dir("atp_2mg_eb2")
            == numbered_root / "03_analysis" / "atp_2mg_eb2" / "compare"
        )

    def test_clean_hybrid_layout(self, tmp_path):
        """Test unnumbered root (structure, simulation, analysis) with numbered sub-stages (01_prep, 02_eq, 03_prod)."""
        root = tmp_path / "hybrid_project"
        paths = ProjectPaths(root=root)
        spec = SystemSpec.from_tag("atp_2mg", variant="WT")

        (root / "structure" / "variants" / "WT").mkdir(parents=True)
        (root / "simulation" / "atp_2mg" / "WT" / "01_prep").mkdir(parents=True)
        (root / "simulation" / "atp_2mg" / "WT" / "02_eq").mkdir(parents=True)
        (root / "simulation" / "atp_2mg" / "WT" / "03_prod" / "r1").mkdir(parents=True)
        (root / "analysis" / "atp_2mg" / "WT").mkdir(parents=True)

        assert paths.structure_dir == root / "structure"
        assert paths.simulation_base == root / "simulation"
        assert paths.analysis_base == root / "analysis"
        assert paths.prep_dir(spec) == root / "simulation" / "atp_2mg" / "WT" / "01_prep"
        assert paths.eq_dir(spec) == root / "simulation" / "atp_2mg" / "WT" / "02_eq"
        assert paths.prod_dir(spec) == root / "simulation" / "atp_2mg" / "WT" / "03_prod"
        assert (
            paths.replica_dir(spec, "r1")
            == root / "simulation" / "atp_2mg" / "WT" / "03_prod" / "r1"
        )
        assert paths.analysis_dir(spec) == root / "analysis" / "atp_2mg" / "WT"
        assert (
            paths.prepared_solvated(spec)
            == root / "simulation" / "atp_2mg" / "WT" / "01_prep" / "solvated.pdb"
        )
        assert (
            paths.prep_json(spec)
            == root / "simulation" / "atp_2mg" / "WT" / "01_prep" / "prep.json"
        )
