from colabmda.amber.cpptraj import generate_cpptraj_strip_script
from colabmda.amber.leap import generate_tleap_script, get_param_path
from colabmda.core.taxonomy import SystemSpec


def test_bundled_parameters_exist():
    assert get_param_path("ATP-B3.prepi").exists()
    assert get_param_path("ATP-B3.frcmod").exists()
    assert get_param_path("ADP-B3.prepi").exists()
    assert get_param_path("GDP-B3.prepi").exists()
    assert get_param_path("GTP-B3.prepi").exists()


def test_tleap_script_apo(tmp_path):
    spec = SystemSpec.from_tag("apo")
    script = generate_tleap_script(
        spec=spec,
        input_pdb=tmp_path / "model.pdb",
        output_prmtop=tmp_path / "system.prmtop",
        output_inpcrd=tmp_path / "system.inpcrd",
    )
    assert "source leaprc.protein.ff19SB" in script
    assert "source leaprc.water.opc" in script
    assert "ATP-B3" not in script
    assert "phosaa19SB" not in script
    assert "solvateBox mol OPCBOX 14.0" in script


def test_tleap_script_atp_2mg(tmp_path):
    spec = SystemSpec.from_tag("atp_2mg")
    script = generate_tleap_script(
        spec=spec,
        input_pdb=tmp_path / "model.pdb",
        output_prmtop=tmp_path / "system.prmtop",
        output_inpcrd=tmp_path / "system.inpcrd",
    )
    assert "addAtomTypes" in script
    assert "ATP-B3.prepi" in script
    assert "ATP-B3.frcmod" in script
    assert "phosaa19SB" not in script


def test_tleap_script_active_phosphorylated_substrate(tmp_path):
    spec = SystemSpec.from_tag("atp_2mg_eb2")
    script = generate_tleap_script(
        spec=spec,
        input_pdb=tmp_path / "model.pdb",
        output_prmtop=tmp_path / "system.prmtop",
        output_inpcrd=tmp_path / "system.inpcrd",
    )
    assert "source leaprc.phosaa19SB" in script  # Loaded for pY171 / PTR
    assert "ATP-B3.prepi" in script
    assert "ATP-B3.frcmod" in script


def test_cpptraj_strip_script(tmp_path):
    top = tmp_path / "system.prmtop"
    trajs = [tmp_path / "chk_001.dcd", tmp_path / "chk_002.dcd"]

    # Test keep_ligand stream
    s_lig = generate_cpptraj_strip_script(
        topology=top,
        trajectories=trajs,
        output_traj=tmp_path / "prod_keep_ligand.dcd",
        mode="keep_ligand",
    )
    assert "autoimage" in s_lig
    assert "strip :WAT,Na+,Cl-" in s_lig
    assert "prod_keep_ligand.dcd" in s_lig

    # Test ca_backbone stream
    s_ca = generate_cpptraj_strip_script(
        topology=top,
        trajectories=trajs,
        output_traj=tmp_path / "prod_ca.dcd",
        mode="ca_backbone",
    )
    assert "strip !@CA" in s_ca
