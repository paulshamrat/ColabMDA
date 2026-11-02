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


def test_tleap_script_with_solvated_pdb(tmp_path):
    spec = SystemSpec.from_tag("atp_2mg")
    out_pdb = tmp_path / "solvated.pdb"
    script = generate_tleap_script(
        spec=spec,
        input_pdb=tmp_path / "model.pdb",
        output_prmtop=tmp_path / "system.prmtop",
        output_inpcrd=tmp_path / "system.inpcrd",
        output_pdb=out_pdb,
    )
    assert f'savepdb mol "{out_pdb.resolve()}"' in script


def test_write_prep_provenance(tmp_path):
    import json

    from colabmda.amber.leap import write_prep_provenance

    spec = SystemSpec.from_tag("atp_2mg", variant="WT")
    dest_json = tmp_path / "prep.json"
    written = write_prep_provenance(spec, dest_json)

    assert written.exists()
    data = json.loads(written.read_text(encoding="utf-8"))
    assert data["spec"]["state"] == "atp_2mg"
    assert data["spec"]["variant"] == "WT"
    assert data["parameters"]["protein_ff"] == "ff19SB"
    assert data["parameters"]["water_model"] == "opc"
    assert data["provenance"]["engine"] == "ColabMDA"


def test_stripped_topology_is_written_by_strip_action(tmp_path):
    script = generate_cpptraj_strip_script(
        tmp_path / "input.prmtop",
        [tmp_path / "frame.dcd"],
        tmp_path / "stripped.dcd",
        output_top=tmp_path / "stripped.prmtop",
    )
    assert 'parmout "' in script
    assert "parmwrite" not in script


def test_opc_preparation_uses_one_complete_ion_potential(tmp_path):
    import shutil

    import pytest

    openmm = pytest.importorskip("openmm")
    from openmm import unit
    from openmm.app import AmberPrmtopFile

    from colabmda.amber.leap import run_tleap

    tleap = shutil.which("tleap")
    if not tleap:
        pytest.skip("AmberTools is required for the ion topology integration test")
    pdb = tmp_path / "mg.pdb"
    pdb.write_text(
        "HETATM    1 MG    MG A   1       0.000   0.000   0.000  1.00  0.00          MG\nEND\n"
    )
    top = tmp_path / "system.prmtop"
    text = generate_tleap_script(
        SystemSpec(ions={"MG": 1}),
        pdb,
        top,
        tmp_path / "system.inpcrd",
        padding_angstroms=8,
        salt_conc_molar=0,
    )
    status, log = run_tleap(text, tmp_path, tleap_binary=tleap)
    assert status == 0 and top.exists(), log
    amber = AmberPrmtopFile(str(top))
    system = amber.createSystem()
    force = next(f for f in system.getForces() if isinstance(f, openmm.NonbondedForce))
    mg = next(a for a in amber.topology.atoms() if a.residue.name == "MG")
    charge, sigma, epsilon = force.getParticleParameters(mg.index)
    # Complete OPC 12-6 Mg values; the former 12-6-4 override gives 1.405/0.01652939.
    assert charge.value_in_unit(unit.elementary_charge) == pytest.approx(2)
    assert sigma.value_in_unit(unit.angstrom) * 2 ** (1 / 6) / 2 == pytest.approx(1.330)
    assert epsilon.value_in_unit(unit.kilocalorie_per_mole) == pytest.approx(0.00716930)


def test_tleap_script_multi_component(tmp_path):
    spec = SystemSpec.from_tag("atpmg")
    script = generate_tleap_script(
        spec=spec,
        input_pdb=tmp_path / "model.pdb",
        output_prmtop=tmp_path / "system.prmtop",
        output_inpcrd=tmp_path / "system.inpcrd",
        component_pdbs={
            "REC": tmp_path / "WT.pdb",
            "LIG": tmp_path / "atp.pdb",
            "MG": tmp_path / "mg.pdb",
        },
    )
    assert 'REC = loadpdb "' in script
    assert 'LIG = loadpdb "' in script
    assert 'MG = loadpdb "' in script
    assert "mol = combine { REC LIG MG }" in script


def test_apply_hu2024_cmap_safe_on_non_nucleotide(tmp_path):
    import pytest

    from colabmda.amber.cmap import apply_hu2024_cmap

    with pytest.raises(FileNotFoundError):
        apply_hu2024_cmap(tmp_path / "nonexistent.prmtop")
