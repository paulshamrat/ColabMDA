import numpy as np
import pytest

from colabmda.core.taxonomy import SystemSpec
from colabmda.modeller.substrate import (
    assemble_system_model,
    build_ptr_phosphate_atoms,
    dock_substrate_peptide,
    kabsch,
    parse_pdb_atoms,
    phosphorylate_residue,
    position_ions,
)


def test_kabsch_identity():
    pts = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    rot, trans = kabsch(pts, pts)
    np.testing.assert_allclose(rot, np.eye(3), atol=1e-6)
    np.testing.assert_allclose(trans, np.zeros(3), atol=1e-6)


def test_kabsch_translation_and_rotation():
    pts = np.array(
        [
            [1.0, 2.0, 3.0],
            [4.0, 5.0, 6.0],
            [7.0, 8.0, 10.0],
            [2.0, 0.0, 1.0],
        ]
    )
    theta = np.pi / 2
    r_true = np.array(
        [
            [np.cos(theta), -np.sin(theta), 0],
            [np.sin(theta), np.cos(theta), 0],
            [0, 0, 1],
        ]
    )
    t_true = np.array([5.0, -3.0, 2.0])

    target = pts @ r_true + t_true
    rot, trans = kabsch(pts, target)

    reconstructed = pts @ rot + trans
    np.testing.assert_allclose(reconstructed, target, atol=1e-5)


def test_build_ptr_phosphate_atoms_geometry():
    cz = np.array([1.660, 8.446, 1.413])
    oh = np.array([0.702, 9.116, 2.067])
    cg = np.array([3.674, 6.928, 0.170])

    tyr_atoms = {"CZ": cz, "OH": oh, "CG": cg}
    p_coord, oxygens = build_ptr_phosphate_atoms(tyr_atoms)

    # Verify P-OH distance is ~1.60 A
    dist_p_oh = np.linalg.norm(p_coord - oh)
    np.testing.assert_allclose(dist_p_oh, 1.60, atol=1e-3)

    # Verify tetrahedral P-O bond distances are ~1.51 A
    assert len(oxygens) == 3
    for oname in ["O1P", "O2P", "O3P"]:
        assert oname in oxygens
        dist_p_o = np.linalg.norm(oxygens[oname] - p_coord)
        np.testing.assert_allclose(dist_p_o, 1.51, atol=1e-3)


def test_phosphorylate_residue_synthetic(tmp_path):
    # Create minimal synthetic PDB with a TYR residue
    test_pdb = tmp_path / "tyr_input.pdb"
    lines = [
        "ATOM      1  N   TYR A 171       3.555   3.970  -0.000  1.00  0.00           N\n",
        "ATOM      2  CA  TYR A 171       4.853   4.614  -0.000  1.00  0.00           C\n",
        "ATOM      3  C   TYR A 171       5.608   4.507   1.328  1.00  0.00           C\n",
        "ATOM      4  O   TYR A 171       6.804   4.376   1.328  1.00  0.00           O\n",
        "ATOM      5  CB  TYR A 171       4.738   6.076  -0.493  1.00  0.00           C\n",
        "ATOM      6  CG  TYR A 171       3.674   6.928   0.170  1.00  0.00           C\n",
        "ATOM      7  CD1 TYR A 171       2.396   7.027  -0.377  1.00  0.00           C\n",
        "ATOM      8  CD2 TYR A 171       3.929   7.648   1.334  1.00  0.00           C\n",
        "ATOM      9  CE1 TYR A 171       1.394   7.764   0.230  1.00  0.00           C\n",
        "ATOM     10  CE2 TYR A 171       2.942   8.398   1.946  1.00  0.00           C\n",
        "ATOM     11  CZ  TYR A 171       1.660   8.446   1.413  1.00  0.00           C\n",
        "ATOM     12  OH  TYR A 171       0.702   9.116   2.067  1.00  0.00           O\n",
        "TER\n",
        "END\n",
    ]
    test_pdb.write_text("".join(lines))

    out_pdb = tmp_path / "ptr_output.pdb"
    phosphorylate_residue(
        input_pdb=test_pdb,
        output_pdb=out_pdb,
        res_num=171,
        chain="A",
        target_res="PTR",
    )

    ptr_atoms = parse_pdb_atoms(out_pdb, res_name="PTR", res_num=171)
    assert "P" in ptr_atoms
    assert "O1P" in ptr_atoms
    assert "O2P" in ptr_atoms
    assert "O3P" in ptr_atoms
    assert "OH" in ptr_atoms
    assert "CA" in ptr_atoms


@pytest.fixture
def mock_catalytic_system(tmp_path):
    """Create minimal synthetic PDBs for receptor and template with ATP, MG, and Chain S."""
    rec_lines = [
        "ATOM      1  N   ASP A 184     -20.000  -4.000  14.000  1.00  0.00           N\n",
        "ATOM      2  CA  ASP A 184     -21.000  -4.000  14.500  1.00  0.00           C\n",
        "ATOM     11  CZ  TYR A 171       1.660   8.446   1.413  1.00  0.00           C\n",
        "ATOM     12  OH  TYR A 171       0.702   9.116   2.067  1.00  0.00           O\n",
        "ATOM     13  CG  TYR A 171       3.674   6.928   0.170  1.00  0.00           C\n",
        "ATOM     14  CA  TYR A 171       4.853   4.614  -0.000  1.00  0.00           C\n",
        "HETATM 2974  PG  ATP A 300     -23.284  -1.633  18.344  1.00  6.41           P\n",
        "HETATM 2975  O1G ATP A 300     -23.539  -0.309  19.041  1.00  7.63           O\n",
        "HETATM 2976  O2G ATP A 300     -22.041  -2.359  18.835  1.00  6.62           O\n",
        "HETATM 2977  O3G ATP A 300     -23.373  -1.572  16.853  1.00  4.20           O\n",
        "HETATM 2978  PB  ATP A 300     -25.716  -3.327  17.941  1.00  7.58           P\n",
        "HETATM 2979  O1B ATP A 300     -25.238  -3.656  16.555  1.00  4.97           O\n",
        "HETATM 2980  O2B ATP A 300     -26.990  -2.516  18.088  1.00  8.68           O\n",
        "HETATM 2981  O3B ATP A 300     -24.518  -2.621  18.823  1.00 10.34           O\n",
        "HETATM 2982  PA  ATP A 300     -25.167  -6.025  18.854  1.00  4.75           P\n",
        "HETATM 2983  O1A ATP A 300     -25.417  -6.891  17.664  1.00  6.44           O\n",
        "HETATM 2984  O2A ATP A 300     -23.758  -5.701  19.278  1.00  5.63           O\n",
        "TER\n",
        "END\n",
    ]
    rec_pdb = tmp_path / "synthetic_receptor.pdb"
    rec_pdb.write_text("".join(rec_lines))

    tmpl_lines = [
        "HETATM 2972 MG    MG A 401     -23.565  -3.087  15.570  1.00  6.73          MG\n",
        "HETATM 2973 MG    MG A 402     -22.325  -4.320  18.980  1.00  4.80          MG\n",
        "HETATM 2974  PG  ATP A 403     -23.284  -1.633  18.344  1.00  6.41           P\n",
        "HETATM 2975  O1G ATP A 403     -23.539  -0.309  19.041  1.00  7.63           O\n",
        "HETATM 2976  O2G ATP A 403     -22.041  -2.359  18.835  1.00  6.62           O\n",
        "HETATM 2977  O3G ATP A 403     -23.373  -1.572  16.853  1.00  4.20           O\n",
        "HETATM 2978  PB  ATP A 403     -25.716  -3.327  17.941  1.00  7.58           P\n",
        "HETATM 2979  O1B ATP A 403     -25.238  -3.656  16.555  1.00  4.97           O\n",
        "HETATM 2980  O2B ATP A 403     -26.990  -2.516  18.088  1.00  8.68           O\n",
        "HETATM 2981  O3B ATP A 403     -24.518  -2.621  18.823  1.00 10.34           O\n",
        "HETATM 2982  PA  ATP A 403     -25.167  -6.025  18.854  1.00  4.75           P\n",
        "HETATM 2983  O1A ATP A 403     -25.417  -6.891  17.664  1.00  6.44           O\n",
        "HETATM 2984  O2A ATP A 403     -23.758  -5.701  19.278  1.00  5.63           O\n",
        "ATOM   2940  N   SER S 621     -22.000   1.500  17.000  1.00  5.00           N\n",
        "ATOM   2941  CA  SER S 621     -22.200   1.200  16.500  1.00  5.00           C\n",
        "ATOM   2942  CB  SER S 621     -22.400   0.800  16.200  1.00  5.00           C\n",
        "ATOM   2943  OG  SER S 621     -22.620   0.919  16.011  1.00  7.79           O\n",
        "ATOM   2944  C   SER S 621     -21.800   2.000  15.500  1.00  5.00           C\n",
        "ATOM   2945  O   SER S 621     -21.500   2.200  15.000  1.00  5.00           O\n",
        "TER\n",
        "END\n",
    ]
    tmpl_pdb = tmp_path / "synthetic_template.pdb"
    tmpl_pdb.write_text("".join(tmpl_lines))

    return rec_pdb, tmpl_pdb


def test_position_ions_synthetic(mock_catalytic_system, tmp_path):
    rec_pdb, tmpl_pdb = mock_catalytic_system
    out_pdb = tmp_path / "ions_out.pdb"

    placed = position_ions(
        receptor_pdb=rec_pdb,
        output_pdb=out_pdb,
        template_pdb=tmpl_pdb,
        ion_type="MG",
        count=2,
    )

    assert len(placed) == 2
    assert placed[0]["ion_type"] == "MG"
    assert placed[1]["ion_type"] == "MG"

    content = out_pdb.read_text()
    assert "HETATM" in content
    assert "MG" in content


def test_dock_substrate_peptide_synthetic(mock_catalytic_system, tmp_path):
    rec_pdb, tmpl_pdb = mock_catalytic_system
    out_pdb = tmp_path / "dock_out.pdb"

    qc = dock_substrate_peptide(
        receptor_pdb=rec_pdb,
        output_pdb=out_pdb,
        template_pdb=tmpl_pdb,
        peptide_seq="KRPSSAPSG",
        acceptor_idx=3,
        chain_id="B",
    )

    assert qc["status"] == "PASS"
    assert qc["distance_og_pg_A"] is not None
    # Distance in template is ~3.52 A
    assert 2.5 <= qc["distance_og_pg_A"] <= 5.0

    chain_b = parse_pdb_atoms(out_pdb, chain="B")
    assert len(chain_b) > 0


def test_assemble_system_model_synthetic(mock_catalytic_system, tmp_path):
    rec_pdb, tmpl_pdb = mock_catalytic_system
    out_pdb = tmp_path / "assembled_model.pdb"

    spec = SystemSpec.from_tag("atp_2mg_eb2", variant="WT")
    result = assemble_system_model(
        spec=spec,
        receptor_pdb=rec_pdb,
        output_pdb=out_pdb,
        template_pdb=tmpl_pdb,
    )

    assert result.exists()
    content = result.read_text()
    assert "ATP" in content
    assert "MG" in content
    assert " B " in content


def test_assemble_system_model_apo_normalization(tmp_path):
    # Starting structure has PTR171, ATP, MG, and Chain B peptide
    lines = [
        "ATOM      1  N   MET A   1       0.000   0.000   0.000  1.00  0.00           N\n",
        "ATOM      2  CA  MET A   1       1.000   0.000   0.000  1.00  0.00           C\n",
        "ATOM     11  CZ  PTR A 171       1.660   8.446   1.413  1.00  0.00           C\n",
        "ATOM     12  OH  PTR A 171       0.702   9.116   2.067  1.00  0.00           O\n",
        "ATOM     13  P   PTR A 171       0.000   9.000   2.000  1.00  0.00           P\n",
        "ATOM     14  O1P PTR A 171       0.100   9.100   2.100  1.00  0.00           O\n",
        "HETATM 2974  PG  ATP A 300     -23.284  -1.633  18.344  1.00  6.41           P\n",
        "HETATM 2975  MG   MG A 301     -22.000  -4.000  19.000  1.00  4.80          MG\n",
        "ATOM   3000  N   LYS B 219     -20.000   1.000  15.000  1.00  5.00           N\n",
        "ATOM   3001  CA  LYS B 219     -21.000   1.200  15.500  1.00  5.00           C\n",
        "TER\n",
        "END\n",
    ]
    inp_pdb = tmp_path / "complex_input.pdb"
    inp_pdb.write_text("".join(lines))

    out_pdb = tmp_path / "apo_output.pdb"
    spec = SystemSpec.from_tag("apo", variant="WT")
    assemble_system_model(spec=spec, receptor_pdb=inp_pdb, output_pdb=out_pdb)

    content = out_pdb.read_text()
    assert "TYR A 171" in content
    assert "PTR" not in content
    assert "O1P" not in content
    assert "ATP" not in content
    assert "MG" not in content
    assert " B " not in content


def test_position_nucleotide(tmp_path):
    from colabmda.modeller.substrate import position_nucleotide

    # Receptor with CA atoms in 3D but NO ATP
    rec_lines = [
        f"ATOM  {i:5d}  CA  ALA A{i:4d}    {float(i):8.3f}{float(i*2):8.3f}{float(i%5):8.3f}  1.00  0.00           C\n"
        for i in range(1, 30)
    ] + ["TER\n", "END\n"]
    rec_pdb = tmp_path / "rec_no_atp.pdb"
    rec_pdb.write_text("".join(rec_lines))

    # Template with matching CA atoms and an ATP molecule
    tmpl_lines = [
        f"ATOM  {i:5d}  CA  ALA A{i:4d}    {float(i):8.3f}{float(i*2):8.3f}{float(i%5):8.3f}  1.00  0.00           C\n"
        for i in range(1, 30)
    ] + [
        "HETATM 2974  PG  ATP A 300     -23.284  -1.633  18.344  1.00  6.41           P\n",
        "HETATM 2975  O1G ATP A 300     -23.539  -0.309  19.041  1.00  7.63           O\n",
        "TER\n",
        "END\n",
    ]
    tmpl_pdb = tmp_path / "tmpl_with_atp.pdb"
    tmpl_pdb.write_text("".join(tmpl_lines))

    out_pdb = tmp_path / "rec_with_placed_atp.pdb"
    success = position_nucleotide(
        receptor_pdb=rec_pdb,
        output_pdb=out_pdb,
        template_pdb=tmpl_pdb,
        ligand_name="ATP",
    )
    assert success is True
    content = out_pdb.read_text()
    assert "ATP" in content
    assert "PG" in content
