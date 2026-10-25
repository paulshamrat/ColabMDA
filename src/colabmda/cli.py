#!/usr/bin/env python3
import argparse
import os
import re
import shutil
from pathlib import Path

DEFAULT_DRIVE_ROOT = "/content/drive/MyDrive/projects"
ENV_ROOT = "COLABMDA_ROOT"
ENV_PROJECT = "COLABMDA_PROJECT_DIR"
ROOT_HELP = "Override project root (defaults to current working directory)"


def _normalize_equil_protocol(protocol: str) -> str:
    """Map user-facing equilibration labels to internal implementation names."""
    p_norm = str(protocol).lower().strip()
    if p_norm in {"default", "full", "staged"}:
        return "varmdyn"
    if p_norm in {"varmdyn", "quick"}:
        return p_norm
    raise SystemExit(
        "ERROR: --protocol must be 'staged' for the full staged protocol "
        "or 'quick' for a short smoke-test equilibration."
    )


def _resolve_root(use_drive: bool, root: str | None) -> str:
    # 1. Use explicit root if provided
    if root:
        return str(Path(root).resolve())

    # 2. Use environment variable if provided
    for env_var in [ENV_PROJECT, ENV_ROOT]:
        val = os.environ.get(env_var)
        if val:
            return str(Path(val).resolve())

    cwd = os.getcwd()
    # 3. If on Colab but in the temporary /content folder, check Drive
    if cwd == "/content" and os.path.exists("/content/drive/MyDrive"):
        if os.path.exists("/content/drive/MyDrive/projects"):
            return "/content/drive/MyDrive/projects"
        if os.path.exists("/content/drive/MyDrive/ColabMDA"):
            return "/content/drive/MyDrive/ColabMDA"
        if os.path.exists("/content/drive/MyDrive/research/projects"):
            return "/content/drive/MyDrive/research/projects"
        return "/content/drive/MyDrive/projects"

    # 4. Otherwise, use the Current Working Directory (standard for HPC/Laptop/Colab projects)
    return cwd


def _ensure_dir(path: str):
    Path(path).mkdir(parents=True, exist_ok=True)


def _prepare_run_inputs(root: str, pdbid: str, name: str):
    prep_dir = Path(root) / pdbid / "prep"
    run_dir = Path(root) / pdbid / "run"
    run_dir.mkdir(parents=True, exist_ok=True)

    prep_clean = prep_dir / f"{name}_cleaned.pdb"
    run_clean = run_dir / f"{name}_cleaned.pdb"
    if not run_clean.exists() and prep_clean.exists():
        shutil.copy2(prep_clean, run_clean)

    prep_raw = prep_dir / f"{name}.pdb"
    run_raw = run_dir / f"{name}.pdb"
    if not run_raw.exists() and prep_raw.exists():
        shutil.copy2(prep_raw, run_raw)

    if not run_clean.exists():
        raise SystemExit(
            f"ERROR: cleaned PDB missing for run.\n"
            f"  Expected at: {run_clean}\n"
            f"Run prep first, e.g.:\n"
            f"  colabmda openmm prep --pdb-id {pdbid}\n"
        )


def _guess_pdbid_from_workdir(workdir: str) -> str | None:
    candidates = list(Path(workdir).glob("*_cleaned.pdb"))
    if not candidates:
        return None
    if len(candidates) > 1:
        names = ", ".join([c.name for c in candidates])
        raise SystemExit(
            f"ERROR: multiple *_cleaned.pdb files found in {workdir}: {names}\n"
            f"Please specify --name."
        )
    return candidates[0].name.replace("_cleaned.pdb", "")


def _default_project_root() -> str:
    return _resolve_root(use_drive=True, root=None) or DEFAULT_DRIVE_ROOT


def _simulation_workdir(root: str, name: str) -> Path:
    """Resolve a named simulation under the current and legacy layouts."""
    root_path = Path(root).resolve()
    candidates = [
        root_path / "data" / "simulation" / name,
        root_path / "data" / "sim" / name,
        root_path / "simulation" / name,
        root_path / "sim" / name,
        root_path / "simulations" / name,
    ]
    return next((path for path in candidates if path.is_dir()), candidates[0])


def _analysis_outdir(root: str, name: str) -> Path:
    """Resolve a named analysis directory under the current and legacy layouts."""
    root_path = Path(root).resolve()
    candidates = [
        root_path / "data" / "analysis" / name,
        root_path / "analysis" / name,
    ]
    return next((path for path in candidates if path.is_dir()), candidates[0])


def _dispatch_prep(args) -> int:
    from colabmda.amber.leap import generate_tleap_script, run_tleap, write_prep_provenance
    from colabmda.core.ions import parse_ion_spec
    from colabmda.core.project import ProjectPaths, resolve_project_root
    from colabmda.core.taxonomy import SystemSpec
    from colabmda.modeller.substrate import assemble_system_model, fetch_rcsb_pdb

    root = resolve_project_root(args.root)
    paths = ProjectPaths(root, state=getattr(args, "state", None))
    paths.init_layout()

    if args.state:
        spec = SystemSpec.from_tag(args.state, variant=args.variant, protein=args.protein)
    else:
        spec = SystemSpec(protein=args.protein, variant=args.variant)

    if args.ligand:
        spec.ligand = args.ligand
    if args.ions:
        spec.ions = parse_ion_spec(args.ions)
    if args.phos:
        spec.phos_residues = [p.strip() for p in args.phos.split(",") if p.strip()]
    if args.peptide_seq:
        spec.peptide_seq = args.peptide_seq
    if args.peptide_name:
        spec.peptide_name = args.peptide_name

    # Receptor resolution
    if paths.state_dir is not None and not args.receptor:
        model = paths.structure_dir / "models" / spec.variant
        candidate = model / ("complex.pdb" if not spec.is_apo else f"{spec.variant}.pdb")
        if not candidate.is_file():
            raise FileNotFoundError(f"Verified state/variant receptor missing: {candidate}")
        args.receptor = str(candidate)
    if args.receptor:
        rec_path = Path(args.receptor).resolve()
    else:
        rec_path = None
        # 1. If variant is specified and not WT, prioritize variant-specific receptor
        if args.variant and args.variant.upper() != "WT":
            var_candidates = [
                f"{args.protein}_{args.variant}.pdb",
                f"{args.protein}_{args.variant.lower()}.pdb",
                f"{args.variant}.pdb",
                f"{args.variant.lower()}.pdb",
                f"target.B99990001_{args.variant}_with_cryst.pdb",
                f"target.B99990001_{args.variant}.pdb",
            ]
            for cand in var_candidates:
                rec_path = paths.find_template(cand)
                if rec_path:
                    break
            if not rec_path and (paths.variants_dir / spec.variant / "model.pdb").exists():
                rec_path = paths.variants_dir / spec.variant / "model.pdb"
            if not rec_path and (paths.variants_dir / f"{spec.variant}.pdb").exists():
                rec_path = paths.variants_dir / f"{spec.variant}.pdb"

        # 2. Otherwise look for WT / base protein receptor
        if not rec_path:
            wt_candidates = [
                f"{args.protein}_wt.pdb",
                f"{args.protein}_WT.pdb",
                f"{args.protein}.pdb",
                "target.B99990001_with_cryst.pdb",
                "target.B99990001.pdb",
                "4bgq.pdb",
            ]
            for cand in wt_candidates:
                rec_path = paths.find_template(cand)
                if rec_path:
                    break

        # 3. Direct fallback to protein_dir
        if not rec_path and paths.protein_dir.exists():
            if (paths.protein_dir / f"{args.protein}.pdb").exists():
                rec_path = paths.protein_dir / f"{args.protein}.pdb"
            else:
                pdb_files = list(paths.protein_dir.glob("*.pdb"))
                if len(pdb_files) == 1:
                    rec_path = pdb_files[0]

    # Template resolution
    tmpl_path = None
    if args.template:
        if len(args.template) == 4 and not Path(args.template).exists():
            tmpl_path = fetch_rcsb_pdb(args.template, paths.template_dir)
        else:
            tmpl_path = Path(args.template).resolve()
    else:
        tmpl_path = paths.find_template("3X2V.pdb")
        if not tmpl_path and (
            spec.peptide_seq or spec.peptide_name or (spec.ions and len(spec.ions) > 1)
        ):
            tmpl_path = fetch_rcsb_pdb("3X2V", paths.template_dir)

    if paths.state_dir is not None:
        if not rec_path or not rec_path.exists():
            raise FileNotFoundError("Provide a reviewed --receptor")
        from colabmda.core.workflow import preparation_request

        if preparation_request(paths, spec, rec_path, tmpl_path):
            return 0
    out_model = paths.prep_dir(spec) / "model.pdb"
    out_model.parent.mkdir(parents=True, exist_ok=True)

    if rec_path and rec_path.exists():
        assemble_system_model(spec, rec_path, out_model, tmpl_path)
        print(f"[OK] Assembled model coordinates: {out_model}")
    else:
        print(f"[WARN] No receptor PDB found for {args.protein}. Please provide --receptor.")
        return 1

    if not args.no_amber:
        try:
            solvated_out = paths.prepared_solvated(spec)
            script_text = generate_tleap_script(
                spec=spec,
                input_pdb=out_model,
                output_prmtop=paths.prepared_prmtop(spec),
                output_inpcrd=paths.prepared_inpcrd(spec),
                output_pdb=solvated_out,
            )
            retcode, log_out = run_tleap(
                script_text,
                workdir=paths.prep_dir(spec),
                keep_script=False,
            )
            if retcode == 0 and paths.prepared_prmtop(spec).exists():
                write_prep_provenance(spec=spec, dest_json=paths.prep_json(spec))
                print(f"[OK] Built Amber topology: {paths.prepared_prmtop(spec)}")
                print(f"[OK] Built Amber coordinates: {paths.prepared_inpcrd(spec)}")
                if solvated_out.exists():
                    print(f"[OK] Built solvated model: {solvated_out}")
                print(f"[OK] Generated provenance metadata: {paths.prep_json(spec)}")
                if paths.state_dir is not None:
                    from colabmda.core.proteins import check_topology

                    check_topology(paths.prep_dir(spec), spec=spec)
                    (paths.prep_dir(spec) / "PREP_DONE").write_text(
                        "Verified canonical coordinates, topology composition, charges and outputs\n"
                    )
            else:
                (paths.prep_dir(spec) / "build_error.log").write_text(log_out)
                raise RuntimeError(
                    f"Amber tleap failed with returncode {retcode}; see build_error.log"
                )
        except Exception as e:
            raise RuntimeError(f"Amber preparation failed: {e}") from e

    return 0


def _dispatch_stage(args) -> int:
    if getattr(args, "example", None):
        from colabmda.core.examples import install_example

        if not args.root:
            raise ValueError("--example requires an explicit --root")
        return install_example(args.example, args.root)

    from colabmda.core.project import ProjectPaths, resolve_project_root
    from colabmda.core.taxonomy import SystemSpec

    root = resolve_project_root(args.root)
    paths = ProjectPaths(root, state=getattr(args, "state", None))
    paths.init_layout()

    spec = SystemSpec.from_tag(args.state, variant=args.variant, protein=args.protein)
    paths.log_dir(spec).mkdir(parents=True, exist_ok=True)
    paths.prep_dir(spec).mkdir(parents=True, exist_ok=True)
    paths.eq_dir(spec).mkdir(parents=True, exist_ok=True)
    paths.analysis_dir(spec).mkdir(parents=True, exist_ok=True)

    replicas = [r.strip() for r in args.replicas.split(",") if r.strip()]
    for rep in replicas:
        paths.replica_dir(spec, rep).mkdir(parents=True, exist_ok=True)

    if not getattr(args, "scripts", False):
        print(f"[OK] Staged {paths.variant_dir(spec)}; use direct cmda commands")
        return 0

    run_sh = paths.scripts_dir / "run.sh"
    run_reps = " ".join(f"r{r}" if r.isdigit() else r for r in replicas)
    run_sh_content = f"""#!/usr/bin/env bash
set -euo pipefail

STATE="{spec.state_tag}"
VARIANT="{spec.variant}"

echo "=== Running ColabMDA Simulation Pipeline ==="
echo "State:   ${{STATE}}"
echo "Variant: ${{VARIANT}}"

# 1. Equilibration (em -> nvt -> npt -> qc)
cmda sim --state "${{STATE}}" --variant "${{VARIANT}}" --stage em,nvt,npt,qc

# 2. Production Replicas
for rep in {run_reps}; do
    echo "[RUN] Starting production replica ${{rep}}..."
    cmda sim --state "${{STATE}}" --variant "${{VARIANT}}" --stage md --replica "${{rep}}"
done

# 3. Trajectory Autoimaging & Dual-Stream Stripping
cmda merge --state "${{STATE}}" --variant "${{VARIANT}}"

# 4. Analysis
cmda analyze --state "${{STATE}}" --variant "${{VARIANT}}"
"""
    run_sh.write_text(run_sh_content)
    run_sh.chmod(0o755)

    # HPC SLURM scripts
    eq_slurm = paths.scripts_dir / "eq.slurm"
    eq_slurm_compat = paths.scripts_dir / "submit_eq.slurm"
    eq_slurm_content = f"""#!/usr/bin/env bash
#SBATCH --job-name=eq_{spec.state_tag}_{spec.variant}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gpus=1
#SBATCH --mem=32gb
#SBATCH --time=24:00:00
#SBATCH --output=logs/eq_%j.out
#SBATCH --error=logs/eq_%j.err

set -euo pipefail

# Ensure execution runs from project root
SCRIPT_DIR="$(cd "$(dirname "${{BASH_SOURCE[0]}}")" && pwd)"
PROJECT_ROOT="$(cd "${{SCRIPT_DIR}}/.." && pwd)"
cd "${{PROJECT_ROOT}}"

# Activate environment
if [ -f "$HOME/miniforge3/etc/profile.d/conda.sh" ]; then
    source "$HOME/miniforge3/etc/profile.d/conda.sh"
    conda activate colabmda 2>/dev/null || conda activate openmm_env 2>/dev/null || true
fi

STATE="{spec.state_tag}"
VARIANT="{spec.variant}"

echo "[SLURM] Running equilibration for ${{STATE}}/${{VARIANT}} on $(hostname)"
mkdir -p logs

cmda sim --state "${{STATE}}" --variant "${{VARIANT}}" --stage eq
"""
    eq_slurm.write_text(eq_slurm_content)
    eq_slurm.chmod(0o755)
    eq_slurm_compat.write_text(eq_slurm_content)
    eq_slurm_compat.chmod(0o755)

    prod_slurm = paths.scripts_dir / "prod.slurm"
    prod_slurm_compat = paths.scripts_dir / "submit_prod.slurm"
    prod_slurm_content = f"""#!/usr/bin/env bash
#SBATCH --job-name=cmda_{spec.state_tag}_{spec.variant}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gpus=1
#SBATCH --mem=32gb
#SBATCH --time=48:00:00
#SBATCH --array=1-{len(replicas)}
#SBATCH --output=logs/prod_%A_%a.out
#SBATCH --error=logs/prod_%A_%a.err

set -euo pipefail

# Ensure execution runs from project root
SCRIPT_DIR="$(cd "$(dirname "${{BASH_SOURCE[0]}}")" && pwd)"
PROJECT_ROOT="$(cd "${{SCRIPT_DIR}}/.." && pwd)"
cd "${{PROJECT_ROOT}}"

# Activate environment
if [ -f "$HOME/miniforge3/etc/profile.d/conda.sh" ]; then
    source "$HOME/miniforge3/etc/profile.d/conda.sh"
    conda activate colabmda 2>/dev/null || conda activate openmm_env 2>/dev/null || true
fi

STATE="{spec.state_tag}"
VARIANT="{spec.variant}"
REPLICA="r${{SLURM_ARRAY_TASK_ID}}"

echo "[SLURM] Running replica ${{REPLICA}} on $(hostname)"
mkdir -p logs

cmda sim --state "${{STATE}}" --variant "${{VARIANT}}" --stage md --replica "${{REPLICA}}"
"""
    prod_slurm.write_text(prod_slurm_content)
    prod_slurm.chmod(0o755)
    prod_slurm_compat.write_text(prod_slurm_content)
    prod_slurm_compat.chmod(0o755)

    submit_slurm = paths.scripts_dir / "submit.slurm"
    submit_slurm.write_text(prod_slurm_content)
    submit_slurm.chmod(0o755)

    chain_sh = paths.scripts_dir / "chain.sh"
    chain_sh_compat = paths.scripts_dir / "submit_chain.sh"
    chain_content = """#!/usr/bin/env bash
# Submit equilibration first, then chain production replicas as a dependency
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${PROJECT_ROOT}"
mkdir -p logs

echo "Submitting equilibration job..."
EQ_ID=$(sbatch scripts/eq.slurm | awk '{print $NF}')
echo "Equilibration Job ID: ${EQ_ID}"

echo "Submitting production replica array (dependent on ${EQ_ID})..."
PROD_ID=$(sbatch --dependency=afterok:${EQ_ID} scripts/prod.slurm | awk '{print $NF}')
echo "Production Array Job ID: ${PROD_ID}"

echo "Watch queue with: squeue -u $USER"
"""
    chain_sh.write_text(chain_content)
    chain_sh.chmod(0o755)
    chain_sh_compat.write_text(chain_content)
    chain_sh_compat.chmod(0o755)

    print(
        f"[OK] Staged project directory tree for {spec.state_tag}/{spec.variant} at: {paths.variant_dir(spec)}"
    )
    print(f"[OK] Generated decoupled runner: {run_sh}")
    print(f"[OK] Generated HPC equilibration SLURM: {eq_slurm}")
    print(f"[OK] Generated HPC production SLURM array: {prod_slurm}")
    print(f"[OK] Generated HPC chained submission helper: {chain_sh}")
    return 0


def _dispatch_sim(args) -> int:
    from colabmda.core.project import ProjectPaths, resolve_project_root
    from colabmda.core.taxonomy import SystemSpec
    from colabmda.openmm.commands import (
        openmm_em,
        openmm_md,
        openmm_npt,
        openmm_nvt,
        openmm_qc,
    )

    root = resolve_project_root(args.root)
    paths = ProjectPaths(root, state=getattr(args, "state", None))
    spec = SystemSpec.from_tag(args.state, variant=args.variant)

    raw_stages = [s.strip().lower() for s in args.stage.split(",") if s.strip()]
    stages = []
    for s in raw_stages:
        if s in {"eq", "equil", "equilibration"}:
            stages.extend(["em", "nvt", "npt", "qc"])
        else:
            stages.append(s)
    eq_dir = str(paths.eq_dir(spec))
    paths.eq_dir(spec).mkdir(parents=True, exist_ok=True)

    model_src = paths.model_pdb(spec)
    name = f"{spec.state_tag}_{spec.variant}"
    eq_clean = paths.eq_dir(spec) / f"{name}_cleaned.pdb"
    if model_src.exists() and not eq_clean.exists():
        shutil.copy2(model_src, eq_clean)

    prmtop = str(paths.prepared_prmtop(spec)) if paths.prepared_prmtop(spec).exists() else None
    inpcrd = str(paths.prepared_inpcrd(spec)) if paths.prepared_inpcrd(spec).exists() else None

    protocol = _normalize_equil_protocol(getattr(args, "protocol", "staged"))
    software_test = getattr(args, "software_test", False)
    if os.environ.get("COLABMDA_TEST_1_24", "").lower() in {"1", "true"} and not software_test:
        raise SystemExit("Reduced 1-24 test steps require --software-test")
    if protocol == "quick" and not software_test:
        raise SystemExit(
            "Quick equilibration requires --software-test; research uses --protocol staged"
        )

    if paths.state_dir is not None and not software_test:
        from colabmda.core.workflow import verify_prepared

        verify_prepared(paths.prep_dir(spec))

    for stage in stages:
        print(f"\n[SIM] Executing stage: {stage.upper()} for {spec.state_tag}/{spec.variant}")
        if stage == "em":
            max_iter = 200 if protocol == "quick" else 0
            openmm_em(
                eq_dir,
                name,
                amber_prmtop=prmtop,
                amber_inpcrd=inpcrd,
                max_iterations=max_iter,
            )
        elif stage == "nvt":
            eq_time = args.equil_time
            if protocol == "quick" and getattr(args, "steps", None) is not None:
                eq_time = min(args.equil_time, max(0.002, args.steps * 0.002))
            openmm_nvt(
                eq_dir,
                name,
                equil_time=eq_time,
                protocol=protocol,
            )
        elif stage == "npt":
            eq_time = args.equil_time
            if protocol == "quick" and getattr(args, "steps", None) is not None:
                eq_time = min(args.equil_time, max(0.002, args.steps * 0.002))
            openmm_npt(
                eq_dir,
                name,
                equil_time=eq_time,
                protocol=protocol,
            )
        elif stage in ["qc", "check-equil"]:
            openmm_qc(eq_dir, strict=not software_test, protocol=protocol)
        elif stage == "md":
            if not software_test:
                import json

                recorded = Path(eq_dir) / "npt_protocol.json"
                if (
                    not recorded.exists()
                    or json.loads(recorded.read_text()).get("protocol") != "varmdyn"
                ):
                    raise SystemExit("Research production requires completed staged equilibration")
                openmm_qc(eq_dir, strict=True, protocol="varmdyn")
            replicas = [r.strip() for r in args.replica.split(",") if r.strip()]
            for rep in replicas:
                rep_dir = str(paths.replica_dir(spec, rep))
                paths.replica_dir(spec, rep).mkdir(parents=True, exist_ok=True)
                for artifact in ["npt.state.xml", "system.xml", "solvated.pdb", "restart.chk"]:
                    eq_art = paths.eq_dir(spec) / artifact
                    rep_art = paths.replica_dir(spec, rep) / artifact
                    if eq_art.exists() and not rep_art.exists():
                        shutil.copy2(eq_art, rep_art)
                rep_clean = paths.replica_dir(spec, rep) / f"{name}_cleaned.pdb"
                if eq_clean.exists() and not rep_clean.exists():
                    shutil.copy2(eq_clean, rep_clean)

                total_ns = args.total_ns
                traj_interval = getattr(args, "traj_interval", 10.0)
                if args.steps is not None:
                    total_ns = (args.steps * 0.002) / 1000.0
                    total_ps = total_ns * 1000.0
                    traj_interval = min(traj_interval, max(0.002, total_ps / 5.0))

                openmm_md(
                    rep_dir,
                    name,
                    total_ns=total_ns,
                    traj_interval=traj_interval,
                    checkpoint_ps=args.checkpoint_ps,
                )
    return 0


def _dispatch_merge(args) -> int:
    from colabmda.amber.cpptraj import generate_cpptraj_strip_script, run_cpptraj
    from colabmda.core.project import ProjectPaths, resolve_project_root
    from colabmda.core.taxonomy import SystemSpec
    from colabmda.openmm.merge import merge_logs

    paths = ProjectPaths(resolve_project_root(args.root), state=args.state)
    spec = SystemSpec.from_tag(args.state, variant=args.variant)
    for rep in [r.strip() for r in args.replicas.split(",") if r.strip()]:
        directory = paths.replica_dir(spec, rep)
        chunks = []
        for path in directory.glob("prod_*to*ps.dcd"):
            match = re.fullmatch(r"prod_([0-9.]+)to([0-9.]+)ps.dcd", path.name)
            if match:
                chunks.append((float(match[1]), float(match[2]), path))
        chunks.sort()
        if not chunks or chunks[0][0] != 0 or any(a[1] != b[0] for a, b in zip(chunks, chunks[1:])):
            raise RuntimeError(
                f"Missing, overlapping or non-contiguous trajectory chunks: {directory}"
            )
        for start, end, path in chunks:
            if end <= start or path.stat().st_size == 0:
                raise RuntimeError(f"Invalid trajectory chunk: {path}")
        for mode, filename, top in [
            ("keep_ligand", "prod_keep_ligand.dcd", "stripped.prmtop"),
            ("ca_backbone", "prod_ca.dcd", "ca.prmtop"),
        ]:
            script = generate_cpptraj_strip_script(
                paths.prepared_prmtop(spec),
                [c[2] for c in chunks],
                directory / filename,
                output_top=directory / top,
                mode=mode,
            )
            rc, output = run_cpptraj(script, directory)
            if rc or not (directory / filename).is_file() or not (directory / top).is_file():
                raise RuntimeError(f"cpptraj {mode} failed: {output[-3000:]}")
        merge_logs(str(directory), str(directory / "prod_full.log"))
        (directory / "MERGE_DONE").write_text(
            "Verified trajectory inputs and both stripped outputs\n"
        )
    return 0


def _dispatch_analyze(args) -> int:
    from colabmda.core.project import ProjectPaths, resolve_project_root
    from colabmda.core.taxonomy import SystemSpec
    from colabmda.openmm.commands import openmm_analysis

    root = resolve_project_root(args.root)
    paths = ProjectPaths(root, state=getattr(args, "state", None))
    spec = SystemSpec.from_tag(args.state, variant=args.variant)
    name = f"{spec.state_tag}_{spec.variant}"
    analysis_dir = paths.analysis_dir(spec)
    analysis_dir.mkdir(parents=True, exist_ok=True)

    print(f"[OK] Running analysis for {name} -> {analysis_dir}")
    rep_dir = paths.replica_dir(spec, "r1")
    traj = None
    top = None
    for cand_traj in ["prod_keep_ligand.dcd", f"{name}_merged.dcd", "prod_ca.dcd"]:
        if (rep_dir / cand_traj).exists():
            traj = str(rep_dir / cand_traj)
            break
    for cand_top in ["stripped.system.prmtop", "stripped.prmtop", "system.prmtop"]:
        if (rep_dir / cand_top).exists():
            top = str(rep_dir / cand_top)
            break
        elif (paths.prep_dir(spec) / cand_top).exists():
            top = str(paths.prep_dir(spec) / cand_top)
            break

    if traj and top:
        openmm_analysis(
            pdbid_dir=str(rep_dir),
            topology=top,
            trajectory=traj,
            interval_ps=10.0,
            outdir=str(analysis_dir),
            reference=str(paths.model_pdb(spec)),
        )
    else:
        raise FileNotFoundError(f"No suitable trajectory/topology found in {rep_dir}")
    return 0


def _dispatch_compare(args) -> int:
    from colabmda.core.project import ProjectPaths, resolve_project_root

    root = resolve_project_root(args.root)
    paths = ProjectPaths(root, state=getattr(args, "state", None))
    variants = [v.strip() for v in args.variants.split(",") if v.strip()]
    print(f"[OK] Comparing variants {variants} for state {args.state} under {paths.root}")
    return 0


def main():
    p = argparse.ArgumentParser(prog="colabmda")
    sub = p.add_subparsers(dest="tool", required=True)

    fetch = sub.add_parser(
        "fetch", help="Retrieve an archive and verify author/UniProt residue mapping"
    )
    fetch.add_argument("--pdb-id", required=True)
    fetch.add_argument("--uniprot-id", required=True)
    fetch.add_argument("--chain", default="A")
    fetch.add_argument("--range", nargs=2, type=int, required=True)
    fetch.add_argument(
        "--mutations", default="", help="Comma-separated WT mutation identities to verify"
    )
    fetch.add_argument("--state", required=True)
    fetch.add_argument("--root", default=None, help=ROOT_HELP)

    model = sub.add_parser("model", help="Build WT or a variant from checked public project inputs")
    model.add_argument("--state", required=True)
    model.add_argument("--variant", default="WT")
    model.add_argument("--seed", type=int, default=13031)
    model.add_argument("--root", default=None, help=ROOT_HELP)

    check = sub.add_parser(
        "check", help="Verify a canonical WT/mutant model against frozen structure inputs"
    )
    check.add_argument("--state", required=True)
    check.add_argument("--variant", default="WT")
    check.add_argument("--model", default=None)
    check.add_argument("--root", default=None, help=ROOT_HELP)

    assemble = sub.add_parser(
        "assemble", help="Assemble ATP/Mg from a public donor or an explicit reviewed pose"
    )
    for name in ["receptor", "reference", "ligand", "ions", "mapping", "out"]:
        assemble.add_argument("--" + name)
    assemble.add_argument("--state", help="Use verified project models for public pose transfer")
    assemble.add_argument("--template", help="Public ATP/Mg donor archive ID (e.g. 8FP5)")
    assemble.add_argument(
        "--relax-sidechains",
        action="store_true",
        help="Locally minimize pocket sidechains with backbone/ATP/Mg fixed",
    )
    assemble.add_argument("--root", default=None, help=ROOT_HELP)
    assemble.add_argument("--chain", default="A")
    assemble.add_argument("--variant", default="WT")
    assemble.add_argument("--core", nargs=2, type=int, default=[30, 220])
    assemble.add_argument("--max-rmsd", type=float, default=1.0)

    # ---------------- OpenMM ----------------
    p_openmm = sub.add_parser("openmm", help="OpenMM protein-water workflow")
    sub_openmm = p_openmm.add_subparsers(dest="cmd", required=True)

    # prep
    p_prep = sub_openmm.add_parser("prep", help="Prepare/clean PDB")
    g = p_prep.add_mutually_exclusive_group(required=True)
    g.add_argument("--pdb-id", help="4-letter PDB id (downloads from RCSB)")
    g.add_argument("--pdb-file", help="Local PDB file path (no download)")
    p_prep.add_argument(
        "--outdir", default=None, help="Output directory for --pdb-file (default: ./<name>)"
    )
    p_prep.add_argument(
        "--name", default=None, help="Prefix name (default: from --pdb-id or file stem)"
    )
    p_prep.add_argument("--ph", type=float, default=7.4, help="Hydrogen pH (default: 7.4)")
    p_prep.add_argument(
        "--sync-dir", default=None, help="Optional: copy cleaned prep outputs to this directory"
    )
    p_prep.add_argument(
        "--drive", action="store_true", help="(compat) Use Drive root (default behavior)"
    )
    p_prep.add_argument(
        "--root",
        "--project",
        dest="root",
        default=None,
        help=ROOT_HELP,
    )

    # run (colab-safe runner)
    p_run = sub_openmm.add_parser("run", help="Run/resume chunked MD (colab-safe)")
    g = p_run.add_mutually_exclusive_group(required=False)
    g.add_argument("--pdb-id", help="Use ./<pdb-id> as workdir (pdb-id download workflow)")
    g.add_argument("--workdir", help="Folder containing <name>_cleaned.pdb")
    p_run.add_argument(
        "--name", default=None, help="Prefix name (default: --pdb-id or inferred from workdir)"
    )
    p_run.add_argument("--total-ns", type=float, default=100.0)
    p_run.add_argument(
        "--traj-interval",
        type=float,
        default=10.0,
        help="ps between saved frames (default: 10.0 ps)",
    )
    p_run.add_argument(
        "--equil-time",
        type=float,
        default=100.0,
        help="NVT and NPT duration in ps for --equil-protocol quick",
    )
    p_run.add_argument("--checkpoint-ps", type=float, default=1000.0, help="ps per chunk")
    p_run.add_argument("--sync-dir", default=None, help="Optional: sync outputs to this directory")
    p_run.add_argument(
        "--replica", default=None, help="Optional: replica subfolder (e.g. r1, r2, r3)"
    )
    p_run.add_argument("--ligand", default=None, help="Optional: Path to ligand SDF/MOL2 file")
    p_run.add_argument(
        "--small-molecule-ff",
        default="gaff-2.11",
        help="Small-molecule force field for ligand setup (default: gaff-2.11; e.g. gaff-2.2.20 if installed)",
    )
    p_run.add_argument(
        "--keep-mg", action="store_true", help="Optional: Keep Mg2+ ion from raw structure"
    )
    p_run.add_argument(
        "--padding-nm",
        default="auto",
        help="Solvent padding in nm, or auto; 1.4 is conservative and 1.0 is compact",
    )
    p_run.add_argument(
        "--protein-ff",
        choices=["ff19SB", "ff14SB"],
        default="ff19SB",
        help="Protein force field for native OpenMM setup (default: ff19SB)",
    )
    p_run.add_argument(
        "--water-model",
        choices=["opc", "tip3p", "tip3pfb", "tip4pew"],
        default="opc",
        help="Water/ion model for native OpenMM setup (default: opc)",
    )
    p_run.add_argument(
        "--equil-protocol",
        default="default",
        metavar="{default,quick}",
        help="Full staged equilibration or short smoke-test equilibration",
    )
    p_run.add_argument("--amber-prmtop", default=None, help="Pre-parameterized AMBER topology")
    p_run.add_argument("--amber-inpcrd", default=None, help="Matching AMBER coordinates")
    p_run.add_argument(
        "--seed", type=int, default=None, help="Optional: random seed for velocity assignment"
    )
    p_run.add_argument(
        "--equil-only",
        action="store_true",
        help="Run Energy Minimization, NVT, and NPT equilibration, then exit before production MD",
    )
    p_run.add_argument(
        "--drive", action="store_true", help="(compat) Use Drive root (default behavior)"
    )
    p_run.add_argument(
        "--root",
        "--project",
        dest="root",
        default=None,
        help=ROOT_HELP,
    )

    # merge
    p_merge = sub_openmm.add_parser("merge", help="Merge chunk DCDs/logs")
    g = p_merge.add_mutually_exclusive_group(required=False)
    g.add_argument("--pdb-id", help="Use ./<pdb-id> as simulation directory")
    g.add_argument("--pdb-dir", help="Simulation directory (e.g. 4ldj_wt)")
    g.add_argument("--name", help="Simulation system name (e.g. 4ldj_wt)")
    p_merge.add_argument("--replica", default="r1", help="Replica subfolder (default: r1)")
    p_merge.add_argument(
        "--topology", default=None, help="Topology PDB (default: <dir>/solvated.pdb)"
    )
    p_merge.add_argument("--out-traj", default="prod_full.dcd")
    p_merge.add_argument("--out-log", default="prod_full.log")
    p_merge.add_argument(
        "--stride", type=int, default=1, help="Keep every Nth frame while merging (default: 1)"
    )
    p_merge.add_argument("--center", action="store_true", help="Center protein in the box")
    p_merge.add_argument(
        "--wrap", action="store_true", help="Wrap solvent molecules (image_molecules)"
    )
    p_merge.add_argument(
        "--mda", action="store_true", help="Use MDAnalysis for merging and PBC correction"
    )
    p_merge.add_argument(
        "--selection", default=None, help="MDAnalysis selection string (implies --mda)"
    )
    p_merge.add_argument(
        "--ca-only", action="store_true", help="Extract Cα atoms only (implies --mda)"
    )
    p_merge.add_argument(
        "--protein-only", action="store_true", help="Extract protein-only atoms (implies --mda)"
    )
    p_merge.add_argument(
        "--drive", action="store_true", help="(compat) Use Drive root (default behavior)"
    )
    p_merge.add_argument(
        "--root",
        "--project",
        dest="root",
        default=None,
        help=ROOT_HELP,
    )

    # analysis
    p_ana = sub_openmm.add_parser("analysis", help="RMSD/Rg/RMSF analysis")
    g = p_ana.add_mutually_exclusive_group(required=False)
    g.add_argument("--pdb-id", help="Use ./<pdb-id> as simulation directory")
    g.add_argument("--pdb-dir", help="Simulation directory (e.g. 4ldj_wt)")
    g.add_argument("--name", help="Simulation system name (e.g. 4ldj_wt)")
    p_ana.add_argument(
        "--replica",
        default=None,
        help="Replica subfolder (default: all replicas if present, or r1)",
    )
    p_ana.add_argument("--topology", default=None)
    p_ana.add_argument("--trajectory", default=None)
    p_ana.add_argument(
        "--reference", default=None, help="Native-numbered reference structure for Amber topology"
    )
    p_ana.add_argument(
        "--interval", type=float, default=None, help="ps per frame (if not auto-detected)"
    )
    p_ana.add_argument("--outdir", default=None)
    p_ana.add_argument(
        "--drive", action="store_true", help="(compat) Use Drive root (default behavior)"
    )
    p_ana.add_argument(
        "--root",
        "--project",
        dest="root",
        default=None,
        help=ROOT_HELP,
    )

    # status
    p_stat = sub_openmm.add_parser("status", help="Sanity-check frames/time/resume readiness")
    g = p_stat.add_mutually_exclusive_group(required=False)
    g.add_argument("--pdb-id", help="Use ./<pdb-id> as simulation directory")
    g.add_argument("--pdb-dir", help="Simulation directory (e.g. 4ldj_wt)")
    g.add_argument("--name", help="Simulation system name (e.g. 4ldj_wt)")
    p_stat.add_argument("--replica", default="r1", help="Replica subfolder (default: r1)")
    p_stat.add_argument(
        "--drive", action="store_true", help="(compat) Use Drive root (default behavior)"
    )
    p_stat.add_argument(
        "--root",
        "--project",
        dest="root",
        default=None,
        help=ROOT_HELP,
    )

    # compare
    p_comp = sub_openmm.add_parser("compare", help="Compare multiple systems (aggregate replicas)")
    p_comp.add_argument(
        "--series",
        action="append",
        required=True,
        help="LABEL=DIR1,DIR2 (e.g. WT=data/analysis/wt/r1,data/analysis/wt/r2)",
    )
    p_comp.add_argument("--outdir", required=True, help="Output directory for plots")

    # view
    p_view = sub_openmm.add_parser("view", help="View trajectory in PyMOL")
    p_view.add_argument(
        "--pdb-dir",
        default=None,
        help="Directory containing the merged files (default: current directory)",
    )
    p_view.add_argument(
        "-t", "--topology", default=None, help="Topology file (default: prod_full.pdb)"
    )
    p_view.add_argument(
        "-x", "--trajectory", default=None, help="Trajectory file (default: prod_full.dcd)"
    )
    p_view.add_argument(
        "-r", "--resi", type=int, default=12, help="Highlight residue index (default: 12)"
    )

    # snapshots
    p_snap = sub_openmm.add_parser(
        "snapshots", help="Generate comparative snapshot grid (e.g. 3x8 WT vs Mutants)"
    )
    p_snap.add_argument(
        "-c", "--config", default=None, help="Path to custom JSON configuration file"
    )
    p_snap.add_argument(
        "-o",
        "--output",
        default="data/analysis/snapshots/master_3x8_snapshots.png",
        help="Path to output image grid file",
    )
    p_snap.add_argument(
        "--temp-dir", default="scratch/master_grid_3x8", help="Temporary folder for panel rendering"
    )
    p_snap.add_argument("--font-path", default=None, help="Custom TrueType font path")

    # stage
    p_stage = sub_openmm.add_parser(
        "stage", help="Stage a WT/mutant structure into simulations/<name>"
    )
    p_stage.add_argument(
        "--pdb-file",
        required=True,
        help="Input structure PDB file (typically from data/structure/)",
    )
    p_stage.add_argument("--name", required=True, help="Simulation name (e.g. 4ldj_wt, 4ldj_G12C)")
    p_stage.add_argument(
        "--replica",
        default=None,
        help="Optional: create nested replica subfolder (e.g. r1, r2, r3)",
    )
    p_stage.add_argument("--ph", type=float, default=7.4, help="Hydrogen pH (default: 7.4)")
    p_stage.add_argument("--root", "--project", dest="root", default=None, help=ROOT_HELP)

    # em/nvt/npt/check-equil/md (individual modular steps)
    p_em = sub_openmm.add_parser("em", help="Modular: Minimization")
    p_em.add_argument("--name", required=True)
    p_em.add_argument("--workdir", default=None)
    p_em.add_argument("--root", "--project", dest="root", default=None)
    p_em.add_argument("--ligand", default=None, help="Optional: Path to ligand SDF/MOL2 file")
    p_em.add_argument(
        "--small-molecule-ff",
        default="gaff-2.11",
        help="Small-molecule force field for ligand setup",
    )
    p_em.add_argument(
        "--keep-mg", action="store_true", help="Optional: Keep Mg2+ ion from raw structure"
    )
    p_em.add_argument("--padding-nm", default="auto")
    p_em.add_argument("--protein-ff", choices=["ff19SB", "ff14SB"], default="ff19SB")
    p_em.add_argument(
        "--water-model", choices=["opc", "tip3p", "tip3pfb", "tip4pew"], default="opc"
    )
    p_em.add_argument("--amber-prmtop", default=None)
    p_em.add_argument("--amber-inpcrd", default=None)

    p_nvt = sub_openmm.add_parser("nvt", help="Modular: NVT Equilibration")
    p_nvt.add_argument("--name", required=True)
    p_nvt.add_argument("--equil-time", type=float, default=100.0)
    p_nvt.add_argument("--seed", type=int, default=None)
    p_nvt.add_argument("--equil-protocol", default="default", metavar="{default,quick}")
    p_nvt.add_argument("--workdir", default=None)
    p_nvt.add_argument("--root", "--project", dest="root", default=None)

    p_npt = sub_openmm.add_parser("npt", help="Modular: NPT Equilibration")
    p_npt.add_argument("--name", required=True)
    p_npt.add_argument("--equil-time", type=float, default=100.0)
    p_npt.add_argument("--seed", type=int, default=None)
    p_npt.add_argument("--equil-protocol", default="default", metavar="{default,quick}")
    p_npt.add_argument("--workdir", default=None)
    p_npt.add_argument("--root", "--project", dest="root", default=None)

    p_qc = sub_openmm.add_parser(
        "qc", aliases=["check-equil"], help="Modular: Stability Check & QC Plots"
    )
    p_qc.add_argument("--name", required=True)
    p_qc.add_argument("--workdir", default=None)
    p_qc.add_argument("--root", "--project", dest="root", default=None)
    p_qc.add_argument(
        "--warn-only", action="store_true", help="Report failed QC without a nonzero exit"
    )

    p_md = sub_openmm.add_parser("md", help="Modular: Production MD")
    p_md.add_argument("--name", required=True)
    p_md.add_argument("--total-ns", type=float, default=100.0)
    p_md.add_argument("--traj-interval", type=float, default=10.0)
    p_md.add_argument("--checkpoint-ps", type=float, default=1000.0)
    p_md.add_argument("--sync-dir", default=None)
    p_md.add_argument("--workdir", default=None)
    p_md.add_argument("--root", "--project", dest="root", default=None)
    p_md.add_argument("--replica", default="r1", help="Production replica subfolder (default: r1)")
    p_md.add_argument("--seed", type=int, default=None, help="Random seed for velocity assignment")

    # ---------------- Modeller ----------------
    p_mod = sub.add_parser("modeller", help="Modeller workflows")
    sub_mod = p_mod.add_subparsers(dest="cmd", required=True)

    # build
    p_build = sub_mod.add_parser("build", help="Build homology model")
    p_build.add_argument("--pdb-id", required=True, help="Template PDB ID (e.g. 1d6o)")
    p_build.add_argument("--uniprot-id", required=True, help="UniProt ID (e.g. O76039)")
    p_build.add_argument("--chain", default="A", help="Chain ID (default: A)")
    p_build.add_argument("--range", nargs=2, type=int, metavar=("START", "END"))
    p_build.add_argument("--truncate", action="store_true")
    p_build.add_argument("--uniprot-numbering", action="store_true")
    p_build.add_argument("--template", help="Reviewed local PDB template")
    p_build.add_argument("--sequence", help="Reviewed local UniProt FASTA")
    p_build.add_argument("--mut", default=None)
    p_build.add_argument("--list", default=None, help="File with one mutation per line")
    p_build.add_argument("--outdir", default=None)
    p_build.add_argument("--outdir-mut", default=None)
    p_build.add_argument("--seed", type=int, default=None)
    p_build.add_argument("--logfile", default=None)
    p_build.add_argument("--verbose", action="store_true")

    # mutate
    p_mut = sub_mod.add_parser("mutate", help="Mutate an existing PDB")
    p_mut.add_argument("--pdb-in", required=True, help="Input PDB to mutate")
    p_mut.add_argument("--chain", default="A", help="Chain ID (default: A)")
    p_mut.add_argument("--mut", default=None)
    p_mut.add_argument("--list", default=None, help="File with one mutation per line")
    p_mut.add_argument("--outdir-mut", default=None)
    p_mut.add_argument("--seed", type=int, default=None)
    p_mut.add_argument("--logfile", default=None)
    p_mut.add_argument("--verbose", action="store_true")

    # ---------------- Top-level Modern Commands (cmda / colabmda) ----------------
    # 1. colab
    p_colab = sub.add_parser("colab", help="Remote Google Colab VM provisioning & session manager")
    p_colab.add_argument(
        "colab_args", nargs=argparse.REMAINDER, help="Arguments passed to colab manager"
    )

    # 2. prep
    p_top_prep = sub.add_parser("prep", help="Unified model assembly & parameterization")
    p_top_prep.add_argument(
        "--state", default=None, help="State tag (e.g. atp_2mg_eb2, apo, atp_2mg)"
    )
    p_top_prep.add_argument("--variant", default="WT", help="Variant name (default: WT)")
    p_top_prep.add_argument("--protein", default="cdkl5", help="Protein name (default: cdkl5)")
    p_top_prep.add_argument("--receptor", default=None, help="Path to receptor PDB")
    p_top_prep.add_argument(
        "--template", default=None, help="Template PDB path or 4-letter ID (e.g. 3X2V)"
    )
    p_top_prep.add_argument("--ligand", default=None, help="Ligand name or path (e.g. ATP, ADP)")
    p_top_prep.add_argument("--ions", default=None, help="Ion spec (e.g. MG:2, MN:2)")
    p_top_prep.add_argument("--phos", default=None, help="Phosphorylation residue(s) (e.g. Y171)")
    p_top_prep.add_argument(
        "--peptide-seq", default=None, help="Substrate peptide sequence (e.g. KRPSSAPSG)"
    )
    p_top_prep.add_argument(
        "--peptide-name", default=None, help="Substrate peptide name (e.g. eb2, eb1)"
    )
    p_top_prep.add_argument("--root", "--project", dest="root", default=None, help=ROOT_HELP)
    p_top_prep.add_argument(
        "--no-amber", action="store_true", help="Skip automated tleap topology generation"
    )

    # 3. stage
    p_top_stage = sub.add_parser(
        "stage", help="Stage standardized directory trees & decoupled scripts"
    )
    p_top_stage.add_argument(
        "--state", default="atp_2mg_eb2", help="State tag (default: atp_2mg_eb2)"
    )
    p_top_stage.add_argument(
        "--example", choices=["trpcage", "fkbp", "hras"], help="Install a complete runnable example"
    )
    p_top_stage.add_argument(
        "--scripts", action="store_true", help="Optionally generate legacy batch wrappers"
    )
    p_top_stage.add_argument("--variant", default="WT", help="Variant name (default: WT)")
    p_top_stage.add_argument("--protein", default="cdkl5", help="Protein name (default: cdkl5)")
    p_top_stage.add_argument(
        "--replicas", default="1,2,3", help="Replicas to stage (default: 1,2,3)"
    )
    p_top_stage.add_argument("--root", "--project", dest="root", default=None, help=ROOT_HELP)

    # 4. sim
    p_top_sim = sub.add_parser("sim", help="Execute simulation stages (em, nvt, npt, qc, md)")
    p_top_sim.add_argument("--state", required=True, help="State tag (e.g. atp_2mg_eb2, apo)")
    p_top_sim.add_argument("--variant", default="WT", help="Variant name (default: WT)")
    p_top_sim.add_argument(
        "--stage",
        default="em,nvt,npt,qc",
        help="Comma-separated stage list (default: em,nvt,npt,qc)",
    )
    p_top_sim.add_argument(
        "--engine",
        choices=["openmm", "amber"],
        default="openmm",
        help="MD engine (default: openmm)",
    )
    p_top_sim.add_argument(
        "--replica",
        "--replicas",
        dest="replica",
        default="r1",
        help="Production replica(s) (comma-separated, e.g. r1 or r1,r2,r3; default: r1)",
    )
    p_top_sim.add_argument(
        "--total-ns", type=float, default=100.0, help="Production duration in ns (default: 100.0)"
    )
    p_top_sim.add_argument(
        "--traj-interval",
        type=float,
        default=10.0,
        help="Trajectory save interval in ps (default: 10.0)",
    )
    p_top_sim.add_argument(
        "--steps", type=int, default=None, help="Override step count for testing"
    )
    p_top_sim.add_argument("--checkpoint-ps", type=float, default=1000.0)
    p_top_sim.add_argument(
        "--equil-time", type=float, default=100.0, help="Equilibration duration in ps"
    )
    p_top_sim.add_argument(
        "--protocol",
        choices=["staged", "quick", "default", "varmdyn"],
        default="staged",
        help="Equilibration protocol: staged (rigorous 24-stage; default) or quick (smoke-test)",
    )
    p_top_sim.add_argument(
        "--software-test",
        action="store_true",
        help="Explicit short test; permits quick protocol and reports physical QC without a research gate",
    )
    p_top_sim.add_argument("--root", "--project", dest="root", default=None, help=ROOT_HELP)

    # 5. merge
    p_top_merge = sub.add_parser("merge", help="Autoimage, wrap PBC, and strip trajectories")
    p_top_merge.add_argument("--state", required=True, help="State tag")
    p_top_merge.add_argument("--variant", default="WT", help="Variant name")
    p_top_merge.add_argument(
        "--replicas",
        "--replica",
        dest="replicas",
        default="r1,r2,r3",
        help="Replicas to merge (default: r1,r2,r3)",
    )
    p_top_merge.add_argument("--root", "--project", dest="root", default=None, help=ROOT_HELP)

    # 6. analyze (alias: analysis)
    p_top_ana = sub.add_parser(
        "analyze",
        aliases=["analysis"],
        help="Calculate RMSD, Rg and RMSF, then render saved tables",
    )
    p_top_ana.add_argument("--state", required=True, help="State tag")
    p_top_ana.add_argument("--variant", default="WT", help="Variant name")
    p_top_ana.add_argument(
        "--metrics",
        default="rmsd,rg,rmsf",
        help="Compatibility option; current analysis calculates RMSD, Rg and RMSF",
    )
    p_top_ana.add_argument("--root", "--project", dest="root", default=None, help=ROOT_HELP)

    # 7. compare
    p_top_comp = sub.add_parser(
        "compare", help="Report a comparison request (plotting: openmm compare --series)"
    )
    p_top_comp.add_argument("--state", required=True, help="State tag")
    p_top_comp.add_argument(
        "--variants", required=True, help="Comma-separated variants list (e.g. WT,R59W)"
    )
    p_top_comp.add_argument("--metrics", default="rmsf,rmsd", help="Metrics to compare")
    p_top_comp.add_argument("--root", "--project", dest="root", default=None, help=ROOT_HELP)

    # 8. view
    p_top_view = sub.add_parser("view", help="View trajectory in PyMOL")
    p_top_view.add_argument(
        "--pdb-dir",
        default=None,
        help="Directory containing the merged files (default: current directory)",
    )
    p_top_view.add_argument(
        "-t", "--topology", default=None, help="Topology file (default: prod_full.pdb)"
    )
    p_top_view.add_argument(
        "-x", "--trajectory", default=None, help="Trajectory file (default: prod_full.dcd)"
    )
    p_top_view.add_argument(
        "-r", "--resi", type=int, default=12, help="Highlight residue index (default: 12)"
    )

    # 9. snapshots
    p_top_snap = sub.add_parser("snapshots", help="Generate comparative snapshot grid")
    p_top_snap.add_argument(
        "-c", "--config", required=True, help="Path to custom JSON configuration file"
    )
    p_top_snap.add_argument(
        "-o",
        "--output",
        default="data/analysis/snapshots/master_3x8_snapshots.png",
        help="Path to output image grid file",
    )
    p_top_snap.add_argument(
        "--temp-dir", default="scratch/master_grid_3x8", help="Temporary folder for panel rendering"
    )
    p_top_snap.add_argument("--font-path", default=None, help="Custom TrueType font path")

    args = p.parse_args()

    if args.tool in {"fetch", "model", "check", "assemble"}:
        from colabmda.core import workflow

        raise SystemExit(getattr(workflow, args.tool)(args))

    if args.tool == "colab":
        from colabmda.remote.colab import colab_manager_cli

        raise SystemExit(colab_manager_cli(args.colab_args))

    elif args.tool == "prep":
        raise SystemExit(_dispatch_prep(args))

    elif args.tool == "stage":
        raise SystemExit(_dispatch_stage(args))

    elif args.tool == "sim":
        raise SystemExit(_dispatch_sim(args))

    elif args.tool == "merge":
        raise SystemExit(_dispatch_merge(args))

    elif args.tool in ["analyze", "analysis"]:
        raise SystemExit(_dispatch_analyze(args))

    elif args.tool == "compare":
        raise SystemExit(_dispatch_compare(args))

    elif args.tool == "view":
        from colabmda.openmm.commands import openmm_view

        raise SystemExit(
            openmm_view(args.pdb_dir, args.topology, args.trajectory, resi=args.resi) or 0
        )

    elif args.tool == "snapshots":
        from colabmda.openmm.commands import openmm_snapshots

        raise SystemExit(
            openmm_snapshots(
                config_path=args.config,
                output=args.output,
                temp_dir=args.temp_dir,
                font_path=args.font_path,
            )
            or 0
        )

    if args.tool == "openmm":
        from colabmda.openmm.commands import (
            openmm_analysis,
            openmm_check_equil,
            openmm_compare,
            openmm_em,
            openmm_md,
            openmm_merge,
            openmm_npt,
            openmm_nvt,
            openmm_prep_from_file,
            openmm_prep_from_pdbid,
            openmm_qc,
            openmm_snapshots,
            openmm_status,
            openmm_view,
        )

        if args.cmd == "prep":
            if args.pdb_id:
                root = _resolve_root(args.drive, args.root)
                if root:
                    _ensure_dir(root)
                openmm_prep_from_pdbid(
                    args.pdb_id, root_dir=root, ph=args.ph, sync_dir=args.sync_dir
                )
            else:
                name = args.name or Path(args.pdb_file).stem
                root = _resolve_root(args.drive, args.root)
                if root:
                    _ensure_dir(root)
                outdir = args.outdir or (str(Path(root) / name / "prep") if root else name)
                openmm_prep_from_file(
                    args.pdb_file, outdir, pdbid=name, ph=args.ph, sync_dir=args.sync_dir
                )

        elif args.cmd == "stage":
            name = args.name or Path(args.pdb_file).stem
            root_path = Path(
                _resolve_root(getattr(args, "drive", False), getattr(args, "root", None))
            ).resolve()
            _ensure_dir(str(root_path))

            # Priority: data/simulation > data/sim > simulation > sim > simulations (legacy)
            if (root_path / "data" / "simulation").exists():
                sim_base = root_path / "data" / "simulation"
            elif (root_path / "data" / "sim").exists():
                sim_base = root_path / "data" / "sim"
            elif (root_path / "simulation").exists():
                sim_base = root_path / "simulation"
            elif (root_path / "sim").exists():
                sim_base = root_path / "sim"
            elif (root_path / "simulations").exists():
                sim_base = root_path / "simulations"
            else:
                sim_base = root_path / "data" / "simulation"

            sim_dir = sim_base / name
            sim_dir.mkdir(parents=True, exist_ok=True)

            target_clean = sim_dir / f"{name}_cleaned.pdb"
            in_pdb = Path(args.pdb_file).resolve()

            if in_pdb.exists():
                shutil.copy2(in_pdb, target_clean)
                summary_json = in_pdb.parent / "protonation_summary.json"
                if summary_json.exists():
                    shutil.copy2(summary_json, sim_dir / "protonation_summary.json")
            else:
                openmm_prep_from_file(str(in_pdb), str(sim_dir), pdbid=name, ph=args.ph)

            # Always stage the shared equilibration hub 'eq'
            eq_dir = sim_dir / "eq"
            eq_dir.mkdir(parents=True, exist_ok=True)
            if target_clean.exists():
                shutil.copy2(target_clean, eq_dir / f"{name}_cleaned.pdb")
            if (sim_dir / "protonation_summary.json").exists():
                shutil.copy2(
                    sim_dir / "protonation_summary.json", eq_dir / "protonation_summary.json"
                )

            if args.replica:
                rep_dir = sim_dir / args.replica
                rep_dir.mkdir(parents=True, exist_ok=True)
                if target_clean.exists():
                    shutil.copy2(target_clean, rep_dir / f"{name}_cleaned.pdb")
                if (sim_dir / "protonation_summary.json").exists():
                    shutil.copy2(
                        sim_dir / "protonation_summary.json", rep_dir / "protonation_summary.json"
                    )

            print(f"[INFO] Staged simulation folder: {sim_dir}")

        elif args.cmd == "run":
            if args.pdb_id:
                root = _resolve_root(args.drive, args.root)
                if root:
                    _ensure_dir(root)
                workdir = (
                    str(Path(root) / args.pdb_id / "run")
                    if root
                    else str(Path(args.pdb_id) / "run")
                )
                name = args.name or args.pdb_id
                if root:
                    _prepare_run_inputs(root=root, pdbid=args.pdb_id, name=name)
            elif args.workdir:
                workdir = str(Path(args.workdir).resolve())
                name = args.name or _guess_pdbid_from_workdir(workdir)
            else:
                root = Path(_resolve_root(args.drive, args.root)).resolve()
                name = args.name
                cwd = Path.cwd().resolve()

                print(f"[DEBUG] root: {root}")
                print(f"[DEBUG] cwd:  {cwd}")
                print(f"[DEBUG] name: {name}")

                # Search strategy: prioritize data/simulation > data/sim > simulation > sim > simulations
                search_paths = [
                    cwd / "data" / "simulation" / name if name else None,
                    cwd / "data" / "sim" / name if name else None,
                    cwd / "simulation" / name if name else None,
                    cwd / "sim" / name if name else None,
                    cwd / "simulations" / name if name else None,
                    cwd / name if name else None,
                    cwd if (name and name in str(cwd)) else None,  # Check if we are already inside
                    root / "data" / "simulation" / name if name else None,
                    root / "data" / "sim" / name if name else None,
                    root / "simulation" / name if name else None,
                    root / "sim" / name if name else None,
                    root / "simulations" / name if name else None,
                    root / name if name else None,
                ]

                workdir = None
                for p in search_paths:
                    if p:
                        exists = p.exists()
                        is_dir = p.is_dir() if exists else False
                        print(f"[DEBUG] Checking: {p} (exists={exists}, is_dir={is_dir})")
                        if exists and is_dir:
                            # If we found the replica folder itself, use its parent as the simulation base
                            if p.name == args.replica:
                                workdir = p.parent
                            else:
                                workdir = p
                            break

                if workdir is None:
                    # If we are in a folder named 'name', use it
                    if name and cwd.name == name:
                        workdir = cwd
                    else:
                        print(f"[DEBUG] No folder found, defaulting to cwd: {cwd}")
                        workdir = cwd

                if not name:
                    name = _guess_pdbid_from_workdir(str(workdir))

            workdir = Path(workdir).resolve()
            if args.replica:
                # If we are already inside the replica folder, don't double-append
                if workdir.name != args.replica:
                    workdir = workdir / args.replica

                _ensure_dir(str(workdir))

                # Ensure cleaned PDB is available in the replica folder
                local_clean = workdir / f"{name}_cleaned.pdb"
                parent_clean = workdir.parent / f"{name}_cleaned.pdb"
                if not local_clean.exists() and parent_clean.exists():
                    shutil.copy2(parent_clean, local_clean)
                elif not local_clean.exists():
                    # Fallback: find any cleaned PDB in parent to infer name
                    guesses = list(workdir.parent.glob("*_cleaned.pdb"))
                    if guesses and not name:
                        name = guesses[0].name.replace("_cleaned.pdb", "")
                        shutil.copy2(guesses[0], local_clean)

            workdir_str = str(workdir.resolve())
            if not name:
                name = _guess_pdbid_from_workdir(workdir_str)

            if not name:
                raise SystemExit("ERROR: could not infer pdbid; please specify --name.")

            # Auto-discover cleaned starting structure if missing from workdir
            workdir_path = Path(workdir_str)
            clean_target = workdir_path / f"{name}_cleaned.pdb"
            if not clean_target.exists():
                candidates = [
                    workdir_path / "r1" / f"{name}_cleaned.pdb",
                    workdir_path / "eq" / f"{name}_cleaned.pdb",
                    workdir_path / "equil" / f"{name}_cleaned.pdb",
                    workdir_path / "prep" / f"{name}_cleaned.pdb",
                    workdir_path.parent / f"{name}_cleaned.pdb",
                    workdir_path.parent / "eq" / f"{name}_cleaned.pdb",
                    workdir_path.parent / "equil" / f"{name}_cleaned.pdb",
                ]
                glob_cands = list(workdir_path.glob("*/*_cleaned.pdb")) + list(
                    workdir_path.glob("*_cleaned.pdb")
                )
                str_cands = (
                    list(Path.cwd().glob(f"data/structure/**/{name}*/*_cleaned.pdb"))
                    + list(Path.cwd().glob(f"data/structures/**/{name}*/*_cleaned.pdb"))
                    + list(Path.cwd().glob(f"data/str/**/{name}*/*_cleaned.pdb"))
                    + list(Path.cwd().glob(f"structure/**/{name}*/*_cleaned.pdb"))
                    + list(Path.cwd().glob(f"str/**/{name}*/*_cleaned.pdb"))
                )
                for c in glob_cands + str_cands:
                    candidates.append(c)

                for c in candidates:
                    if c.exists() and c.is_file():
                        print(f"[INFO] Auto-discovered cleaned starting structure: {c}")
                        shutil.copy2(c, clean_target)
                        eq_dir = workdir_path / "eq"
                        if eq_dir.exists() and eq_dir.is_dir():
                            eq_clean = eq_dir / f"{name}_cleaned.pdb"
                            if not eq_clean.exists():
                                shutil.copy2(c, eq_clean)
                        r1_dir = workdir_path / "r1"
                        if r1_dir.exists() and r1_dir.is_dir():
                            r1_clean = r1_dir / f"{name}_cleaned.pdb"
                            if not r1_clean.exists():
                                shutil.copy2(c, r1_clean)
                        break

            # Copy ligand file to workdir if provided
            if args.ligand:
                ligand_path = Path(args.ligand)
                local_ligand = Path(workdir_str) / ligand_path.name
                if ligand_path.exists() and ligand_path.resolve() != local_ligand.resolve():
                    shutil.copy2(ligand_path, local_ligand)
                args.ligand = str(local_ligand.resolve())

            # MODULAR RUN: EM -> NVT -> NPT -> Check -> MD
            from colabmda.openmm.common import derive_replica_seed

            replica_seed = derive_replica_seed(args.seed, args.replica)

            # Skip equilibration for non-r1 replicas (e.g. r2, r3) as they inherit from r1/parent
            is_non_r1_replica = False
            if args.replica and args.replica != "r1":
                is_non_r1_replica = True

            if is_non_r1_replica:
                print(
                    f"[INFO] Replica '{args.replica}' detected (seed={replica_seed}). Skipping EM/NVT/NPT and inheriting equilibrated structure."
                )
            else:
                equil_protocol = _normalize_equil_protocol(args.equil_protocol)
                openmm_em(
                    workdir_str,
                    name,
                    ligand=args.ligand,
                    keep_mg=args.keep_mg,
                    amber_prmtop=args.amber_prmtop,
                    amber_inpcrd=args.amber_inpcrd,
                    padding_nm=args.padding_nm,
                    protein_ff=args.protein_ff,
                    water_model=args.water_model,
                    small_molecule_ff=args.small_molecule_ff,
                )
                openmm_nvt(
                    workdir_str, name, args.equil_time, seed=replica_seed, protocol=equil_protocol
                )
                openmm_npt(
                    workdir_str, name, args.equil_time, seed=replica_seed, protocol=equil_protocol
                )
                openmm_check_equil(workdir_str, strict=equil_protocol == "varmdyn")

            if getattr(args, "equil_only", False):
                print(
                    "✔ Shared equilibration complete. Stopping before production MD as requested by --equil-only."
                )
                return

            openmm_md(
                workdir=workdir_str,
                pdbid=name,
                total_ns=args.total_ns,
                traj_interval=args.traj_interval,
                checkpoint_ps=args.checkpoint_ps,
                sync_dir=args.sync_dir,
                seed=replica_seed,
            )

        elif args.cmd == "merge":
            root = _resolve_root(args.drive, args.root)
            name = getattr(args, "name", None)
            replica = getattr(args, "replica", "r1") or "r1"
            if name:
                sim_base = _simulation_workdir(root or ".", name)
                pdbid_dir = str((sim_base / replica).resolve())
            elif args.pdb_id and root:
                pdbid_dir = str(Path(root) / args.pdb_id / "run")
            elif args.pdb_id:
                pdbid_dir = str(Path(args.pdb_id) / "run")
            elif args.pdb_dir:
                pdb_path = Path(args.pdb_dir)
                if not pdb_path.is_absolute() and root and (Path(root) / pdb_path).exists():
                    pdbid_dir = str((Path(root) / pdb_path).resolve())
                else:
                    pdbid_dir = args.pdb_dir
            else:
                pdbid_dir = "."

            # Resolve output paths: bare filenames go inside pdbid_dir, not CWD
            _pdbid_path = Path(pdbid_dir).resolve()
            out_traj = (
                args.out_traj
                if Path(args.out_traj).is_absolute()
                else str(_pdbid_path / args.out_traj)
            )
            out_log = (
                args.out_log
                if Path(args.out_log).is_absolute()
                else str(_pdbid_path / args.out_log)
            )

            openmm_merge(
                pdbid_dir,
                args.topology,
                out_traj,
                out_log,
                stride=args.stride,
                center=args.center,
                wrap=args.wrap,
                mda=args.mda,
                selection=args.selection,
                ca_only=args.ca_only,
                protein_only=args.protein_only,
            )

        elif args.cmd == "analysis":
            root = _resolve_root(args.drive, args.root)
            name = getattr(args, "name", None)
            replica = getattr(args, "replica", None)
            if name:
                sim_base = _simulation_workdir(root or ".", name)
                sim_dir = (
                    str((sim_base / replica).resolve()) if replica else str(sim_base.resolve())
                )
                if args.outdir:
                    out_base = Path(args.outdir)
                elif root:
                    out_base = _analysis_outdir(root, name)
                    if replica:
                        out_base = out_base / replica
                else:
                    out_base = Path("data/analysis") / name
                    if replica:
                        out_base = out_base / replica
            elif args.pdb_id and root:
                sim_path = Path(root) / "data" / "sim" / args.pdb_id
                if not sim_path.exists():
                    sim_path = Path(root) / "sim" / args.pdb_id
                if not sim_path.exists():
                    sim_path = Path(root) / "simulations" / args.pdb_id
                sim_dir = str(sim_path)
                out_base = Path(root) / "data" / "analysis" / args.pdb_id
            elif args.pdb_id:
                sim_path = Path("data") / "sim" / args.pdb_id
                if not sim_path.exists():
                    sim_path = Path("sim") / args.pdb_id
                if not sim_path.exists():
                    sim_path = Path("simulations") / args.pdb_id
                sim_dir = str(sim_path)
                out_base = Path("data") / "analysis" / args.pdb_id
            elif args.pdb_dir:
                pdb_path = Path(args.pdb_dir)
                if not pdb_path.is_absolute() and root and (Path(root) / pdb_path).exists():
                    sim_dir = str((Path(root) / pdb_path).resolve())
                else:
                    sim_dir = args.pdb_dir
                parent_name = (
                    Path(sim_dir).parent.name
                    if Path(sim_dir).name.startswith("r")
                    else Path(sim_dir).name
                )
                if args.outdir:
                    out_base = Path(args.outdir)
                elif root:
                    out_base = _analysis_outdir(root, parent_name)
                else:
                    out_base = Path("data/analysis") / parent_name
            else:
                sim_dir = "."
                out_base = Path(args.outdir) if args.outdir else Path("data/analysis/current")

            openmm_analysis(
                sim_dir,
                args.topology,
                args.trajectory,
                args.interval,
                str(out_base),
                reference=args.reference,
            )

        elif args.cmd == "compare":
            openmm_compare(args.series, args.outdir)

        elif args.cmd == "view":
            openmm_view(args.pdb_dir, args.topology, args.trajectory, resi=args.resi)

        elif args.cmd == "snapshots":
            openmm_snapshots(
                config_path=args.config,
                output=args.output,
                temp_dir=args.temp_dir,
                font_path=args.font_path,
            )

        elif args.cmd == "status":
            root = _resolve_root(args.drive, args.root)
            name = getattr(args, "name", None)
            replica = getattr(args, "replica", "r1") or "r1"
            if name:
                sim_base = _simulation_workdir(root or ".", name)
                pdbid_dir = str((sim_base / replica).resolve())
            elif args.pdb_id and root:
                pdbid_dir = str(Path(root) / args.pdb_id / "run")
            elif args.pdb_id:
                pdbid_dir = str(Path(args.pdb_id) / "run")
            elif args.pdb_dir:
                pdb_path = Path(args.pdb_dir)
                if not pdb_path.is_absolute() and root and (Path(root) / pdb_path).exists():
                    pdbid_dir = str((Path(root) / pdb_path).resolve())
                else:
                    pdbid_dir = args.pdb_dir
            else:
                pdbid_dir = "."
            openmm_status(pdbid_dir)

        elif args.cmd in ["em", "nvt", "npt", "qc", "check-equil", "md"]:
            name = args.name
            if args.workdir:
                workdir = str(Path(args.workdir).resolve())
            else:
                root_val = args.root or _resolve_root(getattr(args, "drive", False), None)
                sim_base = _simulation_workdir(root_val, name)
                if args.cmd in ["em", "nvt", "npt", "qc", "check-equil"]:
                    if (sim_base / "eq").exists():
                        workdir_path = sim_base / "eq"
                    elif (sim_base / "equil").exists():
                        workdir_path = sim_base / "equil"
                    else:
                        workdir_path = sim_base / "eq"
                elif args.cmd == "md":
                    rep = getattr(args, "replica", "r1") or "r1"
                    workdir_path = sim_base / rep
                else:
                    workdir_path = sim_base
                workdir_path.mkdir(parents=True, exist_ok=True)
                workdir = str(workdir_path.resolve())

                # If cleaned PDB exists in sim_base, eq, or equil, ensure it's available in workdir
                workdir_clean = workdir_path / f"{name}_cleaned.pdb"
                if not workdir_clean.exists():
                    for src in [
                        sim_base / f"{name}_cleaned.pdb",
                        sim_base / "eq" / f"{name}_cleaned.pdb",
                        sim_base / "equil" / f"{name}_cleaned.pdb",
                    ]:
                        if src.exists():
                            shutil.copy2(src, workdir_clean)
                            break
            if args.cmd == "em":
                if args.ligand:
                    ligand_path = Path(args.ligand)
                    local_ligand = Path(workdir) / ligand_path.name
                    if ligand_path.exists() and ligand_path.resolve() != local_ligand.resolve():
                        shutil.copy2(ligand_path, local_ligand)
                    args.ligand = str(local_ligand.resolve())
                openmm_em(
                    workdir,
                    name,
                    ligand=args.ligand,
                    keep_mg=args.keep_mg,
                    amber_prmtop=args.amber_prmtop,
                    amber_inpcrd=args.amber_inpcrd,
                    padding_nm=args.padding_nm,
                    protein_ff=args.protein_ff,
                    water_model=args.water_model,
                    small_molecule_ff=args.small_molecule_ff,
                )
            elif args.cmd == "nvt":
                openmm_nvt(
                    workdir,
                    name,
                    args.equil_time,
                    args.seed,
                    _normalize_equil_protocol(args.equil_protocol),
                )
            elif args.cmd == "npt":
                openmm_npt(
                    workdir,
                    name,
                    args.equil_time,
                    seed=args.seed,
                    protocol=_normalize_equil_protocol(args.equil_protocol),
                )
            elif args.cmd in ["qc", "check-equil"]:
                openmm_qc(workdir, strict=not args.warn_only)
            elif args.cmd == "md":
                openmm_md(
                    workdir,
                    name,
                    args.total_ns,
                    args.traj_interval,
                    args.checkpoint_ps,
                    args.sync_dir,
                    seed=args.seed,
                )

    elif args.tool == "modeller":
        from colabmda.modeller.commands import (
            modeller_build,
            modeller_mutate,
        )

        if args.cmd == "build":
            modeller_build(args)
        elif args.cmd == "mutate":
            if not (args.mut or args.list):
                raise SystemExit("ERROR: mutate requires --mut or --list")
            modeller_mutate(args)


if __name__ == "__main__":
    main()
