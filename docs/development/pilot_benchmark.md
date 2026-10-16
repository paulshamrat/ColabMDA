# Developer & CI Benchmark Pilot Run: 1 ns Equilibration + 1 ns Production

> 🧑‍💻 **Internal Developer & CI Protocol**  
> This document records the standardized 1 ns Equilibration + 1 ns Production pilot execution protocol used for developer testing, automated benchmarking, active-site coordination validation, and CI verification. It is isolated from the main user-facing README to keep user workflows focused on 100 ns production campaigns.

---

## 1. Pilot Run Benchmark Parameters

* **Equilibration Duration:** 1 ns total (`--equil-time 1000` = 1,000 ps across NVT and NPT).
* **Production Duration:** 1 ns total (`--total-ns 1.0` = 500,000 steps at 2.0 fs timestep).
* **Trajectory Frame Interval:** `--traj-interval 10` (10 ps = 100 frames).
* **Expected Walltime:** ~1.5–2.5 minutes per replica on modern NVIDIA GPUs (T4 / V100 / A100 / RTX).

---

## 2. Google Colab End-to-End Pilot Execution

```bash
# ==============================================================================
# Google Colab Pilot Execution (WT & L119R across Apo and Holo atp_mg)
# ==============================================================================
PROJECT_DIR="/content/drive/MyDrive/projects/test-project"

# 1. Stage project directories and scripts
cmda stage --state apo --variant WT --root "${PROJECT_DIR}"
cmda stage --state apo --variant L119R --root "${PROJECT_DIR}"
cmda stage --state atp_mg --variant WT --root "${PROJECT_DIR}"
cmda stage --state atp_mg --variant L119R --root "${PROJECT_DIR}"

# 2. Build parameter topologies (Amber ff19SB / Hu 2024 / OPC 4-point water)
cmda prep --state apo --variant WT --root "${PROJECT_DIR}"
cmda prep --state apo --variant L119R --root "${PROJECT_DIR}"
cmda prep --state atp_mg --variant WT --root "${PROJECT_DIR}"
cmda prep --state atp_mg --variant L119R --root "${PROJECT_DIR}"

# 3. 1 ns Equilibration Hub (EM -> NVT -> NPT -> QC)
cmda sim --state apo --variant WT --stage eq --equil-time 1000 --root "${PROJECT_DIR}"
cmda sim --state apo --variant L119R --stage eq --equil-time 1000 --root "${PROJECT_DIR}"
cmda sim --state atp_mg --variant WT --stage eq --equil-time 1000 --root "${PROJECT_DIR}"
cmda sim --state atp_mg --variant L119R --stage eq --equil-time 1000 --root "${PROJECT_DIR}"

# 4. 1 ns Production MD Pilot Replicas (r1)
cmda sim --state apo --variant WT --stage md --replica r1 --total-ns 1.0 --traj-interval 10 --root "${PROJECT_DIR}"
cmda sim --state apo --variant L119R --stage md --replica r1 --total-ns 1.0 --traj-interval 10 --root "${PROJECT_DIR}"
cmda sim --state atp_mg --variant WT --stage md --replica r1 --total-ns 1.0 --traj-interval 10 --root "${PROJECT_DIR}"
cmda sim --state atp_mg --variant L119R --stage md --replica r1 --total-ns 1.0 --traj-interval 10 --root "${PROJECT_DIR}"

# 5. Automated Autoimaging & Trajectory Post-Processing
cmda merge --state apo --variant WT --replicas r1 --root "${PROJECT_DIR}"
cmda merge --state apo --variant L119R --replicas r1 --root "${PROJECT_DIR}"
cmda merge --state atp_mg --variant WT --replicas r1 --root "${PROJECT_DIR}"
cmda merge --state atp_mg --variant L119R --replicas r1 --root "${PROJECT_DIR}"

# 6. Biophysical Analysis & Comparative Overlays
cmda analyze --state apo --variant WT --root "${PROJECT_DIR}"
cmda analyze --state apo --variant L119R --root "${PROJECT_DIR}"
cmda analyze --state atp_mg --variant WT --root "${PROJECT_DIR}"
cmda analyze --state atp_mg --variant L119R --root "${PROJECT_DIR}"
cmda compare --state apo --variants WT,L119R --root "${PROJECT_DIR}"
cmda compare --state atp_mg --variants WT,L119R --root "${PROJECT_DIR}"
```

---

## 3. HPC Cluster Scratch End-to-End Pilot Execution (SLURM)

```bash
# ==============================================================================
# HPC Cluster Scratch Pilot Execution (SLURM)
# ==============================================================================
export SCRATCH_PROJECT="/scratch/$USER/projects/test-project"

# 1. Stage project directories and generate SLURM runner scripts
cmda stage --state apo --variant WT --root "${SCRATCH_PROJECT}"
cmda stage --state apo --variant L119R --root "${SCRATCH_PROJECT}"
cmda stage --state atp_mg --variant WT --root "${SCRATCH_PROJECT}"
cmda stage --state atp_mg --variant L119R --root "${SCRATCH_PROJECT}"

# 2. Build parameter topologies
cmda prep --state apo --variant WT --root "${SCRATCH_PROJECT}"
cmda prep --state apo --variant L119R --root "${SCRATCH_PROJECT}"
cmda prep --state atp_mg --variant WT --root "${SCRATCH_PROJECT}"
cmda prep --state atp_mg --variant L119R --root "${SCRATCH_PROJECT}"

# 3. Interactive Pilot Execution on GPU Node (salloc / qsub)
# Request GPU node: salloc --nodes=1 --ntasks=1 --cpus-per-task=8 --gpus=1 --mem=32gb --time=01:00:00
cmda sim --state apo --variant WT --stage eq --equil-time 1000 --root "${SCRATCH_PROJECT}"
cmda sim --state apo --variant WT --stage md --replica r1 --total-ns 1.0 --traj-interval 10 --root "${SCRATCH_PROJECT}"
cmda sim --state atp_mg --variant WT --stage eq --equil-time 1000 --root "${SCRATCH_PROJECT}"
cmda sim --state atp_mg --variant WT --stage md --replica r1 --total-ns 1.0 --traj-interval 10 --root "${SCRATCH_PROJECT}"

# 4. Or Submit via Chained SLURM Batch Helpers (EQ -> Array Prod)
cd "${SCRATCH_PROJECT}"
bash scripts/submit_chain.sh
```
