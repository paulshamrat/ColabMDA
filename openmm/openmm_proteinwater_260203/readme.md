# 260203 new approach openmm
nvidia-smi


mkdir -p /content/md_project
mkdir -p /content/work/4ldj_wt
mkdir -p /content/work/4ldj_G12C
mkdir -p /content/drive/MyDrive/pdbs
mkdir -p /content/drive/MyDrive/openmm_runs/4ldj_wt
mkdir -p /content/drive/MyDrive/openmm_runs/4ldj_G12C




ls -lh /content/drive/MyDrive/pdbs/4LDJ.pdb


/content/md_project/
  pdbfixer_clean_fromfile.py
  openmm_proteinwater_colab.py
  openmm_trajmerge.py
  openmm_trajanalysis.py


mv -v openmm_trajmerge.py openmm_trajanalysis.py /content/md_project/ 2>/dev/null || true


python3 /content/md_project/pdbfixer_clean_fromfile.py \
  --in /content/drive/MyDrive/pdbs/4LDJ.pdb \
  --outdir /content/work/4ldj_wt \
  --pdbid 4ldj



python3 /content/md_project/pdbfixer_clean_fromfile.py \
  --in /content/drive/MyDrive/pdbs/4LDJ.pdb \
  --outdir /content/work/4ldj_G12C \
  --pdbid 4ldj


# 3B) Start WT production MD (resume-safe)
python3 /content/md_project/openmm_proteinwater_colab.py /content/work/4ldj_wt \
  --pdbid 4ldj \
  --total-ns 10 \
  --traj-interval 10 \
  --equil-time 1000 \
  --checkpoint-ps 1000 \
  --sync-dir /content/drive/MyDrive/openmm_runs/4ldj_wt


ls -lh /content/work/4ldj_wt | tail
