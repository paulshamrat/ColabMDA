# Integration validation

The main README owns the complete Trp-cage user guide. Templates live under
examples/trpcage/scripts and are installed with cmda stage --example trpcage.
The private research/projects/colabmda-example project tests those same files.

Run unit/regression tests on Divergent. The project's tests/gpu.py runs a bounded
Colab GPU integration test: prepare, quick 100 ps NVT/NPT, 100 ps smoke, SIGINT
during the next chunk, resume to 1 ns, merge, analyze and check. Require the first
chunk hash unchanged, 100 frames with continuous times, finite numerical outputs,
20 C-alpha atoms, author labels 1–20, reviewed PNGs and completion markers.

Record tool/project/template/input versions and retain physical QC separately.
Trp-cage is one apo system, not validation of CDKL5 mutations, ATP/Mg or the full
staged protocol. Keep historical provenance unchanged during migrations. Heavy
outputs stay on Drive. Stop the dedicated session after verified retrieval.
