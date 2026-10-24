# Interrupt and resume

Complete the basic example's preparation and equilibration first. Start production
with the same project root, state, variant and `COLABMDA_RUN_ID`. After at least one
100 ps chunk has completed, press Ctrl+C while the following chunk is running.
The interrupted command must report failure rather than claim completion.

Rerun the same `simulate` command. It loads `prod.chk` and replays only the
unfinished chunk. Retain the original input, force field, GPU/software environment,
frame interval and checkpoint interval. Never delete completed chunks to force a
resume, and do not assume binary checkpoints work across changed hardware/software.

To extend the basic 100 ps smoke test to 1 ns, retain the run ID and execute:

```bash
COLABMDA_EXAMPLE_NS=1 bash /path/to/ColabMDA/examples/basic/run.sh /path/to/project /path/to/protein.pdb simulate
```

Then run the `analyze` stage. Verify that completed chunk hashes are unchanged and
merged times are continuous, without duplicate or missing frames. The private
`test-colabmda` integration project automates SIGINT, hash comparison and frame-time
checks. A new scientific protocol or input requires a new run ID.
