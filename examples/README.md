# Examples

Start with [basic](basic/README.md) for one prepared protein. Follow
[restart](restart/README.md) to resume the same production run, or
[colab](colab/README.md) to execute on a remote GPU with Drive storage.

Examples call the installed `cmda` CLI. They contain no simulation engine,
private research inputs, trajectories or checkpoints. Write all runtime products
into a separate project with an explicit root and run ID. The basic runner can
be called one stage at a time, allowing integration projects to check outputs
between stages instead of maintaining a second copy of the command sequence.

Quick equilibration and a short trajectory test the software; they do not establish
biological convergence. See the detailed pilot acceptance criteria in
[the developer guide](../docs/development/pilot_benchmark.md).
