# Colab GPU execution

On Divergent, create a dedicated GPU session, mount Drive and deploy committed
ColabMDA code. Complete authentication when requested. Run these operations
serially; wait for each to finish.

```bash
colab new -s example-colabmda --gpu T4
colab drivemount --session example-colabmda
cmda colab setup -s example-colabmda --local-source
colab console -s example-colabmda
```

Inside the Colab console, activate the installed environment and follow the
[basic example](../basic/README.md). The deployed tool is `/content/ColabMDA`;
choose a project root under `/content/drive/MyDrive/research/projects/` and supply
a reviewed input already stored on Drive. Keep the run in the foreground with
live output. Session termination remains possible; use committed checkpoints to
recover. A replacement session needs Drive authentication and the recorded tool
and environment before resuming.

You can also dispatch foreground commands from Divergent with
`cmda colab exec -s example-colabmda --timeout 21600 'YOUR_COMMAND'`.
Keep trajectories and checkpoints on Drive. Retrieve and verify compact reports,
CSV tables, figures and logs. Stop the dedicated session after review:

```bash
colab stop -s example-colabmda
```

Do not stop a session that still owns unrelated user work.
