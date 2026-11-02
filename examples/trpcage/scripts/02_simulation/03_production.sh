#!/usr/bin/env bash
set -euo pipefail
source "$(dirname -- "${BASH_SOURCE[0]}")/common.sh"
# Change duration here, or set COLABMDA_EXAMPLE_NS. Keep run ID/settings to resume.
PRODUCTION_NS=${COLABMDA_EXAMPLE_NS:-1}
CHECKPOINT_PS=100
FRAME_PS=10
test -s "$SYSTEM_DIR/02_eq/EQ_DONE"
export PRODUCTION_NS CHECKPOINT_PS FRAME_PS
python - <<'PYTHON'
import json,os
from pathlib import Path
p=Path(os.environ["SYSTEM_DIR"]).parents[1]/"production.json"
record={"requested_ns":float(os.environ["PRODUCTION_NS"]),"frame_ps":float(os.environ["FRAME_PS"]),"checkpoint_ps":float(os.environ["CHECKPOINT_PS"])}
if record["requested_ns"]<=0 or abs(record["requested_ns"]*10-round(record["requested_ns"]*10))>1e-8:
    raise ValueError("Use a positive duration divisible by 0.1 ns")
history=json.loads(p.read_text()) if p.exists() else []
if not history or history[-1]!=record: history.append(record)
p.write_text(json.dumps(history,indent=2)+"\n")
PYTHON
verify_run
cmda sim --software-test --root "$PROJECT_ROOT" --state apo --variant WT --stage md --replica r1 --total-ns "$PRODUCTION_NS" --checkpoint-ps "$CHECKPOINT_PS" --traj-interval "$FRAME_PS"
test -s "$SYSTEM_DIR/03_md/r1/prod.chk"
printf 'Production command completed for %s ns; frame validation follows\n' "$PRODUCTION_NS" > "$SYSTEM_DIR/03_md/r1/MD_DONE"
