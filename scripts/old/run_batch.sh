#!/bin/bash
# Self-contained batch launcher — sets up ATLAS env then runs all scenarios.
# Launch detached: screen -dmS shower_batch bash run_batch.sh
# Check progress: screen -r shower_batch   (Ctrl+A D to detach again)
# Or watch log:   tail -f /tmp/shower_batch.log

LOG="/tmp/shower_batch.log"
exec > >(tee -a "$LOG") 2>&1

echo "=== shower_batch started at $(date) ==="

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source "$ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh" --quiet
asetup Athena,25.0.40

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

bash run.sh --batch

echo "=== shower_batch finished at $(date) ==="
