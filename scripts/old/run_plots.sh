#!/bin/bash
# Waits for shower_batch screen to finish, then re-runs all plots.
# Launch: screen -dmS shower_plots bash run_plots.sh
# Check:  tail -f /tmp/shower_plots.log

LOG="/tmp/shower_plots.log"
exec > >(tee -a "$LOG") 2>&1

echo "=== shower_plots waiting for shower_batch at $(date) ==="

# Wait until the shower_batch screen session disappears
while screen -ls | grep -q shower_batch; do
    sleep 60
done

echo "=== shower_batch done, starting plots at $(date) ==="

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source "$ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh" --quiet
asetup Athena,25.0.40

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

bash run.sh --batch --plot-only

echo "=== shower_plots finished at $(date) ==="
