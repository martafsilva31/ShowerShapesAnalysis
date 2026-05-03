#!/bin/bash
# run_all_channels.sh — Run full pipeline for all channels
#
# Usage:
#   ./run_all_channels.sh              # all 3 channels, all events
#   ./run_all_channels.sh 100000       # all channels, maxEvents=100000

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh --quiet 2>/dev/null
lsetup "root 6.28.04-x86_64-el9-gcc13-opt" --quiet 2>/dev/null

cd /project/atlas/users/mfernand/QT/ShowerShapes/ShowerShapesAnalysis/scripts/data_mc

MAX_EVENTS="${1:--1}"

echo "Starting at $(date)"
echo "Max events: $MAX_EVENTS"

for ch in Zeeg Zmumug combined; do
    echo ""
    echo "================================================================"
    echo "  Running $ch at $(date)"
    echo "================================================================"
    bash run_pipeline.sh "$ch" "$MAX_EVENTS"
    echo "=== Done: $ch ==="
done

echo ""
echo "=== Fudge factor plots ==="
root -l -b -q "plot_fudge_binned.C(${MAX_EVENTS})"

echo ""
echo "All done at $(date)"
