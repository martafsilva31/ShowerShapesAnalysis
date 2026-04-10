#!/bin/bash
###############################################################################
# run.sh — Driver for the cell-energy reweighting pipeline
#
# Usage:
#   ./run.sh                          # defaults: eegamma, baseline
#   ./run.sh eegamma baseline
#   ./run.sh mumugamma tight_id
#   ./run.sh eegamma baseline --plot-only   # skip fill, only plot
#   ./run.sh --batch                        # all 3 channels × 3 scenarios
#   ./run.sh --batch --plot-only            # batch plots only
#
# Requirements:
#   - ATLAS environment (asetup Athena,25.0.40)
#   - ROOT available on PATH
###############################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="../../output/cell_energy_reweighting_Francisco_method/data24"
PLOT_ONLY=false
BATCH=false

# Parse arguments
for arg in "$@"; do
    case "$arg" in
        --batch)     BATCH=true ;;
        --plot-only) PLOT_ONLY=true ;;
    esac
done

# ── Batch mode: loop over all channels × scenarios ──
if [[ "${BATCH}" == true ]]; then
    CHANNELS=("eegamma" "mumugamma" "llgamma")
    SCENARIOS=("baseline" "converted" "all_conv")
    cd "${SCRIPT_DIR}"
    for ch in "${CHANNELS[@]}"; do
        for sc in "${SCENARIOS[@]}"; do
            echo ""
            echo "========================================"
            echo " Running: ${ch} / ${sc}"
            echo "========================================"
            if [[ "${PLOT_ONLY}" == false ]]; then
                echo ">>> fill_histograms.C"
                root -l -b -q "fill_histograms.C(\"${ch}\", \"${sc}\", \"${BASE_DIR}\")"
            fi
            echo ">>> plot_shower_shapes.C"
            root -l -b -q "plot_shower_shapes.C(\"${ch}\", \"${sc}\", \"${BASE_DIR}\")"
            root -l -b -q "plot_shower_shapes.C(\"${ch}\", \"${sc}\", \"${BASE_DIR}\", 2)"
            echo ">>> plot_cell_profiles.C"
            root -l -b -q "plot_cell_profiles.C(\"${ch}\", \"${sc}\", \"${BASE_DIR}\")"
        done
    done
    echo ""
    echo "========================================"
    echo " Batch complete: ${#CHANNELS[@]} channels × ${#SCENARIOS[@]} scenarios"
    echo "========================================"
    exit 0
fi

# ── Single-run mode ──
CHANNEL="${1:-eegamma}"
SCENARIO="${2:-baseline}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="../../output/cell_energy_reweighting_Francisco_method/data24"

echo "========================================"
echo " Cell-energy reweighting pipeline"
echo "  channel  = ${CHANNEL}"
echo "  scenario = ${SCENARIO}"
echo "  base     = ${BASE_DIR}"
echo "========================================"

cd "${SCRIPT_DIR}"

# Step 1: Fill histograms (two-pass: statistics + corrections)
if [[ "${PLOT_ONLY}" == false ]]; then
    echo ""
    echo ">>> Step 1: fill_histograms.C"
    root -l -b -q "fill_histograms.C(\"${CHANNEL}\", \"${SCENARIO}\", \"${BASE_DIR}\")"
fi

# Step 2: Shower shape plots
echo ""
echo ">>> Step 2: plot_shower_shapes.C"
root -l -b -q "plot_shower_shapes.C(\"${CHANNEL}\", \"${SCENARIO}\", \"${BASE_DIR}\")"
root -l -b -q "plot_shower_shapes.C(\"${CHANNEL}\", \"${SCENARIO}\", \"${BASE_DIR}\", 2)"

# Step 3: Cell profile plots
echo ""
echo ">>> Step 3: plot_cell_profiles.C"
root -l -b -q "plot_cell_profiles.C(\"${CHANNEL}\", \"${SCENARIO}\", \"${BASE_DIR}\")"

OUT_DIR="${BASE_DIR}/${CHANNEL}/${SCENARIO}"
echo ""
echo "========================================"
echo " Done."
echo " Histograms: ${OUT_DIR}/histograms.root"
echo " Plots:      ${OUT_DIR}/plots/"
echo "========================================"
