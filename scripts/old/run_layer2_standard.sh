#!/bin/bash
###############################################################################
# run_layer2_standard.sh — Layer 2, standard 14-bin eta binning
#
# Runs all 9 channel × scenario combinations (3 channels × 3 conversions)
# with peripheral-cell M2 regularisation (kMinFracForM2).
#
# Output: output/Layer_2/standard_eta_binning/{channel}/{scenario}/
#
# Includes: fill, shower-shape plots, cell profiles, chi2, compendium.
# Excludes: kinematics plots (no pT/eta/zero-padding).
#
# Usage:
#   nohup ./run_layer2_standard.sh > run_layer2_standard.log 2>&1 &
#
# Requirements:
#   source /cvmfs/sft.cern.ch/lcg/views/LCG_105a/x86_64-el9-gcc13-opt/setup.sh
###############################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="../../output/Layer_2/standard_eta_binning"
REPORT_DIR="${BASE_DIR}/report"
cd "${SCRIPT_DIR}"

CHANNELS=("eegamma" "mumugamma" "llgamma")
SCENARIOS=("baseline" "converted" "all_conv")

TOTAL=$((${#CHANNELS[@]} * ${#SCENARIOS[@]}))
echo "============================================"
echo "  Layer 2 — Standard eta binning"
echo "  ${TOTAL} scenarios"
echo "  Output: ${BASE_DIR}"
echo "  Started: $(date)"
echo "============================================"

# ── Phase 1: Fill histograms in parallel ──────────────────────────
echo ""
echo ">>> Phase 1: Filling histograms (${TOTAL} parallel jobs)"
echo ""

PIDS=()
LABELS=()
for ch in "${CHANNELS[@]}"; do
    for sc in "${SCENARIOS[@]}"; do
        label="${ch}/${sc}"
        logfile="${BASE_DIR}/${ch}/${sc}/fill.log"
        mkdir -p "${BASE_DIR}/${ch}/${sc}"
        echo "  Starting: ${label}"
        root -l -b -q "fill_histograms.C(\"${ch}\", \"${sc}\", \"${BASE_DIR}\")" \
            > "${logfile}" 2>&1 &
        PIDS+=($!)
        LABELS+=("${label}")
    done
done

echo ""
echo "  Waiting for ${#PIDS[@]} fill jobs..."
FAILED=0
for i in "${!PIDS[@]}"; do
    if wait "${PIDS[$i]}"; then
        echo "  OK:   ${LABELS[$i]}"
    else
        echo "  FAIL: ${LABELS[$i]}"
        ((FAILED++))
    fi
done

if [[ ${FAILED} -gt 0 ]]; then
    echo ""
    echo "  WARNING: ${FAILED} fill job(s) failed. Continuing with plotting..."
fi

# ── Phase 2: Plot (no kinematics) ────────────────────────────────
echo ""
echo ">>> Phase 2: Plotting (shower shapes + cell profiles)"
echo ""

for ch in "${CHANNELS[@]}"; do
    for sc in "${SCENARIOS[@]}"; do
        echo "  Plotting: ${ch}/${sc}"
        root -l -b -q "plot_shower_shapes.C(\"${ch}\", \"${sc}\", \"${BASE_DIR}\")" || \
            echo "  WARN: plot_shower_shapes failed for ${ch}/${sc}"
        root -l -b -q "plot_cell_profiles.C(\"${ch}\", \"${sc}\", \"${BASE_DIR}\")" || \
            echo "  WARN: plot_cell_profiles failed for ${ch}/${sc}"
    done
done

# ── Phase 3: Extract chi2 tables ─────────────────────────────────
echo ""
echo ">>> Phase 3: Extracting chi-squared tables"
echo ""
mkdir -p "${REPORT_DIR}"
root -l -b -q "extract_chi2.C(\"${BASE_DIR}\", \"${REPORT_DIR}\")" || \
    echo "  WARN: extract_chi2 failed"

# ── Phase 4: Generate compendia ──────────────────────────────────
echo ""
echo ">>> Phase 4: Generating compendia (no kinematics)"
echo ""
python3 make_compendium.py --base-dir "${BASE_DIR}" --no-kinematics || \
    echo "  WARN: make_compendium failed"

# ── Done ─────────────────────────────────────────────────────────
echo ""
echo "============================================"
echo "  Layer 2 standard run complete"
echo "  Finished: $(date)"
echo "============================================"
echo ""
echo "  Output:     ${BASE_DIR}/{channel}/{scenario}/"
echo "  Chi2:       ${REPORT_DIR}/chi2_*.tex"
echo "  Compendia:  ${BASE_DIR}/{channel}/{scenario}/result_compendium_*.pdf"
