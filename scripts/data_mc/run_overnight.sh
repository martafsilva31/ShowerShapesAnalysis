#!/bin/bash
###############################################################################
# run_overnight.sh — Parallel driver for all 9 scenarios
#
# Fills all 9 channel × scenario combinations in parallel, then runs
# plotting, chi2 extraction, and compendium generation.
#
# Usage:
#   nohup ./run_overnight.sh > overnight_run.log 2>&1 &
#
# Requirements:
#   source /cvmfs/sft.cern.ch/lcg/views/LCG_105a/x86_64-el9-gcc13-opt/setup.sh
###############################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="../../output/cell_energy_reweighting_Francisco_method/data24"
cd "${SCRIPT_DIR}"

CHANNELS=("eegamma" "mumugamma" "llgamma")
SCENARIOS=("baseline" "converted" "all_conv")

TOTAL=$((${#CHANNELS[@]} * ${#SCENARIOS[@]}))
echo "============================================"
echo "  Overnight run: ${TOTAL} scenarios"
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

# ── Phase 2: Plot all scenarios ───────────────────────────────────
echo ""
echo ">>> Phase 2: Plotting"
echo ""

for ch in "${CHANNELS[@]}"; do
    for sc in "${SCENARIOS[@]}"; do
        echo "  Plotting: ${ch}/${sc}"
        root -l -b -q "plot_shower_shapes.C(\"${ch}\", \"${sc}\", \"${BASE_DIR}\")" || \
            echo "  WARN: plot_shower_shapes failed for ${ch}/${sc}"
        root -l -b -q "plot_kinematics.C(\"${ch}\", \"${sc}\", \"${BASE_DIR}\")" || \
            echo "  WARN: plot_kinematics failed for ${ch}/${sc}"
        root -l -b -q "plot_cell_profiles.C(\"${ch}\", \"${sc}\", \"${BASE_DIR}\")" || \
            echo "  WARN: plot_cell_profiles failed for ${ch}/${sc}"
    done
done

# ── Phase 3: Extract chi2 tables ─────────────────────────────────
echo ""
echo ">>> Phase 3: Extracting chi-squared tables"
echo ""
root -l -b -q "extract_chi2.C(\"${BASE_DIR}\", \"../../report\")" || \
    echo "  WARN: extract_chi2 failed"

# ── Phase 4: Generate compendia ──────────────────────────────────
echo ""
echo ">>> Phase 4: Generating compendia"
echo ""
python3 make_compendium.py || \
    echo "  WARN: make_compendium failed"

# ── Done ─────────────────────────────────────────────────────────
echo ""
echo "============================================"
echo "  Overnight run complete"
echo "  Finished: $(date)"
echo "============================================"
echo ""
echo "  Output: ${BASE_DIR}/{channel}/{scenario}/"
echo "  Chi2:   ../../report/chi2_*.tex"
echo "  Compendia: ${BASE_DIR}/{channel}/{scenario}/result_compendium.pdf"
