#!/bin/bash
###############################################################################
# run_layer2_histsigma.sh — Layer 2, histogram-based sigma test
#
# Runs llgamma × 3 scenarios (unconverted, converted, inclusive) only.
# Uses the new histogram-based sigma for M2 corrections (matching Francisco).
#
# Output: output/Layer_2/hist_sigma_test/{channel}/{scenario}/
#
# Usage:
#   nohup bash run_layer2_histsigma.sh > run_layer2_histsigma.log 2>&1 &
###############################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="../../output/Layer_2/hist_sigma_test"
REPORT_DIR="${BASE_DIR}/report"
cd "${SCRIPT_DIR}"

CHANNELS=("llgamma")
SCENARIOS=("unconverted" "converted" "inclusive")

TOTAL=$((${#CHANNELS[@]} * ${#SCENARIOS[@]}))
echo "============================================"
echo "  Layer 2 — Histogram-based sigma test"
echo "  ${TOTAL} scenarios (llgamma only)"
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
        outroot="${BASE_DIR}/${ch}/${sc}/histograms.root"
        logfile="${BASE_DIR}/${ch}/${sc}/fill.log"
        mkdir -p "${BASE_DIR}/${ch}/${sc}"
        if [[ -f "${outroot}" && $(stat -c%s "${outroot}") -gt 100000 ]]; then
            echo "  Skipping (exists): ${label}"
            # fake a successful PID that exits immediately
            (exit 0) &
        else
            echo "  Starting: ${label}"
            root -l -b -q "fill_histograms.C(\"${ch}\", \"${sc}\", \"${BASE_DIR}\")" \
                > "${logfile}" 2>&1 &
        fi
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
python3 make_compendium.py --base-dir "${BASE_DIR}" --no-kinematics --channels llgamma || \
    echo "  WARN: make_compendium failed"

# ── Done ─────────────────────────────────────────────────────────
echo ""
echo "============================================"
echo "  Histogram-sigma test complete"
echo "  Finished: $(date)"
echo "============================================"
echo ""
echo "  Output:     ${BASE_DIR}/{channel}/{scenario}/"
echo "  Chi2:       ${REPORT_DIR}/chi2_*.tex"
echo "  Compendia:  ${BASE_DIR}/{channel}/{scenario}/result_compendium_*.pdf"
