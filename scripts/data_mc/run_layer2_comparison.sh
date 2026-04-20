#!/bin/bash
###############################################################################
# run_layer2_comparison.sh — Compare wide vs cell-dependent histogram ranges
#
# Runs llgamma × 2 scenarios (unconverted, converted) × 2 range modes.
# Produces inclusive by combining unconverted + converted via hadd.
#
# Output:
#   output/Layer_2/hist_sigma_wide/{channel}/{scenario}/
#   output/Layer_2/hist_sigma_cellrange/{channel}/{scenario}/
#
# Usage:
#   nohup bash run_layer2_comparison.sh > run_layer2_comparison.log 2>&1 &
###############################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

CHANNELS=("llgamma")
SCENARIOS=("unconverted" "converted")
RANGE_MODES=("wide" "cellrange")

BASE_WIDE="../../output/Layer_2/hist_sigma_wide"
BASE_CELL="../../output/Layer_2/hist_sigma_cellrange"

echo "============================================"
echo "  Layer 2 — Wide vs Cell-dependent ranges"
echo "  ${#CHANNELS[@]} channel × ${#SCENARIOS[@]} scenarios × ${#RANGE_MODES[@]} modes"
echo "  Started: $(date)"
echo "============================================"

# ── Phase 1: Fill histograms (4 parallel jobs) ───────────────────
echo ""
echo ">>> Phase 1: Filling histograms (4 parallel jobs)"
echo ""

# Clean compiled objects to ensure fresh compilation
rm -f fill_histograms_C.so fill_histograms_C.d fill_histograms_C_ACLiC_dict_rdict.pcm

PIDS=()
LABELS=()
for rm in "${RANGE_MODES[@]}"; do
    if [[ "${rm}" == "wide" ]]; then
        BASE="${BASE_WIDE}"
    else
        BASE="${BASE_CELL}"
    fi
    for ch in "${CHANNELS[@]}"; do
        for sc in "${SCENARIOS[@]}"; do
            label="${rm}/${ch}/${sc}"
            outroot="${BASE}/${ch}/${sc}/histograms.root"
            logfile="${BASE}/${ch}/${sc}/fill.log"
            mkdir -p "${BASE}/${ch}/${sc}"
            if [[ -f "${outroot}" && $(stat -c%s "${outroot}" 2>/dev/null || echo 0) -gt 100000 ]]; then
                echo "  Skipping (exists): ${label}"
                (exit 0) &
            else
                echo "  Starting: ${label}"
                root -l -b -q "fill_histograms.C(\"${ch}\", \"${sc}\", \"${BASE}\", \"${rm}\")" \
                    > "${logfile}" 2>&1 &
            fi
            PIDS+=($!)
            LABELS+=("${label}")
        done
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
    echo "  WARNING: ${FAILED} fill job(s) failed. Check logs in output dirs."
fi

# ── Phase 2: Combine inclusive via hadd ──────────────────────────
echo ""
echo ">>> Phase 2: Combining unconverted + converted → inclusive"
echo ""

for rm in "${RANGE_MODES[@]}"; do
    if [[ "${rm}" == "wide" ]]; then
        BASE="${BASE_WIDE}"
    else
        BASE="${BASE_CELL}"
    fi
    for ch in "${CHANNELS[@]}"; do
        incdir="${BASE}/${ch}/inclusive"
        mkdir -p "${incdir}"
        uncfile="${BASE}/${ch}/unconverted/histograms.root"
        convfile="${BASE}/${ch}/converted/histograms.root"
        if [[ -f "${uncfile}" && -f "${convfile}" ]]; then
            echo "  hadd: ${rm}/${ch}/inclusive"
            hadd -f "${incdir}/histograms.root" "${uncfile}" "${convfile}" \
                > "${incdir}/hadd.log" 2>&1 || \
                echo "  WARN: hadd failed for ${rm}/${ch}"
        else
            echo "  SKIP: missing input for ${rm}/${ch}/inclusive"
        fi
    done
done

# ── Phase 3: Plots ───────────────────────────────────────────────
echo ""
echo ">>> Phase 3: Plotting (shower shapes + cell profiles)"
echo ""

ALL_SCENARIOS=("unconverted" "converted" "inclusive")
for rm in "${RANGE_MODES[@]}"; do
    if [[ "${rm}" == "wide" ]]; then
        BASE="${BASE_WIDE}"
    else
        BASE="${BASE_CELL}"
    fi
    for ch in "${CHANNELS[@]}"; do
        for sc in "${ALL_SCENARIOS[@]}"; do
            echo "  Plotting: ${rm}/${ch}/${sc}"
            root -l -b -q "plot_shower_shapes.C(\"${ch}\", \"${sc}\", \"${BASE}\")" || \
                echo "  WARN: plot_shower_shapes failed for ${rm}/${ch}/${sc}"
            root -l -b -q "plot_cell_profiles.C(\"${ch}\", \"${sc}\", \"${BASE}\")" || \
                echo "  WARN: plot_cell_profiles failed for ${rm}/${ch}/${sc}"
        done
    done
done

# ── Phase 4: Extract chi2 tables ─────────────────────────────────
echo ""
echo ">>> Phase 4: Extracting chi-squared tables"
echo ""

for rm in "${RANGE_MODES[@]}"; do
    if [[ "${rm}" == "wide" ]]; then
        BASE="${BASE_WIDE}"
    else
        BASE="${BASE_CELL}"
    fi
    REPORT="${BASE}/report"
    mkdir -p "${REPORT}"
    echo "  chi2: ${rm}"
    root -l -b -q "extract_chi2.C(\"${BASE}\", \"${REPORT}\")" || \
        echo "  WARN: extract_chi2 failed for ${rm}"
done

# ── Phase 5: Generate compendia ──────────────────────────────────
echo ""
echo ">>> Phase 5: Generating compendia"
echo ""

for rm in "${RANGE_MODES[@]}"; do
    if [[ "${rm}" == "wide" ]]; then
        BASE="${BASE_WIDE}"
    else
        BASE="${BASE_CELL}"
    fi
    echo "  compendium: ${rm}"
    python3 make_compendium.py --base-dir "${BASE}" --no-kinematics --channels llgamma || \
        echo "  WARN: make_compendium failed for ${rm}"
done

# ── Done ─────────────────────────────────────────────────────────
echo ""
echo "============================================"
echo "  Comparison complete"
echo "  Finished: $(date)"
echo "============================================"
echo ""
echo "  Wide output:      ${BASE_WIDE}/{channel}/{scenario}/"
echo "  Cellrange output: ${BASE_CELL}/{channel}/{scenario}/"
echo "  Chi2:             {BASE}/report/chi2_*.tex"
echo "  Compendia:        {BASE}/{channel}/{scenario}/result_compendium_*.pdf"
