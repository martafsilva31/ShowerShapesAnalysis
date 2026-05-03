#!/bin/bash
###############################################################################
# run_kinematics.sh — Photon pT and |eta| kinematic distributions
#
# Produces normalised and unnormalised kinematics plots for llgamma
# across 3 conversion types × 2 isolation selections (loose, iso_tight).
#
# Output structure:
#   output/kinematics/loose/llgamma/{baseline,converted,all_conv}/plots/
#   output/kinematics/iso_tight/llgamma/{baseline,converted,all_conv}/plots/
#
# Each scenario directory contains:
#   kinematics_pt.pdf       — unnormalised pT  (y = "Events")
#   kinematics_eta.pdf      — unnormalised |η| (y = "Events")
#   kinematics_pt_norm.pdf  — normalised pT    (y = "Normalised")
#   kinematics_eta_norm.pdf — normalised |η|   (y = "Normalised")
#   zero_padding.pdf        — zero-padding diagnostic
#
# After plots, builds a kinematics compendium PDF.
#
# Usage:
#   nohup ./run_kinematics.sh > run_kinematics.log 2>&1 &
#
# Requirements:
#   source /cvmfs/sft.cern.ch/lcg/views/LCG_105a/x86_64-el9-gcc13-opt/setup.sh
###############################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

CHANNEL="llgamma"
CONVERSIONS=("baseline" "converted" "all_conv")

# ── Scenario table ────────────────────────────────────────────────
#   baseDir              scenario(dir)  selectionOverride
# Loose isolation (default — no override needed)
#   kinematics/loose     baseline       (empty)
#   kinematics/loose     converted      (empty)
#   kinematics/loose     all_conv       (empty)
# Tight isolation (override selection to iso_tight variants)
#   kinematics/iso_tight baseline       iso_tight
#   kinematics/iso_tight converted      iso_tight_converted
#   kinematics/iso_tight all_conv       iso_tight_all_conv

echo "============================================"
echo "  Kinematics: llgamma × 3 conv × 2 iso"
echo "  Started: $(date)"
echo "============================================"

# ── Phase 1: Fill histograms in parallel ──────────────────────────
echo ""
echo ">>> Phase 1: Filling histograms (6 parallel jobs)"
echo ""

PIDS=()
LABELS=()

fill_job() {
    local base_dir="$1" scenario="$2" sel_override="$3" label="$4"
    local outdir="${base_dir}/${CHANNEL}/${scenario}"
    mkdir -p "${outdir}"
    local logfile="${outdir}/fill.log"
    echo "  Starting: ${label}"
    if [[ -z "${sel_override}" ]]; then
        root -l -b -q "fill_histograms.C(\"${CHANNEL}\", \"${scenario}\", \"${base_dir}\")" \
            > "${logfile}" 2>&1 &
    else
        root -l -b -q "fill_histograms.C(\"${CHANNEL}\", \"${scenario}\", \"${base_dir}\", \"${sel_override}\")" \
            > "${logfile}" 2>&1 &
    fi
    PIDS+=($!)
    LABELS+=("${label}")
}

# Loose isolation
for conv in "${CONVERSIONS[@]}"; do
    fill_job "../../output/kinematics/loose" "${conv}" "" "loose/${conv}"
done

# Tight isolation
for conv in "${CONVERSIONS[@]}"; do
    sel="iso_tight"
    [[ "${conv}" == "converted" ]] && sel="iso_tight_converted"
    [[ "${conv}" == "all_conv"  ]] && sel="iso_tight_all_conv"
    fill_job "../../output/kinematics/iso_tight" "${conv}" "${sel}" "iso_tight/${conv}"
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

# ── Phase 2: Plot kinematics (normalised + unnormalised) ──────────
echo ""
echo ">>> Phase 2: Plotting kinematics"
echo ""

plot_kin() {
    local base_dir="$1" scenario="$2" sel_override="$3" label="$4"

    echo "  Plotting: ${label} (unnormalised)"
    if [[ -z "${sel_override}" ]]; then
        root -l -b -q "plot_kinematics.C(\"${CHANNEL}\", \"${scenario}\", \"${base_dir}\", 0)" || \
            echo "  WARN: plot_kinematics (unnorm) failed for ${label}"
    else
        root -l -b -q "plot_kinematics.C(\"${CHANNEL}\", \"${scenario}\", \"${base_dir}\", 0, \"${sel_override}\")" || \
            echo "  WARN: plot_kinematics (unnorm) failed for ${label}"
    fi

    echo "  Plotting: ${label} (normalised)"
    if [[ -z "${sel_override}" ]]; then
        root -l -b -q "plot_kinematics.C(\"${CHANNEL}\", \"${scenario}\", \"${base_dir}\", 1)" || \
            echo "  WARN: plot_kinematics (norm) failed for ${label}"
    else
        root -l -b -q "plot_kinematics.C(\"${CHANNEL}\", \"${scenario}\", \"${base_dir}\", 1, \"${sel_override}\")" || \
            echo "  WARN: plot_kinematics (norm) failed for ${label}"
    fi
}

# Loose
for conv in "${CONVERSIONS[@]}"; do
    plot_kin "../../output/kinematics/loose" "${conv}" "" "loose/${conv}"
done

# Tight
for conv in "${CONVERSIONS[@]}"; do
    sel="iso_tight"
    [[ "${conv}" == "converted" ]] && sel="iso_tight_converted"
    [[ "${conv}" == "all_conv"  ]] && sel="iso_tight_all_conv"
    plot_kin "../../output/kinematics/iso_tight" "${conv}" "${sel}" "iso_tight/${conv}"
done

# ── Phase 3: Build kinematics compendium ─────────────────────────
echo ""
echo ">>> Phase 3: Building kinematics compendium"
echo ""
python3 make_kinematics_compendium.py || \
    echo "  WARN: kinematics compendium failed"

# ── Done ─────────────────────────────────────────────────────────
echo ""
echo "============================================"
echo "  Kinematics run complete"
echo "  Finished: $(date)"
echo "============================================"
echo ""
echo "  Output: output/kinematics/{loose,iso_tight}/llgamma/{baseline,converted,all_conv}/plots/"
echo "  Compendium: output/kinematics/kinematics_compendium.pdf"
