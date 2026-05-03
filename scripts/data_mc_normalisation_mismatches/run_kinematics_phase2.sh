#!/bin/bash
###############################################################################
# run_kinematics_phase2.sh — Kinematics Phase 2 (plots) + Phase 3 (compendium)
# Assumes Phase 1 (fill) histograms.root files already exist.
###############################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

CHANNEL="llgamma"
CONVERSIONS=("baseline" "converted" "all_conv")

echo "============================================"
echo "  Kinematics Phase 2+3"
echo "  Started: $(date)"
echo "============================================"

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

echo ""
echo "============================================"
echo "  Phase 2+3 complete"
echo "  Finished: $(date)"
echo "============================================"
echo ""
echo "  Compendium: output/kinematics/kinematics_compendium.pdf"
echo "  Plots:"
find ../../output/kinematics -name "*.pdf" | sort
