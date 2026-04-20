#!/usr/bin/env bash
#
# run_layer2_final.sh
#
# Runs the full Layer 2 cell-energy reweighting pipeline for 4 variants:
#   eta_loose    — 14 eta bins, loose isolation
#   eta_tight    — 14 eta bins, tight isolation
#   eta_pt_loose — 14 eta x 6 pT bins, loose isolation
#   eta_pt_tight — 14 eta x 6 pT bins, tight isolation
#
# Each variant: fill (unc+conv parallel) → hadd inclusive → plot → chi2 → compendium
#
# Usage:
#   cd scripts/data_mc
#   source /cvmfs/sft.cern.ch/lcg/views/LCG_105a/x86_64-el9-gcc13-opt/setup.sh
#   bash run_layer2_final.sh
#
# Set VARIANTS environment variable to run a subset:
#   VARIANTS="eta_loose" bash run_layer2_final.sh
#
set -euo pipefail

# ── Environment setup ─────────────────────────────────────────────
LCG_SETUP="/cvmfs/sft.cern.ch/lcg/views/LCG_105a/x86_64-el9-gcc13-opt/setup.sh"
if [[ -f "${LCG_SETUP}" ]]; then
    set +u  # LCG setup uses unbound variables
    # shellcheck disable=SC1090
    source "${LCG_SETUP}"
    set -u
else
    echo "ERROR: LCG setup not found: ${LCG_SETUP}" >&2
    exit 1
fi

echo "============================================"
echo "  Layer 2: Final cell-energy reweighting"
echo "  Started: $(date)"
echo "============================================"

# ── Configuration ─────────────────────────────────────────────────
SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTBASE="${SCRIPTDIR}/../../output/Layer_2"

CHANNELS=(llgamma)
SCENARIOS_FILL=(unconverted converted)
SCENARIOS_ALL=(unconverted converted inclusive)

# Variant definitions: name → (binning, isolation)
declare -A BINNING
declare -A ISOLATION
BINNING[eta_loose]="eta";      ISOLATION[eta_loose]="loose"
BINNING[eta_tight]="eta";      ISOLATION[eta_tight]="tight"
BINNING[eta_pt_loose]="eta_pt"; ISOLATION[eta_pt_loose]="loose"
BINNING[eta_pt_tight]="eta_pt"; ISOLATION[eta_pt_tight]="tight"

ALL_VARIANTS=(eta_loose eta_tight eta_pt_loose eta_pt_tight)
VARIANTS_RUN=("${VARIANTS:-${ALL_VARIANTS[@]}}")

# Parse VARIANTS from environment if set as a single string
if [[ "${#VARIANTS_RUN[@]}" -eq 1 && "${VARIANTS_RUN[0]}" == *" "* ]]; then
    read -ra VARIANTS_RUN <<< "${VARIANTS_RUN[0]}"
fi

echo ""
echo "  Variants: ${VARIANTS_RUN[*]}"
echo "  Channels: ${CHANNELS[*]}"
echo "  Output:   ${OUTBASE}/<variant>/"
echo ""

cd "${SCRIPTDIR}"

# ── Phase 1: Fill histograms (all variants × scenarios in parallel) ──
echo ">>> Phase 1: Filling histograms (all variants in parallel)"
echo ""

FILL_PIDS=()
for var in "${VARIANTS_RUN[@]}"; do
    BIN="${BINNING[$var]}"
    ISO="${ISOLATION[$var]}"
    BASE="${OUTBASE}/${var}"

    for ch in "${CHANNELS[@]}"; do
        for sc in "${SCENARIOS_FILL[@]}"; do
            outf="${BASE}/${ch}/${sc}/histograms.root"
            if [[ -f "${outf}" ]]; then
                echo "  SKIP (exists): ${var}/${ch}/${sc}"
                continue
            fi
            echo "  Filling: ${var}/${ch}/${sc}  [binning=${BIN}, iso=${ISO}]"
            mkdir -p "${BASE}/${ch}/${sc}"
            root -l -b -q "fill_histograms.C(\"${ch}\", \"${sc}\", \"${BASE}\", \"${BIN}\", \"${ISO}\")" \
                > "${BASE}/${ch}/${sc}/fill.log" 2>&1 &
            FILL_PIDS+=($!)
        done
    done
done

echo "  Waiting for ${#FILL_PIDS[@]} fill jobs..."
for pid in "${FILL_PIDS[@]}"; do
    wait "${pid}" || echo "  WARN: fill job ${pid} failed"
done
echo "  All fill jobs done."

# ── Phase 2: hadd inclusive = unconverted + converted ─────────────
echo ""
echo ">>> Phase 2: Creating inclusive histograms (hadd)"
echo ""

for var in "${VARIANTS_RUN[@]}"; do
    BASE="${OUTBASE}/${var}"
    for ch in "${CHANNELS[@]}"; do
        incdir="${BASE}/${ch}/inclusive"
        mkdir -p "${incdir}"
        uncfile="${BASE}/${ch}/unconverted/histograms.root"
        convfile="${BASE}/${ch}/converted/histograms.root"
        if [[ -f "${uncfile}" && -f "${convfile}" ]]; then
            echo "  hadd: ${var}/${ch}/inclusive"
            hadd -f "${incdir}/histograms.root" "${uncfile}" "${convfile}" \
                > "${incdir}/hadd.log" 2>&1 || \
                echo "  WARN: hadd failed for ${var}/${ch}"
        else
            echo "  SKIP: missing input for ${var}/${ch}/inclusive"
        fi
    done
done

# ── Phase 3: Plots ───────────────────────────────────────────────
echo ""
echo ">>> Phase 3: Plotting (shower shapes + cell profiles)"
echo ""

for var in "${VARIANTS_RUN[@]}"; do
    BIN="${BINNING[$var]}"
    ISO="${ISOLATION[$var]}"
    BASE="${OUTBASE}/${var}"

    for ch in "${CHANNELS[@]}"; do
        for sc in "${SCENARIOS_ALL[@]}"; do
            echo "  Plotting: ${var}/${ch}/${sc}"
            root -l -b -q "plot_shower_shapes.C(\"${ch}\", \"${sc}\", \"${BASE}\", \"${BIN}\", \"${ISO}\")" || \
                echo "  WARN: plot_shower_shapes failed for ${var}/${ch}/${sc}"
            root -l -b -q "plot_cell_profiles.C(\"${ch}\", \"${sc}\", \"${BASE}\", \"${BIN}\", \"${ISO}\")" || \
                echo "  WARN: plot_cell_profiles failed for ${var}/${ch}/${sc}"
        done
    done
done

# ── Phase 4: Extract chi2 tables ─────────────────────────────────
echo ""
echo ">>> Phase 4: Extracting chi-squared tables"
echo ""

for var in "${VARIANTS_RUN[@]}"; do
    BASE="${OUTBASE}/${var}"
    REPORT="${BASE}/report"
    mkdir -p "${REPORT}"
    echo "  chi2: ${var}"
    root -l -b -q "extract_chi2.C(\"${BASE}\", \"${REPORT}\")" || \
        echo "  WARN: extract_chi2 failed for ${var}"
done

# ── Phase 5: Generate compendia ──────────────────────────────────
echo ""
echo ">>> Phase 5: Generating compendia"
echo ""

# Check if updated make_compendiums_final.py exists; fall back to make_compendium.py
COMP_SCRIPT="${OUTBASE}/make_compendiums_final.py"
if [[ ! -f "${COMP_SCRIPT}" ]]; then
    COMP_SCRIPT="${SCRIPTDIR}/make_compendium.py"
fi

for var in "${VARIANTS_RUN[@]}"; do
    BASE="${OUTBASE}/${var}"
    echo "  compendium: ${var}"
    if [[ -f "${OUTBASE}/make_compendiums_final.py" ]]; then
        python3 "${OUTBASE}/make_compendiums_final.py" --variant "${var}" || \
            echo "  WARN: make_compendiums_final failed for ${var}"
    else
        python3 "${SCRIPTDIR}/make_compendium.py" --base-dir "${BASE}" --no-kinematics --channels llgamma || \
            echo "  WARN: make_compendium failed for ${var}"
    fi
done

# ── Done ─────────────────────────────────────────────────────────
echo ""
echo "============================================"
echo "  Pipeline complete"
echo "  Finished: $(date)"
echo "============================================"
echo ""
for var in "${VARIANTS_RUN[@]}"; do
    echo "  ${var}: ${OUTBASE}/${var}/{channel}/{scenario}/"
done
echo ""
echo "  Chi2:       <variant>/report/chi2_*.tex"
echo "  Compendia:  <variant>/{channel}/{scenario}/result_compendium_*.pdf"
