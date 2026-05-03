#!/bin/bash
# run_pipeline.sh — Run the full derive+validate+plot pipeline
# Usage:
#   ./run_pipeline.sh <channel> [maxEvents]
#
# Channels: Zeeg, Zmumug, combined
# maxEvents: optional, limit for testing (default: all)

CHANNEL="${1:?Usage: $0 <Zeeg|Zmumug|combined> [maxEvents]}"
MAX_EVENTS="${2:--1}"

DATADIR=/dcache/atlas/mfernand/qt_ntuples/data24
OUTDIR="../../output/data_mc_reweighting/${CHANNEL}"
SCRIPTDIR="$(cd "$(dirname "$0")" && pwd)"

mkdir -p "$OUTDIR"

# Determine input files
case "$CHANNEL" in
    Zeeg)
        DATA="${DATADIR}/egam3.root"
        MC="${DATADIR}/mc_eegamma.root"
        LABEL='Z#rightarrowee#gamma'
        ;;
    Zmumug)
        DATA="${DATADIR}/egam4.root"
        MC="${DATADIR}/mc_mumugamma.root"
        LABEL='Z#rightarrow#mu#mu#gamma'
        ;;
    combined)
        DATA="${DATADIR}/egam3.root,${DATADIR}/egam4.root"
        MC="${DATADIR}/mc_eegamma.root,${DATADIR}/mc_mumugamma.root"
        LABEL='Z#rightarrowll#gamma'
        ;;
    *)
        echo "ERROR: Unknown channel '$CHANNEL'. Use Zeeg, Zmumug, or combined."
        exit 1
        ;;
esac

cd "$SCRIPTDIR"

echo "========================================"
echo "  Channel:   $CHANNEL"
echo "  Data:      $DATA"
echo "  MC:        $MC"
echo "  Output:    $OUTDIR"
echo "  MaxEvt:    $MAX_EVENTS"
echo "========================================"

echo ""
echo "=== Step 1: Derive corrections ==="
root -l -b -q "derive_corrections.C(\"${DATA}\",\"${MC}\",\"${OUTDIR}/corrections.root\",${MAX_EVENTS})"

echo ""
echo "=== Step 2: Validate ==="
root -l -b -q "validate_data_mc.C(\"${DATA}\",\"${MC}\",\"${OUTDIR}/corrections.root\",\"${OUTDIR}/histos.root\",${MAX_EVENTS})"

echo ""
echo "=== Step 3: Plot shower shapes ==="
root -l -b -q "plot_data_mc.C(\"${OUTDIR}/histos.root\",\"${OUTDIR}/plots/\")"

echo ""
echo "=== Step 4: Plot cell profiles ==="
root -l -b -q "plot_cell_profiles.C(\"${OUTDIR}/corrections.root\",\"${OUTDIR}/plots/\",\"${LABEL}\")"

echo ""
echo "=== Step 5: Plot stored vs computed ==="
root -l -b -q "plot_stored_vs_computed.C(\"${OUTDIR}/histos.root\",\"${OUTDIR}/plots/\",\"${LABEL}\")"

echo ""
echo "=== Done: $CHANNEL ==="
echo "  Output in: $OUTDIR"
