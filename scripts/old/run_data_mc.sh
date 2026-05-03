#!/bin/bash
###############################################################################
# run_data_mc.sh
#
# Driver script for the data-MC cell-energy reweighting pipeline.
# Runs the full 3-step pipeline for each channel:
#   1. derive_corrections.C  — derive M1 (flat shift) + M2 (profiled TProfile)
#   2. validate_data_mc.C    — apply corrections + fill comparison histograms
#   3. plot_data_mc.C        — generate branch_*.pdf (Set A) + rew_*.pdf (Set B)
#
# Usage:
#   ./run_data_mc.sh                    # Both channels, all events
#   ./run_data_mc.sh 500000             # Both channels, 500K events
#   ./run_data_mc.sh -1 eeg             # Z→eeγ only, all events
#   ./run_data_mc.sh 500000 mumug       # Z→μμγ only, 500K events
#
# Prerequisites:
#   - ROOT 6.28+ (via lsetup root)
#   - Merged ntuples at NTUPLE_DIR
###############################################################################

set -e  # Exit on first error

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ======================================================================
# Input files (merged ntuples on dcache)
# ======================================================================
NTUPLE_DIR="/dcache/atlas/mfernand/qt_ntuples/data24"

DATA_EEG="${NTUPLE_DIR}/egam3.root"
MC_EEG="${NTUPLE_DIR}/mc_eegamma.root"

DATA_MUMUG="${NTUPLE_DIR}/egam4.root"
MC_MUMUG="${NTUPLE_DIR}/mc_mumugamma.root"

# ======================================================================
# Arguments
# ======================================================================
N_EVENTS=${1:--1}
CHANNEL=${2:-"both"}   # "eeg", "mumug", or "both"

OUTPUT_BASE="${SCRIPT_DIR}/../../output/data_mc_reweighting"

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

echo -e "${YELLOW}============================================${NC}"
echo -e "${YELLOW}  Data-MC Cell-Energy Reweighting Pipeline${NC}"
echo -e "${YELLOW}============================================${NC}"
echo ""
echo "Script directory:  ${SCRIPT_DIR}"
echo "Ntuple directory:  ${NTUPLE_DIR}"
echo "Output base:       ${OUTPUT_BASE}"
echo "Events per file:   ${N_EVENTS}"
echo "Channel:           ${CHANNEL}"
echo ""

# ======================================================================
# Check ROOT is available
# ======================================================================
if ! command -v root &> /dev/null; then
    echo -e "${RED}ERROR: ROOT not found. Please run:${NC}"
    echo "  source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh --quiet"
    echo "  lsetup \"root 6.28.04-x86_64-el9-gcc13-opt\" --quiet"
    exit 1
fi
echo -e "${GREEN}ROOT found: $(root-config --version)${NC}"

# ======================================================================
# Check input files
# ======================================================================
check_files() {
    for f in "$@"; do
        if [[ ! -f "$f" ]]; then
            echo -e "${RED}ERROR: Input file not found: $f${NC}"
            exit 1
        fi
    done
}

cd "${SCRIPT_DIR}"

# ======================================================================
# Run pipeline for one channel
# ======================================================================
run_channel() {
    local DATA_FILE="$1"
    local MC_FILE="$2"
    local CHANNEL_NAME="$3"
    local CHANNEL_LABEL="$4"
    local OUTPUT_DIR="${OUTPUT_BASE}/${CHANNEL_NAME}"

    echo -e "${YELLOW}--- Channel: ${CHANNEL_NAME} ---${NC}"
    mkdir -p "${OUTPUT_DIR}"

    # Step 0: Diagnostic — cell-computed vs branch comparison
    echo -e "${GREEN}[0/3] Diagnostic: cell-computed vs branch...${NC}"
    root -l -b -q "compare_computed_vs_branch.C(\"${DATA_FILE}\", \"${MC_FILE}\", \"${OUTPUT_DIR}/\", ${N_EVENTS})"
    echo ""

    # Step 1: Derive corrections
    echo -e "${GREEN}[1/3] Deriving corrections...${NC}"
    root -l -b -q "derive_corrections.C(\"${DATA_FILE}\", \"${MC_FILE}\", \"${OUTPUT_DIR}/corrections.root\", ${N_EVENTS})"
    echo ""

    # Step 2: Validate (apply corrections inline + fill histograms)
    echo -e "${GREEN}[2/3] Applying corrections and filling histograms...${NC}"
    root -l -b -q "validate_data_mc.C(\"${DATA_FILE}\", \"${MC_FILE}\", \"${OUTPUT_DIR}/corrections.root\", \"${OUTPUT_DIR}/histos.root\", ${N_EVENTS})"
    echo ""

    # Step 3: Plot
    echo -e "${GREEN}[3/3] Generating plots...${NC}"
    root -l -b -q "plot_data_mc.C(\"${OUTPUT_DIR}/histos.root\", \"${OUTPUT_DIR}/\", \"${CHANNEL_LABEL}\")"
    echo ""

    echo -e "${GREEN}Channel ${CHANNEL_NAME} complete.${NC}"
    echo "  Output: ${OUTPUT_DIR}/"
    ls -lh "${OUTPUT_DIR}/"
    echo ""
}

# ======================================================================
# Execute
# ======================================================================
if [[ "${CHANNEL}" == "eeg" || "${CHANNEL}" == "both" ]]; then
    check_files "${DATA_EEG}" "${MC_EEG}"
    run_channel "${DATA_EEG}" "${MC_EEG}" "Zeeg" "Z#rightarrowee#gamma"
fi

if [[ "${CHANNEL}" == "mumug" || "${CHANNEL}" == "both" ]]; then
    check_files "${DATA_MUMUG}" "${MC_MUMUG}"
    run_channel "${DATA_MUMUG}" "${MC_MUMUG}" "Zmumug" "Z#rightarrow#mu#mu#gamma"
fi

echo -e "${YELLOW}============================================${NC}"
echo -e "${GREEN}  Pipeline Complete!${NC}"
echo -e "${YELLOW}============================================${NC}"
