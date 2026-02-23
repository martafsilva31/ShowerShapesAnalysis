#!/bin/bash
###############################################################################
# run_closure_test.sh
#
# Driver script for the full closure test pipeline.
# Runs all steps sequentially and generates final comparison plots.
#
# Usage:
#   ./run_closure_test.sh [N_EVENTS] [DISTORTION_LEVEL]
#
# Arguments:
#   N_EVENTS          Optional. Number of events to process (default: 500000)
#   DISTORTION_LEVEL  Optional. Gaussian noise sigma (default: 0.03)
#                     Common values: 0.01 (1%), 0.03 (3%), 0.05 (5%)
#
# Prerequisites:
#   - ROOT 6.32+ (via lsetup root)
#   - MC23e ntuple at input path specified in config.h
###############################################################################

set -e  # Exit on first error

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# MC23e input ntuple
INPUT_FILE="${SCRIPT_DIR}/../../ntuples/mc23e/mc23e_700770_Zeeg.root"

N_EVENTS=${1:-500000}
DIST_LEVEL=${2:-0.05}
NOISE_SIGMA=${3:--1}

# Compute percentage label for directory names (e.g., 0.03 → 3pct)
PCT_LABEL=$(echo "${DIST_LEVEL}" | awk '{printf "%g", $1*100}')pct

OUTPUT_DIR="${SCRIPT_DIR}/../../output/closure_test_${PCT_LABEL}"
PLOT_DIR="${SCRIPT_DIR}/../../plots/closure_test_${PCT_LABEL}"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${YELLOW}============================================${NC}"
echo -e "${YELLOW}  MC23e Cell Reweighting Closure Test${NC}"
echo -e "${YELLOW}============================================${NC}"
echo ""
echo "Script directory:   ${SCRIPT_DIR}"
echo "Input file:         ${INPUT_FILE}"
echo "Output directory:   ${OUTPUT_DIR}"
echo "Plot directory:     ${PLOT_DIR}"
echo "Events to process:  ${N_EVENTS}"
echo "Distortion level:   ${DIST_LEVEL} (${PCT_LABEL})"
echo "Noise sigma:        ${NOISE_SIGMA} (-1 = auto: level/2)"
echo ""

# Create output directories if needed
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${PLOT_DIR}"

cd "${SCRIPT_DIR}"

# Check ROOT is available
if ! command -v root &> /dev/null; then
    echo -e "${RED}ERROR: ROOT not found. Run: lsetup 'root 6.32.08-x86_64-el9-gcc13-opt'${NC}"
    exit 1
fi

# Check input file exists
if [[ ! -f "${INPUT_FILE}" ]]; then
    echo -e "${RED}ERROR: Input file not found: ${INPUT_FILE}${NC}"
    exit 1
fi

# =============================================================================
# Step 1: Create pseudo-data
# =============================================================================
echo -e "${GREEN}[1/6] Creating pseudo-data...${NC}"
root -l -b -q "create_pseudodata.C(\"${INPUT_FILE}\", \"${OUTPUT_DIR}/pseudo.root\", ${N_EVENTS}, ${DIST_LEVEL}, ${NOISE_SIGMA})"
echo ""

# =============================================================================
# Step 2: Derive corrections
# =============================================================================
echo -e "${GREEN}[2/6] Deriving corrections...${NC}"
root -l -b -q "derive_corrections.C(\"${INPUT_FILE}\", \"${OUTPUT_DIR}/pseudo.root\", \"${OUTPUT_DIR}/corrections.root\", ${N_EVENTS})"
echo ""

# =============================================================================
# Step 3: Apply corrections (both methods)
# =============================================================================
echo -e "${GREEN}[3/6] Applying corrections (Method 1: Flat Shift)...${NC}"
root -l -b -q "apply_corrections.C(\"${INPUT_FILE}\", \"${OUTPUT_DIR}/corrections.root\", \"${OUTPUT_DIR}/reweighted_m1.root\", ${N_EVENTS}, 1)"
echo ""

echo -e "${GREEN}[4/6] Applying corrections (Method 2: TProfile)...${NC}"
root -l -b -q "apply_corrections.C(\"${INPUT_FILE}\", \"${OUTPUT_DIR}/corrections.root\", \"${OUTPUT_DIR}/reweighted_m2.root\", ${N_EVENTS}, 2)"
echo ""

# =============================================================================
# Step 4: Validate closure (compute histograms and chi2)
# =============================================================================
echo -e "${GREEN}[5/6] Validating closure...${NC}"
root -l -b -q "validate_closure.C(\"${INPUT_FILE}\", \"${OUTPUT_DIR}\", \"${OUTPUT_DIR}/closure_histos.root\", ${N_EVENTS})"
echo ""

# =============================================================================
# Step 5: Generate plots
# =============================================================================
echo -e "${GREEN}[6/6] Generating plots...${NC}"
root -l -b -q "plot_closure.C(\"${OUTPUT_DIR}/closure_histos.root\", \"${PLOT_DIR}/\")"
echo ""

# =============================================================================
# Summary
# =============================================================================
echo -e "${YELLOW}============================================${NC}"
echo -e "${GREEN}  Closure Test Complete!${NC}"
echo -e "${YELLOW}============================================${NC}"
echo ""
echo "Output files:"
ls -lh "${OUTPUT_DIR}"/*.root
echo ""
echo "Plot files:"
ls -lh "${PLOT_DIR}"/*.pdf 2>/dev/null || echo "(no PDF files yet)"
echo ""
echo -e "${GREEN}Done.${NC}"
