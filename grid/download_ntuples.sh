#!/bin/bash
###############################################################################
# download_ntuples.sh
#
# Downloads DAOD_EGAM3 ntuples from completed grid jobs, then merges them
# into single ROOT files for the data-MC comparison pipeline.
#
# Grid tasks:
#   Data24: 49030538  (data24 run 473235 DAOD_EGAM3, 1 file test)
#   MC23e:  49032024  (mc23e 700770 Zeeg DAOD_EGAM3)
#
# Output:
#   ntuples/data24_egam3/data24_egam3.root   (merged data)
#   ntuples/mc23e_egam3/mc23e_egam3.root     (merged MC)
#
# Prerequisites:
#   - lsetup rucio, RUCIO_ACCOUNT=femarta
#   - Valid grid proxy: voms-proxy-init -voms atlas
#
# Usage:
#   ./download_ntuples.sh           # download + merge
#   ./download_ntuples.sh --skip-download   # merge only (if already downloaded)
###############################################################################

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
NTUPLE_DIR="${BASE_DIR}/ntuples"

# Output directories
DATA_DIR="${NTUPLE_DIR}/data24_egam3"
MC_DIR="${NTUPLE_DIR}/mc23e_egam3"
DATA_RAW="${DATA_DIR}/raw"
MC_RAW="${MC_DIR}/raw"

# Grid task IDs
DATA_TASK=49030538
MC_TASK=49032024

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

SKIP_DOWNLOAD=false
if [[ "$1" == "--skip-download" ]]; then
    SKIP_DOWNLOAD=true
fi

echo -e "${YELLOW}============================================${NC}"
echo -e "${YELLOW}  Download DAOD_EGAM3 Ntuples${NC}"
echo -e "${YELLOW}============================================${NC}"
echo ""

# Create directories
mkdir -p "${DATA_RAW}" "${MC_RAW}"

# =============================================================================
# Step 1: Find output dataset names from task IDs
# =============================================================================
if [[ "${SKIP_DOWNLOAD}" == false ]]; then

    echo -e "${GREEN}[1/4] Looking up output datasets...${NC}"

    # Query rucio for the output dataset containers
    DATA_DS=$(rucio list-dids "user.femarta.*" --filter "type=CONTAINER" 2>/dev/null \
        | grep -i "egam3" | grep -i "data24" | tail -1 | awk '{print $2}') || true

    MC_DS=$(rucio list-dids "user.femarta.*" --filter "type=CONTAINER" 2>/dev/null \
        | grep -i "egam3" | grep -i "mc23e\|700770" | tail -1 | awk '{print $2}') || true

    # If automatic lookup fails, try bigpanda-style naming convention
    if [[ -z "${DATA_DS}" ]]; then
        echo -e "${YELLOW}  Auto-lookup for data DS failed. Trying manual lookup...${NC}"
        echo "  Please find the output DS name from: https://bigpanda.cern.ch/task/${DATA_TASK}/"
        echo "  Then set DATA_DS below and re-run, or use:"
        echo "    rucio list-dids 'user.femarta.00473235.*egam3*' --short"
        DATA_DS=$(rucio list-dids "user.femarta.00473235.*egam3*" --short 2>/dev/null \
            | grep "_output" | head -1) || true
    fi

    if [[ -z "${MC_DS}" ]]; then
        echo -e "${YELLOW}  Auto-lookup for MC DS failed. Trying manual lookup...${NC}"
        MC_DS=$(rucio list-dids "user.femarta.700770.*egam3*" --short 2>/dev/null \
            | grep "_output" | head -1) || true
    fi

    echo "  Data DS: ${DATA_DS:-NOT FOUND}"
    echo "  MC DS:   ${MC_DS:-NOT FOUND}"
    echo ""

    # =============================================================================
    # Step 2: Download with rucio
    # =============================================================================
    echo -e "${GREEN}[2/4] Downloading ntuples...${NC}"

    if [[ -n "${DATA_DS}" ]]; then
        echo "  Downloading data to ${DATA_RAW}/ ..."
        rucio download "${DATA_DS}" --dir "${DATA_RAW}" --ndownloader 4
    else
        echo -e "${RED}  ERROR: Data dataset not found. Set DATA_DS manually.${NC}"
    fi
    echo ""

    if [[ -n "${MC_DS}" ]]; then
        echo "  Downloading MC to ${MC_RAW}/ ..."
        rucio download "${MC_DS}" --dir "${MC_RAW}" --ndownloader 4
    else
        echo -e "${RED}  ERROR: MC dataset not found. Set MC_DS manually.${NC}"
    fi
    echo ""

else
    echo -e "${YELLOW}  Skipping download (--skip-download)${NC}"
    echo ""
fi

# =============================================================================
# Step 3: Merge data files
# =============================================================================
echo -e "${GREEN}[3/4] Merging data ntuples...${NC}"

DATA_FILES=$(find "${DATA_RAW}" -name "*.root" -type f 2>/dev/null | sort)
N_DATA=$(echo "${DATA_FILES}" | grep -c "\.root$" 2>/dev/null || echo 0)

if [[ ${N_DATA} -gt 0 ]]; then
    echo "  Found ${N_DATA} data file(s)"
    if [[ ${N_DATA} -eq 1 ]]; then
        cp ${DATA_FILES} "${DATA_DIR}/data24_egam3.root"
    else
        hadd -f "${DATA_DIR}/data24_egam3.root" ${DATA_FILES}
    fi
    echo "  Merged: ${DATA_DIR}/data24_egam3.root"
    echo "  Size: $(du -h "${DATA_DIR}/data24_egam3.root" | cut -f1)"
    root -l -b -q -e "TFile f(\"${DATA_DIR}/data24_egam3.root\"); auto t=(TTree*)f.Get(\"tree\"); std::cout << \"  Entries: \" << t->GetEntries() << std::endl;"
else
    echo -e "${RED}  No data ROOT files found in ${DATA_RAW}${NC}"
fi
echo ""

# =============================================================================
# Step 4: Merge MC files
# =============================================================================
echo -e "${GREEN}[4/4] Merging MC ntuples...${NC}"

MC_FILES=$(find "${MC_RAW}" -name "*.root" -type f 2>/dev/null | sort)
N_MC=$(echo "${MC_FILES}" | grep -c "\.root$" 2>/dev/null || echo 0)

if [[ ${N_MC} -gt 0 ]]; then
    echo "  Found ${N_MC} MC file(s)"
    if [[ ${N_MC} -eq 1 ]]; then
        cp ${MC_FILES} "${MC_DIR}/mc23e_egam3.root"
    else
        hadd -f "${MC_DIR}/mc23e_egam3.root" ${MC_FILES}
    fi
    echo "  Merged: ${MC_DIR}/mc23e_egam3.root"
    echo "  Size: $(du -h "${MC_DIR}/mc23e_egam3.root" | cut -f1)"
    root -l -b -q -e "TFile f(\"${MC_DIR}/mc23e_egam3.root\"); auto t=(TTree*)f.Get(\"tree\"); std::cout << \"  Entries: \" << t->GetEntries() << std::endl;"
else
    echo -e "${RED}  No MC ROOT files found in ${MC_RAW}${NC}"
fi

echo ""
echo -e "${YELLOW}============================================${NC}"
echo -e "${GREEN}  Download & Merge Complete${NC}"
echo -e "${YELLOW}============================================${NC}"
echo ""
echo "Output:"
echo "  Data: ${DATA_DIR}/data24_egam3.root"
echo "  MC:   ${MC_DIR}/mc23e_egam3.root"
echo ""
echo "Next: run the data-MC comparison pipeline:"
echo "  cd ../scripts/data_mc && ./run_data_mc.sh"
