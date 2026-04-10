#!/bin/bash

# Submit data24 DAOD_EGAM4 NTupleMaker jobs to the grid
# Reads dataset list from NTupleMaker/script/data24_13p6TeV.DAOD_EGAM4.25ns.p61xx.txt
# Analysis: mumugamma (Z->mumugamma), year 2024, DAOD cell-links (-d 0)
#
# Usage:
#   bash grid/submit_data24_egam4.sh          # submit all 157 runs
#   bash grid/submit_data24_egam4.sh --test   # submit first run only (1 file)

WORKDIR="/project/atlas/users/mfernand/QT/ShowerShapes/NTupleMaker_workspace"
SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DSLIST="${WORKDIR}/source/NTupleMaker/script/data24_13p6TeV.DAOD_EGAM4.25ns.p61xx.txt"

# --- Parse arguments ---
TEST_MODE=false
if [[ "${1:-}" == "--test" ]]; then
    TEST_MODE=true
    echo "=== TEST MODE: submitting first dataset only (--nFiles 1) ==="
fi

# --- Environment setup ---
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source "${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh" --quiet || true
asetup Athena,25.0.40 || true
# Force re-sourcing local build setup (guard variable may be inherited from caller)
unset UserAnalysis_SET_UP 2>/dev/null
source "${WORKDIR}/build/x86_64-el9-gcc14-opt/setup.sh" 2>/dev/null \
  || source "${WORKDIR}/build/x86_64-el9-gcc13-opt/setup.sh" 2>/dev/null || true
export RUCIO_ACCOUNT=femarta
lsetup panda rucio || true

# Verify local build is in CMAKE_PREFIX_PATH (pathena uses this to find CPackConfig.cmake)
if [[ "${CMAKE_PREFIX_PATH:-}" != *"NTupleMaker_workspace/build"* ]]; then
    echo "ERROR: Local build not found in CMAKE_PREFIX_PATH. pathena will fail."
    echo "  CMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH:-unset}"
    exit 1
fi

# Ensure valid grid proxy
voms-proxy-info --exists 2>/dev/null || voms-proxy-init -voms atlas

# Must submit from the run/ directory inside the workspace
cd "${WORKDIR}/run"

# From here on, fail on errors
set -eo pipefail

# --- Read datasets from file (skip comment and blank lines) ---
mapfile -t DATASETS < <(grep -v '^\s*#\|^\s*$' "${DSLIST}")

NFILESPERJOB=5
OWNER="femarta"
TAG="Zmumug_data24_egam4"

if $TEST_MODE; then
    VERSION="vtest4"
else
    VERSION="v1"
fi

SUBMITTED=0
for dataset in "${DATASETS[@]}"; do
    # Extract run number (field 2: 00XXXXXX)
    runid="$(cut -d'.' -f2 <<<"${dataset}")"
    outDS="user.${OWNER}.${runid}.${TAG}.${VERSION}"

    echo "Submitting: ${dataset}"
    echo "  outDS=${outDS}  nFilesPerJob=${NFILESPERJOB}"

    EXTRA_ARGS=""
    if $TEST_MODE; then
        EXTRA_ARGS="--nFiles 1"
    fi

    pathena \
        --trf "python -m NTupleMaker.jobConfig --inputFile=%IN --outputFile=%OUT.root --MaxEvents=-1 -y 2024 -a mumugamma -d 0" \
        --inDS "${dataset}" \
        --outDS "${outDS}" \
        --nFilesPerJob "${NFILESPERJOB}" \
        ${EXTRA_ARGS}

    SUBMITTED=$((SUBMITTED + 1))

    if $TEST_MODE; then
        echo "=== Test submission done (1 job). Monitor on BigPandaMon. ==="
        break
    fi
done

echo "Submission script finished: ${SUBMITTED} jobs submitted"
