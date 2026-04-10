#!/bin/bash

# Submit MC23e DAOD NTupleMaker jobs to the grid
# Reads dataset list from NTupleMaker/script/mc23_13TeV.DAOD_EGAMx.25ns.p5xxx.txt
# Auto-detects analysis type (eegamma/mumugamma) from sample name
# Year: 2023, DAOD cell-links (-d 0)
#
# Usage:
#   bash grid/submit_mc23e_DAOD.sh          # submit all samples
#   bash grid/submit_mc23e_DAOD.sh --test   # submit first sample only (1 file)

WORKDIR="/project/atlas/users/mfernand/QT/ShowerShapes/NTupleMaker_workspace"
SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DSLIST="${WORKDIR}/source/NTupleMaker/script/mc23_13TeV.DAOD_EGAMx.25ns.p5xxx.txt"

# --- Parse arguments ---
TEST_MODE=false
if [[ "${1:-}" == "--test" ]]; then
    TEST_MODE=true
    echo "=== TEST MODE: submitting first sample only (--nFiles 1) ==="
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

# --- Read datasets from file (skip comment, blank, and header lines) ---
mapfile -t DATASETS < <(grep -v '^\s*#\|^\s*$' "${DSLIST}")

NFILESPERJOB=2
OWNER="femarta"

if $TEST_MODE; then
    VERSION="vtest4"
else
    VERSION="v1"
fi

SUBMITTED=0
for dataset in "${DATASETS[@]}"; do
    # Extract DSID (field 2: e.g. 700770)
    dsid="$(cut -d'.' -f2 <<<"${dataset}")"

    # Auto-detect analysis type from sample name
    if [[ "${dataset}" == *"eegamma"* ]]; then
        ANALYSIS="eegamma"
        TAG="Zeeg_mc23e_DAOD"
    elif [[ "${dataset}" == *"mumugamma"* ]]; then
        ANALYSIS="mumugamma"
        TAG="Zmumug_mc23e_DAOD"
    else
        echo "WARNING: Cannot determine analysis type for ${dataset}, skipping"
        continue
    fi

    outDS="user.${OWNER}.${dsid}.${TAG}.${VERSION}"

    echo "Submitting: ${dataset}"
    echo "  outDS=${outDS}  nFilesPerJob=${NFILESPERJOB}  analysis=${ANALYSIS}"

    EXTRA_ARGS=""
    if $TEST_MODE; then
        EXTRA_ARGS="--nFiles 1"
    fi

    pathena \
        --trf "python -m NTupleMaker.jobConfig --inputFile=%IN --outputFile=%OUT.root --MaxEvents=-1 -y 2023 -a ${ANALYSIS} -d 0" \
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
