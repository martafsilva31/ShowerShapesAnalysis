#!/bin/bash
###############################################################################
# resubmit_477048_egam3.sh
#
# Resubmit the 5 failed files from run 477048 (EGAM3 / eegamma).
# Luca has fixed the corrupted file; this resubmission should succeed.
#
# Usage:
#   bash grid/resubmit_477048_egam3.sh
###############################################################################

set -eo pipefail

WORKDIR="/project/atlas/users/mfernand/QT/ShowerShapes/NTupleMaker_workspace"

# --- Environment setup ---
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source "${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh" --quiet || true
asetup Athena,25.0.40 || true
# Force re-sourcing local build setup
unset UserAnalysis_SET_UP 2>/dev/null
source "${WORKDIR}/build/x86_64-el9-gcc14-opt/setup.sh" 2>/dev/null \
  || source "${WORKDIR}/build/x86_64-el9-gcc13-opt/setup.sh" 2>/dev/null || true
export RUCIO_ACCOUNT=femarta
lsetup panda rucio || true

# Verify local build
if [[ "${CMAKE_PREFIX_PATH:-}" != *"NTupleMaker_workspace/build"* ]]; then
    echo "ERROR: Local build not found in CMAKE_PREFIX_PATH. pathena will fail."
    exit 1
fi

# Ensure valid grid proxy
voms-proxy-info --exists 2>/dev/null || voms-proxy-init -voms atlas

# Must submit from run/ directory
cd "${WORKDIR}/run"

DATASET="data24_13p6TeV.00477048.physics_Main.deriv.DAOD_EGAM3.f1474_m2248_p6143"
VERSION="v$(date +%Y%m%d_%H%M)"
OUTDS="user.femarta.00477048.Zeeg_data24_egam3.${VERSION}"

echo "========================================"
echo " Resubmitting run 477048 (EGAM3)"
echo "  inDS  = ${DATASET}"
echo "  outDS = ${OUTDS}"
echo "  version = ${VERSION}"
echo "========================================"

pathena \
    --trf "python -m NTupleMaker.jobConfig --inputFile=%IN --outputFile=%OUT.root --MaxEvents=-1 -y 2024 -a eegamma -d 0" \
    --inDS "${DATASET}" \
    --outDS "${OUTDS}" \
    --nFilesPerJob 5

echo ""
echo "=== Submission complete. Monitor on BigPandaMon. ==="
echo "=== After completion, re-download and re-merge with grid/download_ntuples.sh ==="
