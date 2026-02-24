#!/bin/bash

# =============================================================================
# Submit NTupleMaker grid jobs for ALL data24 samples (267 runs)
# Reads sample list from NTupleMaker/script/data.txt
#
# Usage (from any directory):
#   bash /project/atlas/users/mfernand/QT/ShowerShapes/ShowerShapesAnalysis/grid/submit_data24_all.sh
#
# The script sets up the full environment, so no prior setup is needed.
# =============================================================================

WORKDIR="/project/atlas/users/mfernand/QT/ShowerShapes/NTupleMaker_workspace"
SAMPLE_FILE="${WORKDIR}/source/NTupleMaker/script/data.txt"

# --- Environment setup (ATLAS scripts may return non-zero, so no set -e here) ---
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source "${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh" --quiet || true
asetup Athena,25.0.40 || true
source "${WORKDIR}/build/x86_64-el9-gcc14-opt/setup.sh" 2>/dev/null \
  || source "${WORKDIR}/build/x86_64-el9-gcc13-opt/setup.sh" 2>/dev/null || true
export RUCIO_ACCOUNT=femarta
lsetup panda rucio || true

# Ensure valid grid proxy
voms-proxy-info --exists 2>/dev/null || voms-proxy-init -voms atlas

# Must submit from run/ directory
cd "${WORKDIR}/run"

# From here on, fail on errors
set -eo pipefail

# --- Configuration ---
NFILESPERJOB=5
OWNER="femarta"
TAG="Zeeg_data24"
VERSION="v$(date +%Y%m%d_%H%M)"

# --- Read data24 samples from data.txt ---
mapfile -t samples < <(grep '^data24' "${SAMPLE_FILE}")

echo "============================================="
echo "Data24 submission: ${#samples[@]} samples"
echo "Tag: ${TAG}   Version: ${VERSION}"
echo "nFilesPerJob: ${NFILESPERJOB}"
echo "============================================="
echo ""

submitted=0
failed=0

for sample in "${samples[@]}"; do
  dsid="$(cut -d'.' -f2 <<<"${sample}")"
  outDS="user.${OWNER}.${dsid}.${TAG}.${VERSION}"

  echo "[${submitted}/${#samples[@]}] Submitting: ${sample}"
  echo "  outDS=${outDS}"

  if pathena \
    --trf "python -m NTupleMaker.jobConfig --inputFile=%IN --outputFile=%OUT.root --MaxEvents=-1 --year 2024 --algo eegamma -d 0" \
    --inDS "${sample}" \
    --outDS "${outDS}" \
    --nFilesPerJob "${NFILESPERJOB}"; then
    submitted=$((submitted + 1))
  else
    echo "  WARNING: Submission failed for ${sample}"
    failed=$((failed + 1))
  fi
  echo ""
done

echo "============================================="
echo "Done. Submitted: ${submitted}  Failed: ${failed}  Total: ${#samples[@]}"
echo "============================================="
