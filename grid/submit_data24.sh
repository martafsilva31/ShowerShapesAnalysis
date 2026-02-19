#!/bin/bash
set -euo pipefail

# Submit data24 eegamma NTupleMaker jobs to the grid
# Run 473235 verified to be in the data24 GRL:
#   GoodRunsLists/data24_13p6TeV/20241118/physics_25ns_data24.xml
#
# Usage: source this script from any directory (it sets up the full environment)
#   bash ShowerShapesAnalysis/grid/submit_data24.sh

WORKDIR="/project/atlas/users/mfernand/QT/ShowerShapes/NTupleMaker_workspace"

# --- Environment setup ---
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source "${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh" --quiet
asetup Athena,25.0.40
source "${WORKDIR}/build/x86_64-el9-gcc14-opt/setup.sh" 2>/dev/null \
  || source "${WORKDIR}/build/x86_64-el9-gcc13-opt/setup.sh"
lsetup panda rucio

# Ensure valid grid proxy
voms-proxy-info --exists 2>/dev/null || voms-proxy-init -voms atlas

# Must submit from the run/ directory inside the workspace
cd "${WORKDIR}/run"

samples_data24=(
  data24_13p6TeV:data24_13p6TeV.00473235.physics_Main.merge.AOD.r15810_p6304
)

NFILESPERJOB=5
OWNER="femarta"
tag="Zeeg_data24"

# Unique version tag so every resubmission gets a different outDS name
version="v$(date +%Y%m%d_%H%M)"

for sample in "${samples_data24[@]}"; do
  echo "Submitting: ${sample}"
  dsid="$(cut -d'.' -f2 <<<"${sample}")"
  outDS="user.${OWNER}.${dsid}.${tag}.${version}"
  echo "  outDS=${outDS}  nFilesPerJob=${NFILESPERJOB}"

  pathena \
    --trf "python -m NTupleMaker.jobConfig --inputFile=%IN --outputFile=%OUT.root --MaxEvents=-1 -y 2024 -a eegamma -d 0" \
    --inDS "${sample}" \
    --outDS "${outDS}" \
    --nFilesPerJob "${NFILESPERJOB}"
done

echo "Submission script finished"
