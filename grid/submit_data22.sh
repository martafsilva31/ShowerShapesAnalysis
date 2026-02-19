#!/bin/bash
set -euo pipefail

# Submit data22 eegamma NTupleMaker jobs to the grid
#
# NOTE: The runs below (428xxx) are NOT in the data22 GRL (which starts at 431810).
#       These will produce 0 events. Update run numbers before resubmitting.

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

samples_data22=(
  data22_13p6TeV:data22_13p6TeV.00428648.physics_Main.merge.AOD.r15869_p6304_tid40703502_00
  data22_13p6TeV:data22_13p6TeV.00428700.physics_Main.merge.AOD.r15869_p6304_tid40703504_00
  data22_13p6TeV:data22_13p6TeV.00428747.physics_Main.merge.AOD.r15869_p6304_tid40703506_00
  data22_13p6TeV:data22_13p6TeV.00428759.physics_Main.merge.AOD.r15869_p6304_tid40703508_00
  data22_13p6TeV:data22_13p6TeV.00428770.physics_Main.merge.AOD.r15869_p6304_tid40703510_00
  data22_13p6TeV:data22_13p6TeV.00428776.physics_Main.merge.AOD.r15869_p6304_tid40703512_00
  data22_13p6TeV:data22_13p6TeV.00428777.physics_Main.merge.AOD.r15869_p6304_tid40703514_00
  data22_13p6TeV:data22_13p6TeV.00428855.physics_Main.merge.AOD.r15869_p6304_tid40703516_00
  data22_13p6TeV:data22_13p6TeV.00429018.physics_Main.merge.AOD.r15869_p6304_tid40703518_00
  data22_13p6TeV:data22_13p6TeV.00429027.physics_Main.merge.AOD.r15869_p6304_tid40703520_00
)

NFILESPERJOB=5
OWNER="femarta"
tag="Zeeg_data22"

# Unique version tag so every resubmission gets a different outDS name
# Example: v20260126_2115
version="v$(date +%Y%m%d_%H%M)"


for sample in "${samples_data22[@]}"; do
  echo "Submitting: ${sample}"
  dsid="$(cut -d'.' -f2 <<<"${sample}")"
  outDS="user.${OWNER}.${dsid}.${tag}.${version}"
  echo "  outDS=${outDS}  nFilesPerJob=${NFILESPERJOB}"

  pathena \
    --trf "python -m NTupleMaker.jobConfig --inputFile=%IN --outputFile=%OUT.root --MaxEvents=-1 -y 2022 -a eegamma -d 0" \
    --inDS "${sample}" \
    --outDS "${outDS}" \
    --nFilesPerJob "${NFILESPERJOB}"
done

echo "Submission script finished"