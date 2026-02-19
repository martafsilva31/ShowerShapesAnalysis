#!/bin/bash
set -euo pipefail

# Always submit from the workspace root (so pathena ships the right sandbox)
WORKDIR="/project/atlas/users/mfernand/QT/NTupleMaker_workspace"
cd "${WORKDIR}"

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