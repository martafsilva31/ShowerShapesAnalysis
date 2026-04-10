#!/bin/bash
set -euo pipefail

# Submit mc23e Zllg NTupleMaker jobs to the grid

WORKDIR="/project/atlas/users/mfernand/QT/ShowerShapes/NTupleMaker_workspace"

# --- Environment setup ---
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source "${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh" --quiet
asetup Athena,25.0.40
# Force re-sourcing local build setup (guard variable may be inherited from caller)
unset UserAnalysis_SET_UP 2>/dev/null
source "${WORKDIR}/build/x86_64-el9-gcc14-opt/setup.sh" 2>/dev/null \
  || source "${WORKDIR}/build/x86_64-el9-gcc13-opt/setup.sh"
lsetup panda rucio

# Ensure valid grid proxy
voms-proxy-info --exists 2>/dev/null || voms-proxy-init -voms atlas

# Must submit from the run/ directory inside the workspace
cd "${WORKDIR}/run"

samples_mc23e=(
	mc23_13p6TeV:mc23_13p6TeV.700770.Sh_2214_eegamma.merge.AOD.e8514_e8586_s4369_s4370_r16083_r15970
	mc23_13p6TeV:mc23_13p6TeV.700771.Sh_2214_mumugamma.merge.AOD.e8514_e8586_s4369_s4370_r16083_r15970
)

NFILESPERJOB=1

OWNER="femarta"

for sample in "${samples_mc23e[@]}"
do
	echo "Submitting: $sample"
	dsid="$(cut -d'.' -f2 <<<${sample})"
	# Determine analysis type based on sample name
	if [[ "$sample" == *"eegamma"* ]]; then
		analysis="eegamma"
		tag="Zeeg_mc23e"
	elif [[ "$sample" == *"mumugamma"* ]]; then
		analysis="mumugamma"
		tag="Zmumug_mc23e"
	fi
	version="v2"
	outDS="user.${OWNER}.${dsid}.${tag}.${version}"
	echo "  outDS=${outDS}  nFilesPerJob=${NFILESPERJOB}  analysis=${analysis}"
	pathena --trf "python -m NTupleMaker.jobConfig --inputFile=%IN --outputFile=%OUT.root --MaxEvents=-1 -y 2023 -a ${analysis} -d 0" \
		   --inDS $sample --outDS ${outDS} --nFilesPerJob ${NFILESPERJOB}
done

echo "Submission script finished"
