#!/bin/bash

samples_mc23e=(
	mc23_13p6TeV:mc23_13p6TeV.700770.Sh_2214_eegamma.merge.AOD.e8514_e8586_s4369_s4370_r16083_r15970
	mc23_13p6TeV:mc23_13p6TeV.700771.Sh_2214_mumugamma.merge.AOD.e8514_e8586_s4369_s4370_r16083_r15970
)

NFILESPERJOB=5

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
	version="v1"
	outDS="user.${OWNER}.${dsid}.${tag}.${version}"
	echo "  outDS=${outDS}  nFilesPerJob=${NFILESPERJOB}  analysis=${analysis}"
	pathena --trf "python -m NTupleMaker.jobConfig --inputFile=%IN --outputFile=%OUT.root --MaxEvents=-1 -y 2023 -a ${analysis} -d 0" \
		   --inDS $sample --outDS ${outDS} --nFilesPerJob ${NFILESPERJOB}
done

echo "Submission script finished"
