#!/bin/bash
set -euo pipefail

WORKDIR="/project/atlas/users/mfernand/QT/NTupleMaker_workspace"

INDS="data22_13p6TeV.00428770.physics_Main.merge.AOD.r15869_p6304_tid40703510_00"

OWNER="femarta"
MAXEVENTS=200
NFILES_TOTAL=1
NFILESPERJOB=1

command -v pathena >/dev/null || { echo "ERROR: pathena not found (run setupATLAS; lsetup panda)"; exit 1; }

cd "$WORKDIR"

# Local sanity check: jobConfig must be importable
export PYTHONPATH="$WORKDIR/python:$PYTHONPATH"
python -c "import jobConfig; print('Local jobConfig OK:', jobConfig.__file__)" || {
  echo "ERROR: jobConfig not importable locally. Do not submit."
  exit 1
}

# Short, policy-safe outDS
OUTDS="user.${OWNER}.zeeg22_canary_${RANDOM}"

echo "INDS  = '${INDS}'"
echo "OUTDS = '${OUTDS}'"

# Dummy jobOptions file required by pathena
JO="jobO_dummy.py"
cat > "$JO" <<'EOF'
# Dummy jobOptions for pathena --trf submissions
EOF

# Single pathena call: build tarball AND submit
pathena \
  --outTarBall=workarea.tgz \
  --trf "export PYTHONPATH=\$PWD/python:\$PYTHONPATH; python -m jobConfig --inputFile=%IN --outputFile=%OUT.root --MaxEvents=${MAXEVENTS} -y 2022 -a eegamma -d 0" \
  --inDS "${INDS}" \
  --outDS "${OUTDS}" \
  --nFiles "${NFILES_TOTAL}" \
  --nFilesPerJob "${NFILESPERJOB}" \
  --extOutFile "OUT.root" \
  "$JO"