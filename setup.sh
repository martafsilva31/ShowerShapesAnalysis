#!/bin/bash
# Setup environment for shower shape analysis
# Source this from the ShowerShapes/ directory:
#   cd /project/atlas/users/mfernand/QT/ShowerShapes
#   source analysis/setup.sh

export SHOWERSHAPES_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export ANALYSIS_DIR="${SHOWERSHAPES_DIR}/analysis"
export NTUPLE_WORKSPACE="${SHOWERSHAPES_DIR}/NTupleMaker_workspace"

echo "ShowerShapes directory: ${SHOWERSHAPES_DIR}"
echo "Analysis directory:     ${ANALYSIS_DIR}"
echo "NTupleMaker workspace:  ${NTUPLE_WORKSPACE}"

# Setup ATLAS environment
setupATLAS 2>/dev/null || echo "Warning: setupATLAS not available (not on lxplus/cluster?)"
asetup Athena,25.0.40 2>/dev/null || echo "Warning: asetup failed"

# Source NTupleMaker build if available
if [[ -f "${NTUPLE_WORKSPACE}/build/x86_64-el9-gcc13-opt/setup.sh" ]]; then
    source "${NTUPLE_WORKSPACE}/build/x86_64-el9-gcc13-opt/setup.sh"
    echo "NTupleMaker build sourced"
elif [[ -f "${NTUPLE_WORKSPACE}/build/x86_64-el9-gcc14-opt/setup.sh" ]]; then
    source "${NTUPLE_WORKSPACE}/build/x86_64-el9-gcc14-opt/setup.sh"
    echo "NTupleMaker build sourced"
else
    echo "Warning: NTupleMaker build not found. Run cmake/make in ${NTUPLE_WORKSPACE}/build/ first."
fi

# Convenience aliases
alias cdanalysis="cd ${ANALYSIS_DIR}"
alias cdntuples="cd ${ANALYSIS_DIR}/ntuples"
alias cdplots="cd ${ANALYSIS_DIR}/plots"
alias cdntuple="cd ${NTUPLE_WORKSPACE}"

echo ""
echo "Environment ready. Shortcuts: cdanalysis, cdntuples, cdplots, cdntuple"
