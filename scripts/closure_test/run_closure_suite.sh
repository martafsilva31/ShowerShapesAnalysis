#!/bin/bash
###############################################################################
# run_closure_suite.sh
#
# Runs the full closure test pipeline at multiple distortion levels.
# This script calls run_closure_test.sh for each level.
#
# Usage:
#   ./run_closure_suite.sh [N_EVENTS]
#
# Arguments:
#   N_EVENTS  Optional. Number of events per level (default: 500000)
#
# Distortion levels: 1%, 3%, 5%
###############################################################################

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
N_EVENTS=${1:-500000}

LEVELS=(0.01 0.03 0.05)

# Colors for output
CYAN='\033[0;36m'
YELLOW='\033[1;33m'
GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m'

echo -e "${CYAN}=====================================================${NC}"
echo -e "${CYAN}  Closure Test Suite — Multiple Distortion Levels${NC}"
echo -e "${CYAN}=====================================================${NC}"
echo ""
echo "Events per level: ${N_EVENTS}"
echo "Distortion levels: ${LEVELS[*]}"
echo ""

FAILED=()
for LEVEL in "${LEVELS[@]}"; do
    PCT=$(echo "${LEVEL}" | awk '{printf "%g", $1*100}')
    echo -e "${YELLOW}>>> Starting ${PCT}% distortion level (sigma=${LEVEL}) <<<${NC}"
    echo ""

    if "${SCRIPT_DIR}/run_closure_test.sh" "${N_EVENTS}" "${LEVEL}"; then
        echo -e "${GREEN}>>> ${PCT}% distortion level completed successfully <<<${NC}"
    else
        echo -e "${RED}>>> ${PCT}% distortion level FAILED <<<${NC}"
        FAILED+=("${LEVEL}")
    fi
    echo ""
done

echo -e "${CYAN}=====================================================${NC}"
echo -e "${CYAN}  Suite Summary${NC}"
echo -e "${CYAN}=====================================================${NC}"

if [[ ${#FAILED[@]} -eq 0 ]]; then
    echo -e "${GREEN}All ${#LEVELS[@]} distortion levels completed successfully.${NC}"
else
    echo -e "${RED}Failed levels: ${FAILED[*]}${NC}"
fi

echo ""
echo "Output directories:"
for LEVEL in "${LEVELS[@]}"; do
    PCT=$(echo "${LEVEL}" | awk '{printf "%g", $1*100}')pct
    echo "  output/closure_test_${PCT}/"
    echo "  plots/closure_test_${PCT}/"
done
