#!/usr/bin/env bash

set -euo pipefail

# --- CONFIGURATION ---

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
PROJECT_DIR="$( cd "${SCRIPT_DIR}/.." &> /dev/null && pwd )"
R_SCRIPT="scripts/split_corrs.R"
LOG_DIR="${PROJECT_DIR}/logs"
#LOG_DIR="/scratch/axin/logs"
mkdir -p "$LOG_DIR"
cd "$PROJECT_DIR"

# --- PARAMETER GRID ---

# use command line arguments if passed
# else, use the sample list

REGIONS=("V1" "V2" "V4" "IT")
# -------------------------

for REGION in "${REGIONS[@]}"; do

    # Create a unique log file name using both arguments and timestamp
    TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
    LOG_FILE="${LOG_DIR}/split_corrs_${REGION}_${TIMESTAMP}.log"

    echo "========================================================="
    echo "Running: ${R_SCRIPT}"
    echo "Params:  [region: ${REGION}]"
    echo "Time:    $(date '+%Y-%m-%d %H:%M:%S')"
    echo "Log:     ${LOG_FILE}"
    echo "========================================================="

    # Execute Rscript, passing both arguments in order
    if Rscript "$R_SCRIPT" "$REGION" > "$LOG_FILE" 2>&1; then
        echo "SUCCESS: split-half correlations calculated for ${REGION}."
    else
        echo "ERROR: Region ${REGION} failed. Log:" >&2
        tail -n 20 "$LOG_FILE" >&2

        # comment to let loop continue on error
        exit 1
    fi

done

echo "========================================================="
echo "Split-half correlation calculations complete."
