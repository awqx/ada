#!/usr/bin/env bash

set -euo pipefail

# --- CONFIGURATION ---

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
PROJECT_DIR="$( cd "${SCRIPT_DIR}/.." &> /dev/null && pwd )"
R_SCRIPT="scripts/dM_window.R"
LOG_DIR="${PROJECT_DIR}/logs"
#LOG_DIR="/scratch/axin/logs"
mkdir -p "$LOG_DIR"
cd "$PROJECT_DIR"

# --- PARAMETER GRID ---

REGIONS=("V1" "V2" "V4")
ONSETS=(1 2 3 4 5 6 7 8 9 10)
WINDOWS=(5 10 15)
#ONSETS=(6 7)
#WINDOWS=(5)
# -------------------------

echo "========================================================="
echo "Running: ${R_SCRIPT}"
echo "========================================================="

for REGION in "${REGIONS[@]}"; do
    for ONSET in "${ONSETS[@]}"; do
        for L in ${WINDOWS[@]}; do
            TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
            LOG_FILE="${LOG_DIR}/${REGION}_${ONSET}_${L}_${TIMESTAMP}.log"

            echo "Params:  (region: ${REGION}, onset: ${ONSET}, l: ${L})"
            echo "Time:    $(date '+%Y-%m-%d %H:%M:%S')"
            echo "---------------------------------------------------------"

            if Rscript "$R_SCRIPT" "$REGION" "$ONSET" "$L"> "$LOG_FILE" 2>&1; then
                echo "SUCCESS."
            else
                echo "ERROR. Log:" >&2
                tail -n 20 "$LOG_FILE" >&2

                # comment out below if want loop to continue if error
                exit 1
            fi
            echo "---------------------------------------------------------"
        done
    done
done

echo "========================================================="
echo "Sliding window calculations complete."
