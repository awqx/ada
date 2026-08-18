#!/usr/bin/env bash

set -euo pipefail

# --- CONFIGURATION ---

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
PROJECT_DIR="$( cd "${SCRIPT_DIR}/.." &> /dev/null && pwd )"
R_SCRIPT="scripts/rd_multi_d.R"
LOG_DIR="${PROJECT_DIR}/logs"

mkdir -p "$LOG_DIR"
cd "$PROJECT_DIR"

# --- PARAMETER GRID ---

# use command line arguments if passed
# else, use the sample list

N_SITES=(10 20 30 40)
FAM_IDS=(13  18  19  23  30  31  38  43  48  52  56  60  68  71  77  95  99 113 145 196 250 251 291 313 327 330 335 336 375 393 402 438 445 493 542)
# FAM_IDS=(13 18 19)

# -------------------------

echo "Starting randomized dichotomies within families..."
echo "Site numbers: ${N_SITES[*]}"

for N in "${N_SITES[@]}"; do
    for FAM_ID in "${FAM_IDS[@]}"; do

        # Create a unique log file name using both arguments and timestamp
        TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
        LOG_FILE="${LOG_DIR}/rd_n${N}_fam${FAM_ID}_${TIMESTAMP}.log"

        echo "========================================================="
        echo "Running: ${R_SCRIPT}"
        echo "Params:  [Sites: ${N}]  [Family: ${FAM_ID}]"
        echo "Time:    $(date '+%Y-%m-%d %H:%M:%S')"
        echo "Log:     ${LOG_FILE}"
        echo "========================================================="

        # Execute Rscript, passing both arguments in order
        if Rscript "$R_SCRIPT" "$N" "$FAM_ID" > "$LOG_FILE" 2>&1; then
            echo "SUCCESS: Randomized dichotomies for (${N}, ${FAM_ID})"
        else
            echo "ERROR: Combination (${N}, ${FAM_ID}) failed! Check log:" >&2
            tail -n 20 "$LOG_FILE" >&2

            # comment to let loop continue on error
            exit 1
        fi

    done
done

echo "========================================================="
echo "Within-family randomized dichotomies complete."
