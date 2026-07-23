#!/usr/bin/env bash

set -euo pipefail

# --- CONFIGURATION ---

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
PROJECT_DIR="$( cd "${SCRIPT_DIR}/.." &> /dev/null && pwd )"
R_SCRIPT="scripts/pf_lda.R"
LOG_DIR="${PROJECT_DIR}/logs"

mkdir -p "$LOG_DIR"
cd "$PROJECT_DIR"

# --- ARGUMENT HANDLING ---

# use command line arguments if passed
# else, use the sample list

if [ "$#" -gt 0 ]; then
    TARGET_ARGS=("$@")
else
    TARGET_ARGS=(10 20 30 40)
fi

# -------------------------

echo "Pairwise comparisons with sample size(s): ${TARGET_ARGS[*]}"

# --- LOOP THROUGH EACH ARGUMENT ---
for ARG in "${TARGET_ARGS[@]}"; do

    # Create a unique log file per argument and timestamp
    TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
    LOG_FILE="${LOG_DIR}/pf_lda_${ARG}_${TIMESTAMP}.log"

    echo "========================================================="
    echo "Running: ${R_SCRIPT} with parameter: ${ARG}"
    echo "Time:    $(date '+%Y-%m-%d %H:%M:%S')"
    echo "Log:     ${LOG_FILE}"
    echo "========================================================="

    # Execute Rscript passing the current argument ($ARG)
    if Rscript "$R_SCRIPT" "$ARG" > "$LOG_FILE" 2>&1; then
        echo "SUCCESS: Finished evaluation with ${ARG} subsamples"
    else
        echo "ERROR: Subsamples of size ${ARG} failed. Log:" >&2
        tail -n 20 "$LOG_FILE" >&2

        # Uncomment to stop on error
        # exit 1
    fi

done

echo "========================================================="
echo "COMPLETE"
