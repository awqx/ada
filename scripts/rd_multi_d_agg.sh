#!/usr/bin/env bash

set -euo pipefail

# --- CONFIGURATION ---

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
PROJECT_DIR="$( cd "${SCRIPT_DIR}/.." &> /dev/null && pwd )"
R_SCRIPT="scripts/rd_multi_d_agg.R"
LOG_DIR="${PROJECT_DIR}/logs"

mkdir -p "$LOG_DIR"
cd "$PROJECT_DIR"

# -------------------------

# Create a unique log file per argument and timestamp
TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
LOG_FILE="${LOG_DIR}/rd_agg_${TIMESTAMP}.log"

echo "========================================================="
echo "Running: ${R_SCRIPT}"
echo "Time:    $(date '+%Y-%m-%d %H:%M:%S')"
echo "Log:     ${LOG_FILE}"
echo "========================================================="

# Execute Rscript passing the current argument ($ARG)
if Rscript "$R_SCRIPT" > "$LOG_FILE" 2>&1; then
    echo "SUCCESS."
else
    echo "ERROR. Log:" >&2
    tail -n 20 "$LOG_FILE" >&2
fi

