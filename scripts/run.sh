#!/usr/bin/env bash
# Run an R script from the project root with a timestamped log.
#
# Usage:
#   scripts/run.sh <script.R> [script args...]
#
# Examples:
#   scripts/run.sh scripts/split_corrs.R
#   scripts/run.sh scripts/split_corrs.R --regions IT --out /tmp/it.csv
#
#   # independent jobs in parallel (here, 4 at a time)
#   printf '%s\n' V1 V2 V4 IT | xargs -P 4 -I{} \
#     scripts/run.sh scripts/split_corrs.R --regions {} \
#       --out data/processed/20261104/split_corrs_{}.csv
#
# Logs go to $ADA_LOG_DIR (default: <project>/logs), so the same command works
# locally and on the server (e.g., export ADA_LOG_DIR=/scratch/$USER/logs).

set -euo pipefail

if [[ $# -lt 1 ]]; then
    sed -n '2,17p' "$0"
    exit 1
fi

PROJECT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/.." &> /dev/null && pwd )"
LOG_DIR="${ADA_LOG_DIR:-${PROJECT_DIR}/logs}"
mkdir -p "$LOG_DIR"
cd "$PROJECT_DIR"

SCRIPT="$1"
shift
NAME="$(basename "$SCRIPT" .R)"
TAG="$(printf '%s' "$*" | tr -c 'A-Za-z0-9,._-' '_' | cut -c1-60)"
LOG_FILE="${LOG_DIR}/${NAME}${TAG:+_${TAG}}_$(date +%Y%m%d_%H%M%S).log"

echo "[$(date '+%F %T')] ${SCRIPT} $* (log: ${LOG_FILE})"
if Rscript "$SCRIPT" "$@" > "$LOG_FILE" 2>&1; then
    echo "[$(date '+%F %T')] SUCCESS: ${SCRIPT} $*"
else
    echo "[$(date '+%F %T')] ERROR: ${SCRIPT} $*. Last lines of log:" >&2
    tail -n 20 "$LOG_FILE" >&2
    exit 1
fi
