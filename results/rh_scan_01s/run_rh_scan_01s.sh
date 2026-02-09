#!/usr/bin/env bash
set -euo pipefail

# Parametric RH scan at fixed Te and total simulation time.
# Saves each run output to a dedicated file to avoid overwrite/resume issues.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
SOLVER="${REPO_ROOT}/src/airGM2.1"
OUTPUT_DIR="${SCRIPT_DIR}/outputs"

TE_EV="4.8"
TOTAL_TIME_S="0.1"
RH_VALUES=(0 10 20 30 40 50 60 70 80 90 100)

mkdir -p "${OUTPUT_DIR}"

if [[ ! -x "${SOLVER}" ]]; then
  echo "Solver not found at ${SOLVER}. Building in src/..."
  (cd "${REPO_ROOT}/src" && make)
fi

for RH in "${RH_VALUES[@]}"; do
  RH_PADDED="$(printf "%03d" "${RH}")"
  OUT_FILE="${OUTPUT_DIR}/run_RH${RH_PADDED}.csv"

  echo "Running RH=${RH}% -> ${OUT_FILE}"

  # Prevent solver from resuming from previous output.csv.
  rm -f "${REPO_ROOT}/output.csv"

  "${SOLVER}" -Te "${TE_EV}" -totaltime "${TOTAL_TIME_S}" -RH "${RH}"

  if [[ ! -f "${REPO_ROOT}/output.csv" ]]; then
    echo "ERROR: expected ${REPO_ROOT}/output.csv was not produced."
    exit 1
  fi

  mv "${REPO_ROOT}/output.csv" "${OUT_FILE}"
done

echo "Completed RH scan. Outputs saved in ${OUTPUT_DIR}"
