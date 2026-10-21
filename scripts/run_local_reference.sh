#!/usr/bin/env bash
# SPDX-License-Identifier: MIT
# Repeatable local reference experiment for the cartesian vector EFG prototype.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${BUILD_DIR:-${ROOT}/build}"
PROFILE="${PROFILE:-local}"  # local: 9x7,13x11,17x15; stargate: 9x7,13x11
JOBS="${JOBS:-2}"
CASE_TIMEOUT="${CASE_TIMEOUT:-900}"
CTEST_BIN="${CTEST_BIN:-ctest}"
OUT_DIR="${OUT_DIR:-${ROOT}/out/m4/${PROFILE}}"
if ! [[ "${JOBS}" =~ ^[1-9][0-9]*$ ]]; then echo "JOBS must be positive" >&2; exit 2; fi
if ! [[ "${CASE_TIMEOUT}" =~ ^[1-9][0-9]*$ ]]; then echo "CASE_TIMEOUT must be positive" >&2; exit 2; fi
case "${PROFILE}" in local|stargate|smoke|ci) ;; *) echo "Invalid PROFILE" >&2; exit 2;; esac
echo "=== Al Maleh M4: local vector EFG reference experiment ==="
echo "profile=${PROFILE} build=${BUILD_DIR} output=${OUT_DIR} jobs=${JOBS}"
if command -v free >/dev/null; then free -h; fi
cmake -S "${ROOT}" -B "${BUILD_DIR}" -DCMAKE_BUILD_TYPE=Release
cmake --build "${BUILD_DIR}" --parallel "${JOBS}"
"${CTEST_BIN}" --test-dir "${BUILD_DIR}" --output-on-failure
python3 "${ROOT}/scripts/m4_campaign.py" --binary "${BUILD_DIR}/tm_validation" \
    --out "${OUT_DIR}" --profile "${PROFILE}" --timeout "${CASE_TIMEOUT}"
echo "=== M4 files ==="
echo "${OUT_DIR}/REPORT.md"
echo "${OUT_DIR}/tm_matches_all.csv"
echo "${OUT_DIR}/runs.csv"
