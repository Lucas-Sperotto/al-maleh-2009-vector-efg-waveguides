#!/usr/bin/env bash
# Rebuild both independent solvers, retain every TM/EFG eigenpair and compare fields.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
FEM_REPO="${FEM_REPO:-$(dirname "${ROOT}")/fem-eigen-em-2d3d-whitney}"
BUILD_DIR="${BUILD_DIR:-${ROOT}/build}"
JOBS="${JOBS:-2}"
export PYTHONNOUSERSITE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
export CTEST_BIN="${CTEST_BIN:-$(dirname "$(command -v cmake)")/ctest}"
if ! [[ "${JOBS}" =~ ^[1-9][0-9]*$ ]]; then echo 'JOBS must be positive' >&2; exit 2; fi
if ! test -f "${FEM_REPO}/src/core/lapack_eig.hpp"; then
  echo 'Set FEM_REPO to the local fem-eigen-em-2d3d-whitney checkout.' >&2; exit 2
fi
CAMPAIGN_DIR="${CAMPAIGN_DIR:-${ROOT}/docs/results/independent_validation/$(TZ=America/Sao_Paulo date +%F)_$(git -C "${ROOT}" rev-parse --short HEAD)}"
if test -e "${CAMPAIGN_DIR}"; then
  echo "Campaign already exists; choose a new CAMPAIGN_DIR: ${CAMPAIGN_DIR}" >&2; exit 2
fi
mkdir -p "${CAMPAIGN_DIR}/efg"
CAMPAIGN_DIR="$(cd "${CAMPAIGN_DIR}" && pwd)"
STAGING="$(mktemp -d "${TMPDIR:-/tmp}/m6a-fem.XXXXXX")"
echo "Campaign: ${CAMPAIGN_DIR}"
echo "Recoverable native FEM outputs: ${STAGING}"
# The original public M4 command runs the complete CTest suite first.
JOBS="${JOBS}" PROFILE=local BUILD_DIR="${BUILD_DIR}" OUT_DIR="${CAMPAIGN_DIR}/m4" \
  bash "${ROOT}/scripts/run_local_reference.sh" >"${CAMPAIGN_DIR}/efg_campaign.log" 2>&1
for shape in 9x7 13x11 17x15; do
  mv "${CAMPAIGN_DIR}/m4/local_${shape}" "${CAMPAIGN_DIR}/efg/${shape}"
done
cmake -S "${FEM_REPO}" -B "${STAGING}/build" -DCMAKE_BUILD_TYPE=Release >"${CAMPAIGN_DIR}/fem_configure.log" 2>&1
cmake --build "${STAGING}/build" --target helm10_rect --parallel "${JOBS}" >"${CAMPAIGN_DIR}/fem_build.log" 2>&1
cc -O2 -shared -fPIC "${ROOT}/scripts/fem_dsygv_capture.c" -o "${STAGING}/capture_lapack.so" -ldl -lm
for specification in 'coarse 16 8' 'medium 30 15' 'fine 40 20'; do
  read -r label nx ny <<<"${specification}"
  mkdir -p "${STAGING}/${label}/capture" "${STAGING}/${label}/raw"
  TP3485_OUT_DIR="${STAGING}/${label}/raw" \
  M6A_LAPACK_CAPTURE_DIR="${STAGING}/${label}/capture" \
  LD_PRELOAD="${STAGING}/capture_lapack.so" \
    "${STAGING}/build/helm10_rect" --ar-m 2.0 --nx "${nx}" --ny "${ny}" \
      --nmodos 12 --backend closed-form \
      >"${STAGING}/${label}/stdout.log" 2>"${STAGING}/${label}/stderr.log"
done
mkdir -p "${STAGING}/uninstrumented/raw"
TP3485_OUT_DIR="${STAGING}/uninstrumented/raw" \
  "${STAGING}/build/helm10_rect" --ar-m 2.0 --nx 16 --ny 8 --nmodos 12 --backend closed-form \
  >"${STAGING}/uninstrumented/stdout.log" 2>"${STAGING}/uninstrumented/stderr.log"
python3 "${ROOT}/scripts/prepare_fem_exports.py" --staging "${STAGING}" \
  --out "${CAMPAIGN_DIR}/fem" --fem-repo "${FEM_REPO}" | tee "${CAMPAIGN_DIR}/fem_export_checks.log"
python3 "${ROOT}/scripts/independent_validation.py" --campaign "${CAMPAIGN_DIR}" | tee "${CAMPAIGN_DIR}/analysis.log"
python3 "${ROOT}/scripts/independent_validation.py" --campaign "${CAMPAIGN_DIR}" --check-only | tee "${CAMPAIGN_DIR}/verification.log"
python3 "${ROOT}/scripts/validation_manifest.py" --campaign "${CAMPAIGN_DIR}" --fem-repo "${FEM_REPO}"
python3 "${ROOT}/scripts/validation_manifest.py" --campaign "${CAMPAIGN_DIR}" --verify
