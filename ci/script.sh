#!/usr/bin/env bash

set -e
set -o pipefail

pwd
cd "${BUILD_DIR}"
pwd

if [[ "${CI_TARGET}" == lint ]]; then
  exit
fi

if [[ "${CI_TARGET}" == unit ]]; then
  make -j2 unit
fi

if [[ "${CI_TARGET}" == test ]]; then
  make -j2 && ctest --output-on-failure
fi

