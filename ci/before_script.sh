#!/usr/bin/env bash

set -e
set -o pipefail

if [[ "${CI_TARGET}" == lint ]]; then
  exit
fi

echo "${TRAVIS_BUILD_DIR}"
echo "${BUILD_DIR}"
echo "${CMAKE_FLAGS}"
mkdir build
cd build
cmake "${CMAKE_FLAGS}" ..
