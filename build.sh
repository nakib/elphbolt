#!/usr/bin/env bash
set -euo pipefail

# adjust these to the paths
MINIMAL_DIR="/home/sissa/Sally_apps/elphbolt/V3gpu"
CUDA_LIB="/usr/local/cuda/lib64"
RPATH="${MINIMAL_DIR}:${CUDA_LIB}"

echo "fpm clean"
fpm clean || true    # ignore error if nothing to clean

echo "removing build/"
rm -rf build

echo "fpm build"
fpm build \
  --link-flag "-L${MINIMAL_DIR}" \
  --link-flag "-L${CUDA_LIB}" \
  --link-flag "-Wl,-rpath,${RPATH}"

echo "fpm install"
fpm install \
  --link-flag "-L${MINIMAL_DIR}" \
  --link-flag "-L${CUDA_LIB}" \
  --link-flag "-Wl,-rpath,${RPATH}"

echo "Done."

