#!/usr/bin/env bash
set -euo pipefail

MINIMAL_DIR="/home/sissa/Sally_apps/elphbolt/V3gpu"
CUDA_LIB="/usr/local/cuda/lib64"
NVHPC_COMPILERS="/home/sissa/nvidia/hpc_sdk/Linux_x86_64/25.7/compilers/lib"
NVHPC_CUDA="/home/sissa/nvidia/hpc_sdk/Linux_x86_64/25.7/cuda/12.9/lib64"
NVHPC_MATH="/home/sissa/nvidia/hpc_sdk/Linux_x86_64/25.7/math_libs/12.9/lib64"
GCC_LIBGOMP_DIR="/usr/lib/gcc/x86_64-linux-gnu/12"
SYSTEM_LIB_DIR="/usr/lib/x86_64-linux-gnu"
RPATH="${MINIMAL_DIR}:${CUDA_LIB}:${NVHPC_COMPILERS}:${NVHPC_CUDA}:${NVHPC_MATH}"

echo "fpm clean"
fpm clean || true

echo "removing build/"
rm -rf build

echo "fpm build"
fpm build \
  --link-flag "-L${GCC_LIBGOMP_DIR}" \
  --link-flag "-L${SYSTEM_LIB_DIR}" \
  --link-flag "-L${MINIMAL_DIR}" \
  --link-flag "-L${CUDA_LIB}" \
  --link-flag "-L${NVHPC_COMPILERS}" \
  --link-flag "-L${NVHPC_CUDA}" \
  --link-flag "-L${NVHPC_MATH}" \
  --link-flag "-Wl,-rpath,${RPATH}" \
  --link-flag "-Wl,--copy-dt-needed-entries"

echo "fpm install"
fpm install \
  --link-flag "-L${GCC_LIBGOMP_DIR}" \
  --link-flag "-L${SYSTEM_LIB_DIR}" \
  --link-flag "-L${MINIMAL_DIR}" \
  --link-flag "-L${CUDA_LIB}" \
  --link-flag "-L${NVHPC_COMPILERS}" \
  --link-flag "-L${NVHPC_CUDA}" \
  --link-flag "-L${NVHPC_MATH}" \
  --link-flag "-Wl,-rpath,${RPATH}" \
  --link-flag "-Wl,--copy-dt-needed-entries"

echo "Done."

