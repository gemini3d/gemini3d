#!/usr/bin/env bash
# System-wide prerequisites are installed only by explicit user request.
set -euo pipefail

case "$(uname -s)" in
  Linux)
    if ! command -v apt-get >/dev/null; then
      echo "Automatic setup supports Debian/Ubuntu (including WSL). On other Linux systems install GCC, MPI, HDF5, LAPACK/ScaLAPACK and Python 3.11+ with venv, then rerun without --system-deps." >&2
      exit 1
    fi
    elevate=()
    if [[ "$(id -u)" != 0 ]]; then
      command -v sudo >/dev/null || { echo "sudo is required for system packages." >&2; exit 1; }
      elevate=(sudo)
    fi
    "${elevate[@]}" apt-get update
    "${elevate[@]}" apt-get install -y --no-install-recommends \
      build-essential gfortran git python3 python3-venv python3-dev \
      libopenmpi-dev openmpi-bin libhdf5-dev libopenblas-dev \
      libscalapack-openmpi-dev zlib1g-dev
    ;;
  Darwin)
    command -v brew >/dev/null || { echo "Install Homebrew from https://brew.sh first." >&2; exit 1; }
    xcode-select -p >/dev/null || { echo "Run xcode-select --install first." >&2; exit 1; }
    brew install gcc git python@3.12 open-mpi hdf5 openblas scalapack
    ;;
  *)
    echo "Use Windows WSL with Ubuntu, or provision a supported Unix toolchain manually." >&2
    exit 1
    ;;
esac
