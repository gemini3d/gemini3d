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
    apt_update=("${elevate[@]}" apt-get -o Acquire::Retries=3 -o APT::Update::Error-Mode=any update)
    apt_update_ok=0
    for attempt in 1 2 3; do
      if "${apt_update[@]}"; then
        apt_update_ok=1
        break
      fi
      if [[ "$attempt" -lt 3 ]]; then
        sleep "$((attempt * 5))"
      fi
    done
    if [[ "$apt_update_ok" -eq 0 ]]; then
      allow_stale_apt_index="${GEMINI_ALLOW_STALE_APT_INDEX:-0}"
      allow_stale_apt_index="${allow_stale_apt_index//[[:space:]]/}"
      if [[ "$allow_stale_apt_index" != 1 ]]; then
        echo "apt-get update failed after retries; set GEMINI_ALLOW_STALE_APT_INDEX=1 to continue with existing package indexes" >&2
        exit 1
      fi
      has_packages=0
      has_inrelease=0
      has_release=0
      has_release_gpg=0
      compgen -G "/var/lib/apt/lists/*_Packages*" >/dev/null && has_packages=1
      compgen -G "/var/lib/apt/lists/*_InRelease" >/dev/null && has_inrelease=1
      compgen -G "/var/lib/apt/lists/*_Release" >/dev/null && has_release=1
      compgen -G "/var/lib/apt/lists/*_Release.gpg" >/dev/null && has_release_gpg=1
      has_release_metadata=0
      if [[ "$has_inrelease" -eq 1 || ( "$has_release" -eq 1 && "$has_release_gpg" -eq 1 ) ]]; then
        has_release_metadata=1
      fi
      if [[ "$has_packages" -eq 0 || "$has_release_metadata" -eq 0 ]]; then
        echo "apt-get update failed after retries and cached apt indexes are incomplete" >&2
        exit 1
      fi
      echo "warning: apt-get update failed after retries; proceeding with existing package indexes because GEMINI_ALLOW_STALE_APT_INDEX=1" >&2
    fi
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
