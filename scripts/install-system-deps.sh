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
    has_cached_apt_index_set() {
      local package_index
      local release_prefix
      local found_release_metadata
      local has_packages=0
      local has_consistent_metadata=1
      while IFS= read -r package_index; do
        release_prefix="${package_index##*/}"
        case "$release_prefix" in
          *_Packages|*_Packages.lz4|*_Packages.xz|*_Packages.gz|*_Packages.bz2|*_Packages.zst) ;;
          *) continue ;;
        esac
        has_packages=1
        release_prefix="${release_prefix%_Packages*}"
        release_prefix="${release_prefix%_binary-*}"
        found_release_metadata=0
        while [[ "$release_prefix" == *"_dists_"* && "$release_prefix" == *_* ]]; do
          if [[ -f "/var/lib/apt/lists/${release_prefix}_InRelease" ]]; then
            found_release_metadata=1
            break
          fi
          if [[ -f "/var/lib/apt/lists/${release_prefix}_Release" ]]; then
            found_release_metadata=1
            break
          fi
          release_prefix="${release_prefix%_*}"
        done
        if [[ "$found_release_metadata" -eq 0 ]]; then
          has_consistent_metadata=0
          break
        fi
      done < <(compgen -G "/var/lib/apt/lists/*_Packages*")
      if [[ "$has_packages" -eq 0 || "$has_consistent_metadata" -eq 0 ]]; then
        return 1
      fi
      return 0
    }
    run_apt_update_with_fallback() {
      local -a elevate_cmd=("$@")
      local -a apt_update=("${elevate_cmd[@]}" apt-get -o APT::Update::Error-Mode=any update)
      local apt_update_ok=0
      local attempt
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
        local allow_stale_apt_index
        allow_stale_apt_index="${GEMINI_ALLOW_STALE_APT_INDEX:-0}"
        allow_stale_apt_index="${allow_stale_apt_index//[[:space:]]/}"
        if [[ "$allow_stale_apt_index" != 1 ]]; then
          echo "apt-get update failed after retries; set GEMINI_ALLOW_STALE_APT_INDEX=1 (whitespace is ignored) to continue with existing package indexes" >&2
          return 1
        fi
        if ! has_cached_apt_index_set; then
          echo "apt-get update failed after retries and cached apt indexes are incomplete" >&2
          return 1
        fi
        echo "warning: apt-get update failed after retries; proceeding with existing package indexes because GEMINI_ALLOW_STALE_APT_INDEX=1" >&2
      fi
    }
    run_apt_update_with_fallback "${elevate[@]}"
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
