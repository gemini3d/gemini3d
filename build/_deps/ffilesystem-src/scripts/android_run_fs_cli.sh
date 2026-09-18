#!/usr/bin/env bash
set -euo pipefail

# Build and run fs_cli inside an Android emulator with interactive stdin/stdout.

usage() {
  cat <<'EOF'
Usage:
  scripts/android_run_fs_cli.sh [options] -- [fs_cli args...]

Options:
  -s, --serial <serial>      adb device serial (default: first emulator-* device)
  -b, --build-dir <dir>      build directory (default: build-android)
      --remote-path <path>   destination on device (default: /data/local/tmp/fs_cli)
      --skip-build           do not run cmake --build
  -h, --help                 show this help

Examples:
  scripts/android_run_fs_cli.sh
  scripts/android_run_fs_cli.sh -- --help
  scripts/android_run_fs_cli.sh -s emulator-5556 -- --version
EOF
}

case "${OSTYPE:-}" in
  darwin*) ANDROID_SDK_ROOT_DEFAULT="$HOME/Library/Android/sdk" ;;
  linux*)  ANDROID_SDK_ROOT_DEFAULT="$HOME/Android/Sdk" ;;
  *) echo "Unsupported OS: ${OSTYPE:-unknown}" >&2; exit 1 ;;
esac

export ANDROID_SDK_ROOT="${ANDROID_SDK_ROOT:-$ANDROID_SDK_ROOT_DEFAULT}"
export PATH="$ANDROID_SDK_ROOT/emulator:$ANDROID_SDK_ROOT/platform-tools:$ANDROID_SDK_ROOT/cmdline-tools/latest/bin:$PATH"

SERIAL=""
BUILD_DIR="build-android"
REMOTE_BIN="/data/local/tmp/app/fs_cli"
SKIP_BUILD=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    -s|--serial)
      SERIAL="$2"
      shift 2
      ;;
    -b|--build-dir)
      BUILD_DIR="$2"
      shift 2
      ;;
    --remote-path)
      REMOTE_BIN="$2"
      shift 2
      ;;
    --skip-build)
      SKIP_BUILD=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    --)
      shift
      break
      ;;
    *)
      break
      ;;
  esac
done

if ! command -v adb >/dev/null 2>&1; then
  echo "adb not found in PATH. Set ANDROID_SDK_ROOT or install platform-tools." >&2
  exit 1
fi

if [[ -z "$SERIAL" ]]; then
  SERIAL="$(adb devices | awk '/^emulator-[0-9]+[[:space:]]+device$/ {print $1; exit}')"
fi

if [[ -z "$SERIAL" ]]; then
  echo "No running emulator found. Start one, or pass --serial <device>." >&2
  exit 1
fi

if [[ "$SKIP_BUILD" -eq 0 ]]; then
  cmake --build "$BUILD_DIR" --target fs_cli
fi

LOCAL_BIN="$BUILD_DIR/app/fs_cli"
if [[ ! -x "$LOCAL_BIN" ]]; then
  if [[ -x "$BUILD_DIR/bin/app/fs_cli" ]]; then
    LOCAL_BIN="$BUILD_DIR/bin/app/fs_cli"
  else
    echo "Could not find built fs_cli in $BUILD_DIR." >&2
    exit 1
  fi
fi

adb -s "$SERIAL" push "$LOCAL_BIN" "$REMOTE_BIN" >/dev/null
adb -s "$SERIAL" shell chmod 755 "$REMOTE_BIN"

if [[ $# -eq 0 ]]; then
  echo "Binary pushed to $SERIAL:$REMOTE_BIN"
  echo "Opening interactive adb shell. Run: $REMOTE_BIN"
  exec adb -s "$SERIAL" shell
fi

quoted_args=""
for arg in "$@"; do
  quoted_args+=" $(printf '%q' "$arg")"
done

# Ask for a PTY when available so interactive stdin/stdout works reliably.
if adb -s "$SERIAL" shell -t true >/dev/null 2>&1; then
  exec adb -s "$SERIAL" shell -t "$REMOTE_BIN$quoted_args"
else
  exec adb -s "$SERIAL" shell "$REMOTE_BIN$quoted_args"
fi
