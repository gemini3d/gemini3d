#!/bin/sh
# a wrapper for CMAKE_CROSSCOMPILING_EMULATOR
# exit 77 means "skip test" (mapped in CMake via SKIP_RETURN_CODE)
#
# before using this script / tests, the Android SDK emulator must be running and connected via adb
# e.g. "adb devices" should show a device.


ADB="${ADB:-$ANDROID_SDK_ROOT/platform-tools/adb}"
RUN_DIR="${RUN_DIR:-/data/local/tmp/ffilesystem-tests}"

if [ ! -x "$ADB" ]; then
  echo "adb wrapper: adb not found at $ADB" >&2
  exit 77
fi

state="$("$ADB" get-state 2>/dev/null || true)"
if [ "$state" != "device" ]; then
  echo "adb wrapper: no connected Android device/emulator, skipping test." >&2
  exit 77
fi

binary="$1"; shift
base="$(basename "$binary")"
device_path="$RUN_DIR/$base"

"$ADB" shell "mkdir -p '$RUN_DIR'" >/dev/null 2>&1 || exit 77
"$ADB" push "$binary" "$device_path" >/dev/null || exit 77
"$ADB" shell "cd '$RUN_DIR' && export TMPDIR='$RUN_DIR' TMP='$RUN_DIR' TEMP='$RUN_DIR' HOME='$RUN_DIR' && chmod +x '$device_path' && '$device_path' $*"
