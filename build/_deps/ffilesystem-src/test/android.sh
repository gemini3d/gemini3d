
case "$OSTYPE" in
  darwin*) ANDROID_SDK_ROOT_DEFAULT="$HOME/Library/Android/sdk" ;;
  linux*)  ANDROID_SDK_ROOT_DEFAULT="$HOME/Android/Sdk" ;;
  *) echo "Unsupported OS: $OSTYPE" >&2; return 1 2>/dev/null || exit 1 ;;
esac

export ANDROID_SDK_ROOT="${ANDROID_SDK_ROOT:-$ANDROID_SDK_ROOT_DEFAULT}"
export PATH="$ANDROID_SDK_ROOT/emulator:$ANDROID_SDK_ROOT/platform-tools:$ANDROID_SDK_ROOT/cmdline-tools/latest/bin:$PATH"
