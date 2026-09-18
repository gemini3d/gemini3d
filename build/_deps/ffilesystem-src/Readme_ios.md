# Ffilesystem for iPhone / iOS

This project can be cross-compiled for iOS using CMake on macOS.
For a tech demo build, compile library and optional CLI targets without packaging/signing.

## Build presets

### Device build (iphoneos, arm64)

```sh
cmake --workflow ios
```

### Simulator build (iphonesimulator, arm64)

```sh
cmake --workflow ios-sim
```

## Why direct execution can fail

If you build with `iphoneos` (device SDK), running the binary directly from macOS shell can fail, for example:

```
zsh: killed build-ios/fs_cli.app/fs_cli
```

A device-target iOS binary is not a normal macOS executable.

Use:

1. `ios` preset for iPhone/iPad device targets.
2. `ios-sim` preset for Simulator targets.
3. host `default` preset to run directly in Terminal.

## Simulator runtime prerequisite

If this shows no available devices:

```sh
xcrun simctl list devices available
```

install an iOS Simulator runtime first:

1. Open Xcode.
2. Go to Xcode > Settings > Platforms.
3. Download at least one iOS Simulator Runtime.
4. Re-run `xcrun simctl list devices available`.

If needed, also ensure CLI tools point to full Xcode and first-launch tasks are complete:

```sh
xcode-select -s /Applications/Xcode.app/Contents/Developer
xcodebuild -license accept
xcodebuild -runFirstLaunch
```

Recent Xcode may also support:

```sh
xcodebuild -downloadPlatform iOS
```

## Minimal simulator launch demo

1. Build simulator app:

```sh
cmake --workflow ios-sim
```

2. List available simulator devices:

```sh
xcrun simctl list devices available
```

3. Boot a device (example name):

```sh
xcrun simctl boot "iPhone 16"
xcrun simctl bootstatus "iPhone 16" -b
```

4. Install app bundle:

```sh
xcrun simctl install booted build-ios-sim/fs_cli.app
```

5. Get bundle identifier:

```sh
/usr/libexec/PlistBuddy -c "Print :CFBundleIdentifier" build-ios-sim/fs_cli.app/Info.plist
```

6. Launch app and attach console:

```sh
xcrun simctl launch --console-pty booted <bundle-id>
```

Note: simulator launch is app-style execution, so interactive REPL behavior can differ from a normal macOS terminal executable.
