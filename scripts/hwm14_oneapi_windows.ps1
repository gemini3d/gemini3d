# first build with Intel oneAPI like
#  cmake --preset hwm14
#  cmake --build build
# then use this script
#  .\scripts\hwm14_oneapi_windows.ps1 build\gemini3d.run.exe ~\Downloads\eq_ground_hwm

$gemrun=(Resolve-Path $Args[0]).Path
$datadir=(Resolve-Path $Args[1]).Path


if (!(Test-Path $datadir -PathType Container)) {
    Write-Host "Data directory does not exist: $datadir"
    exit 1
}

if (!(Test-Path $gemrun -PathType Leaf)) {
    Write-Host "gemini3d.run program does not exist: $gemrun"
    exit 1
}

$builddir = (Split-Path -Parent $gemrun)

# we assume you've built HDF5 and ZLIB
$Env:Path = "$Env:Path;$builddir/_deps/hdf5_zlib-build;$builddir/_deps/hdf5-build/bin"

# for HWM14 need to execute from build directory so that the data files are found
Push-Location $builddir
try {
    Write-Host (Get-Location)
    .\gemini3d.run.exe
    .\gemini3d.run.exe $datadir
}
finally {
    Pop-Location
}
