param(
    [string]$Distribution = "Ubuntu",
    [switch]$SystemDeps,
    [string]$Root,
    [ValidateRange(1, 2147483647)]
    [int]$Jobs = 2,
    [ValidateSet("Debug", "Release")]
    [string]$BuildType = "Release",
    [switch]$ReferenceTests,
    [string]$SourceCache
)
$ErrorActionPreference = "Stop"
if (-not (Get-Command wsl.exe -ErrorAction SilentlyContinue)) {
    throw "Install WSL and Ubuntu first: wsl --install -d Ubuntu"
}
$scriptPath = Join-Path $PSScriptRoot "install-local.sh"
$linuxPath = & wsl.exe --distribution $Distribution --exec wslpath -a $scriptPath
if ($LASTEXITCODE -ne 0) { throw "Cannot locate this checkout in $Distribution." }
$arguments = @("--distribution", $Distribution, "--exec", "bash", $linuxPath.Trim())
$arguments += @("--jobs", "$Jobs", "--build-type", $BuildType)
if ($SystemDeps) { $arguments += "--system-deps" }
if ($Root) { $arguments += @("--root", $Root) }
if ($ReferenceTests) { $arguments += "--reference-tests" }
if ($SourceCache) { $arguments += @("--source-cache", $SourceCache) }
& wsl.exe @arguments
if ($LASTEXITCODE -ne 0) { throw "Local installation failed (exit $LASTEXITCODE)." }
