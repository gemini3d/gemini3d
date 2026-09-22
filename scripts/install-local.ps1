param(
    [string]$Distribution = "Ubuntu",
    [switch]$SystemDeps
)
$ErrorActionPreference = "Stop"
if (-not (Get-Command wsl.exe -ErrorAction SilentlyContinue)) {
    throw "Install WSL and Ubuntu first: wsl --install -d Ubuntu"
}
$scriptPath = Join-Path $PSScriptRoot "install-local.sh"
$linuxPath = & wsl.exe --distribution $Distribution --exec wslpath -a $scriptPath
if ($LASTEXITCODE -ne 0) { throw "Cannot locate this checkout in $Distribution." }
$arguments = @("--distribution", $Distribution, "--exec", "bash", $linuxPath.Trim())
if ($SystemDeps) { $arguments += "--system-deps" }
& wsl.exe @arguments
if ($LASTEXITCODE -ne 0) { throw "Local installation failed (exit $LASTEXITCODE)." }
