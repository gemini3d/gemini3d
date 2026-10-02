"""Require a controlled, stage-specific LAPACK failure, not a crash or timeout."""
import subprocess
import sys

stage = {"euler": "backEuler1D", "tr": "TRBDF21D TR stage", "bdf2": "TRBDF21D BDF2 stage"}[sys.argv[2]]
result = subprocess.run(sys.argv[1:], capture_output=True, text=True, timeout=20)
output = result.stdout + result.stderr
assert result.returncode > 0, (result.returncode, output)
assert f"PDEparabolic: {stage} gbsv INFO=2" in output, output
assert "PDEparabolic: singular diffusion matrix" in output, output
print(f"Singular {stage} rejected with LAPACK INFO=2")
