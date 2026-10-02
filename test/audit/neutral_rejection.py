"""Require the physical-policy diagnostic, not a bounds error or floating-point trap."""
import subprocess
import sys

messages = {
    "lower_coverage": "insufficient lower-end altitude coverage",
    "upper_coverage": "insufficient upper-end altitude coverage",
    "zero_denominator": "zero lower extrapolation density",
    "negative_density": "negative extrapolation density",
    "underground": "field line is entirely below ground",
    "temperature": "temperatures must be finite and positive",
    "singleton_uncovered": "insufficient altitude coverage",
}
exe, mode = sys.argv[1:]
result = subprocess.run([exe, mode], capture_output=True, text=True, timeout=30)
output = result.stdout + result.stderr
assert result.returncode > 0, (mode, result.returncode, output)
assert messages[mode] in output, (mode, output)
assert "Fortran runtime error" not in output, (mode, output)
print(f"neutral background rejected {mode} with the expected diagnostic")
