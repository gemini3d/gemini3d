#!/usr/bin/env python3
"""
Count the number of Fortran procedures (functions and subroutines) in a given Fortran source file.

example:
    python scripts/count_fortran_procedures.py src/fortran/filesystem.F90
"""

import re
from pathlib import Path

def count_fortran_procedures(file_path: str | Path) -> int:
    content = Path(file_path).read_text()
    # Match function/subroutine declarations, allowing type/attribute prefixes.
    # Example: "logical(C_BOOL) function is_rosetta() bind(C, ...)"
    proc_pattern = re.compile(
        r'^\s*'
        r'(?!\s*!)(?!\s*end\s+(?:function|subroutine)\b)'
        r'.*?\b(function|subroutine)\s+([A-Za-z][A-Za-z0-9_]*)\b',
        re.IGNORECASE,
    )
    matches: list[tuple[str, str]] = []

    for line in content.splitlines():
        # Ignore trailing comments when checking procedure starts.
        code = line.split('!', 1)[0]

        match = proc_pattern.match(code)
        if match:
            matches.append((match.group(1), match.group(2)))

    # Count unique procedure names (case-insensitive) to avoid double-counting
    # interface declarations and repeated signatures.
    unique_matches: list[tuple[str, str]] = []
    seen_names: set[str] = set()
    for kind, name in matches:
        key = name.lower()
        if key in seen_names:
            continue
        seen_names.add(key)
        unique_matches.append((kind, name))

    # print procedure names for debugging
    for kind, name in unique_matches:
        print(f"{kind} {name}")
    return len(unique_matches)

if __name__ == "__main__":
    print(count_fortran_procedures('src/fortran/filesystem.F90'))
