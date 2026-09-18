import argparse
import csv
import os
from pathlib import Path
import subprocess
import sys
import tempfile
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts/qualification'))
from check_transport import check_continuity
p=argparse.ArgumentParser();p.add_argument('--exe',type=Path,required=True);a=p.parse_args()
with tempfile.TemporaryDirectory() as d:
    proc=subprocess.run([str(a.exe),d],env=dict(os.environ,GEMINI_NUMERICAL_AUDIT='1'),capture_output=True,text=True)
    assert proc.returncode==0,proc.stdout+proc.stderr
    result=check_continuity(Path(d))
    assert result['passed'],result
    assert set(result['total_cleanup_number_by_species'].values())=={10.5},result
    print('Native cleanup adds 10.5 particles to each ion species across seven physical cells; null and ghost cells excluded; balance closes.')
