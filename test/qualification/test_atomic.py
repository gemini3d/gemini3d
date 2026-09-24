import argparse
from pathlib import Path
import subprocess
import tempfile
p=argparse.ArgumentParser();p.add_argument('--exe',required=True,type=Path);a=p.parse_args()
with tempfile.TemporaryDirectory() as d:
    root=Path(d);final=root/'checkpoint.h5';final.write_text('old valid checkpoint\n')
    (root/'directory').mkdir()
    for mode in ['leave','fail_publish','missing_parent']:
        r=subprocess.run([str(a.exe),d,mode],capture_output=True,text=True)
        assert (r.returncode==0)==(mode=='leave'),(mode,r.stderr)
        assert final.read_text()=='old valid checkpoint\n',mode
        if mode=='fail_publish':assert 'Atomic checkpoint publication failed' in r.stderr
    assert len(list(root.glob('checkpoint.h5.partial.*')))==2
    r=subprocess.run([str(a.exe),d,'success'],capture_output=True,text=True)
    assert r.returncode==0,r.stderr
    assert final.read_text()=='complete new checkpoint\n'
print('Atomic visibility: interrupted staging and failed publication preserve prior checkpoint; successful publish replaces it.')
