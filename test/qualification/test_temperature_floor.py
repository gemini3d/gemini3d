"""Measure the conversion floor against known physical volume and heat capacity."""
import argparse,csv,os
from pathlib import Path
import subprocess,tempfile

p=argparse.ArgumentParser();p.add_argument('--exe',type=Path,required=True);a=p.parse_args()
with tempfile.TemporaryDirectory() as d:
    p=Path(d);audited=p/'audited';plain=p/'plain';audited.mkdir();plain.mkdir()
    outputs=[]
    for directory,flag in [(audited,'1'),(plain,'0')]:
        run=subprocess.run([str(a.exe),str(directory)],env=dict(os.environ,GEMINI_NUMERICAL_AUDIT=flag),capture_output=True,text=True)
        assert run.returncode==0,run.stdout+run.stderr
        outputs.append(run.stdout)
    assert outputs[0]==outputs[1],'Diagnostics changed native oracle output'
    assert (audited/'temperature.bin').read_bytes()==(plain/'temperature.bin').read_bytes(),'Diagnostics changed native state'
    expected={int(line.split()[0]):float(line.split()[1]) for line in outputs[0].splitlines()}
    with (audited/'temperature-floor-r00000000.csv').open() as f:rows=list(csv.DictReader(f))
    assert len(rows)==7 and {int(r['stage']) for r in rows}=={3},rows
    reported={int(r['species']):float(r['energy_added_J']) for r in rows}
    assert set(reported)==set(expected)==set(range(1,8)),rows
    assert all(abs(reported[s]-e)<=1e-11*abs(e) for s,e in expected.items()),rows
    print('PASS: seven species; 50-to-100 K conversion floor; 183 m3 physical volume; null/ghost exclusions; identical state with diagnostics enabled and disabled.')
