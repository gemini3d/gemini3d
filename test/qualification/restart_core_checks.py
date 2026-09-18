"""Adversarial native restart checks, including legacy-data poisoning.

Uses a completed 300-second mini2dns_fang case. All cases and logs are isolated.
The acceptance bounds test serialization and parser behavior, not physics.
"""
import argparse
from pathlib import Path
import json
import re
import shutil
import subprocess
import h5py
import numpy as np
p=argparse.ArgumentParser();p.add_argument('--exe',type=Path,required=True);p.add_argument('--case',type=Path,required=True)
p.add_argument('--work',type=Path,required=True);p.add_argument('--mpiexec',default='mpiexec');a=p.parse_args()
a.work.mkdir(parents=True,exist_ok=False)
frames=sorted(a.case.glob('????????_?????.??????.h5'));assert len(frames)>=4
saved=frames[2].name;final=frames[3].name
# Truncate at an interior output time with a subsequent valid forcing frame.
with h5py.File(frames[3]) as f:final_ut=float(f['time/UThour'][()])*3600
text=(a.case/'inputs/config.nml').read_text();start=float(re.search(r'(?im)^UTsec0\s*=\s*([\d.]+)',text)[1])
duration=final_ut-start
records=[];outputs={}
for label in ['clean','poison_legacy','incomplete_core','schema','potential_rank','false_precision']:
    case=a.work/label;shutil.copytree(a.case/'inputs',case/'inputs')
    for f in frames[:3]:shutil.copy2(f,case/f.name)
    cfg=case/'inputs/config.nml';cfg.write_text(re.sub(r'(?im)^(tdur\s*=)[^!\n]*',r'\g<1> '+str(duration)+' ',text))
    with h5py.File(case/saved,'r+') as f:
        if label=='poison_legacy':
            for key in ['nsall','Tsall','vs1all','Phiall']:f[key][...]=np.nan
        elif label=='incomplete_core':del f['restart_core/complete']
        elif label=='schema':f['restart_core/schema'][...]=999
        elif label=='potential_rank':
            data=f['Phiall'][...];del f['restart_core/Phi'];f.create_dataset('restart_core/Phi',data=data)
        elif label=='false_precision':
            data=f['restart_core/ns'][...].astype('float32')
            del f['restart_core/ns'];f.create_dataset('restart_core/ns',data=data)
    proc=subprocess.run([a.mpiexec,'-n','2',str(a.exe.resolve()),str(case.resolve())],text=True,capture_output=True,timeout=240)
    (a.work/(label+'.log')).write_text(proc.stdout+proc.stderr)
    expected={'incomplete_core':'Incomplete core restart','schema':'Unsupported core restart',
              'potential_rank':'must be full 3D','false_precision':'must contain float64'}
    passed=proc.returncode==0 if label not in expected else proc.returncode!=0 and expected[label] in proc.stdout+proc.stderr
    records.append(dict(case=label,returncode=proc.returncode,passed=passed))
    if label not in expected and passed:
        with h5py.File(case/final) as f:outputs[label]={k:f['restart_core/'+k][...] for k in ['ns','Ts','vs1','Phi']}
exact=False
if len(outputs)==2:exact=all(np.array_equal(outputs['clean'][k],outputs['poison_legacy'][k]) for k in outputs['clean'])
with h5py.File(a.case/saved) as f:
    precision={k:str(f['restart_core/'+k].dtype) for k in ['ns','Ts','vs1','Phi']}
result=dict(schema='gemini.qualification.restart_core.1',runs=records,legacy_poison_bitwise_equal=exact,dtypes=precision,
            passed=all(r['passed'] for r in records) and exact and set(precision.values())=={'float64'},
            scope='core record precision, preference, and invalid-file rejection; full auxiliary state not certified')
(a.work/'results.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
raise SystemExit(0 if result['passed'] else 1)
