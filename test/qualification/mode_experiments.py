"""Run bounded optional-mode checks without altering regression goldens.

Capacitance needs a 3D field-integrated domain and gradient boundary conditions.
For this smoke experiment its copied 3D input has zero imposed current/gradient
boundaries (flagdirich=0); it is a synthetic case, not a physical reference.
"""
import argparse
from pathlib import Path
import json
import re
import shutil
import subprocess
import time
import h5py
import numpy as np

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--build',type=Path,required=True)
p.add_argument('--inputs',type=Path,required=True)
p.add_argument('--inputs3d',type=Path)
p.add_argument('--work',type=Path,required=True)
p.add_argument('--mpiexec',default='mpiexec')
p.add_argument('--hwm',action='store_true')
a=p.parse_args()
if not a.hwm and a.inputs3d is None:p.error('--inputs3d is required for capacitance experiments')
a.work.mkdir(parents=True,exist_ok=False)
profile=[('fortran','gemini.bin',1,0,True),('cpp','gemini_c.bin',1,0,True)]
if not a.hwm:
    profile += [('altenergy','gemini.altenergy.bin',1,0,True),('field3_2d','gemini.bin',3,0,True),
                ('capacitance1_unsupported2d','gemini.bin',1,1,False),
                ('capacitance2_unsupported2d','gemini.bin',1,2,False),
                ('capacitance1_3d','gemini.bin',1,1,True),('capacitance2_3d','gemini.bin',1,2,True)]
results=[];outputs={}
for name,exe,potential,cap,succeeds in profile:
    case=a.work/name
    cap3d=cap and succeeds
    shutil.copytree(a.inputs3d if cap3d else a.inputs,case/'inputs')
    cfg=case/'inputs/config.nml';text=cfg.read_text()
    text=re.sub(r'(?im)^(tdur\s*=)[^!\n]*',r'\g<1> 120 ',text)
    text=re.sub(r'(?im)^(potsolve\s*=)[^!\n]*',r'\g<1> '+str(potential)+' ',text)
    if cap:text+='\n&capacitance\nflagcap='+str(cap)+'\nmagcap=5\n/\n'
    cfg.write_text(text)
    if cap3d:
        for path in (case/'inputs/Efield').glob('????????_?????.??????.h5'):
            with h5py.File(path,'r+') as f:
                f['flagdirich'][...]=0
                for key in f:
                    if key.startswith(('Vmin','Vmax')):f[key][...]=0
    start=time.monotonic()
    try:
        proc=subprocess.run([a.mpiexec,'-n','2',str(a.build.resolve()/exe),str(case.resolve())],
                            cwd=a.build,capture_output=True,text=True,timeout=240)
        code=proc.returncode;log=proc.stdout+proc.stderr
    except subprocess.TimeoutExpired as e:
        code=124;log=str(e)
    elapsed=time.monotonic()-start;(a.work/(name+'.log')).write_text(log)
    files=sorted(case.glob('????????_?????.??????.h5'));finite=False
    if code==0 and files:
        with h5py.File(files[-1]) as f:
            data={k:f['restart_core/'+k][...] for k in ['ns','Ts','vs1','Phi']}
        finite=all(np.isfinite(v).all() for v in data.values());outputs[name]=data
    passed=(code==0 and finite and len(files)==3) if succeeds else (code!=0 and code!=124 and 'capacitance: requires a 3D domain' in log)
    results.append(dict(mode=name,executable=exe,potential=potential,capacitance=cap,expected_success=succeeds,
                        returncode=code,finite=finite,elapsed_seconds=elapsed,frames=len(files),
                        passed=passed,independently_validated=False))
metrics=[]
if 'fortran' in outputs and 'cpp' in outputs:
    for key,av in outputs['fortran'].items():
        bv=outputs['cpp'][key]
        metrics.append(dict(field=key,relative_l2=float(np.linalg.norm(av-bv)/max(np.linalg.norm(av),1e-30)),
                            bitwise_equal=bool(np.array_equal(av,bv))))
payload=dict(schema='gemini.qualification.modes.2',hwm_enabled=a.hwm,results=results,c_cpp_comparison=metrics,
             scientific_mode_qualification=False,scope='120 s, two ranks on one node; execution and C/Fortran consistency only')
(a.work/'results.json').write_text(json.dumps(payload,indent=2)+'\n');print(json.dumps(payload,indent=2))
raise SystemExit(0 if all(r['passed'] for r in results) and all(m['bitwise_equal'] for m in metrics) else 1)
