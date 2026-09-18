"""Compare independent spherical harmonics to official IAGA/NOAA Fortran."""
import argparse
import importlib.util
import itertools
import json
from pathlib import Path
import subprocess
import math

p=argparse.ArgumentParser();p.add_argument('--exe',required=True);p.add_argument('--output',type=Path,required=True)
a=p.parse_args();root=Path(__file__).resolve().parents[2]
spec=importlib.util.spec_from_file_location('igrf14',root/'scripts/qualification/igrf14.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
rows=module.coefficients(root/'test/qualification/data/igrf14coeffs.txt')
cases=list(itertools.product([2025,2026,2027.5,2029.75,2030], [6371.2,6471.2,6671.2,7371.2],
                             [-89,-70,-30,0,30,70,89],[-180,-90,0,45,179,360]))
text=''.join(f'{y} {r} {90-lat} {lon%360}\n' for y,r,lat,lon in cases)
proc=subprocess.run([a.exe],input=text,text=True,capture_output=True,check=True)
ref=[list(map(float,line.split())) for line in proc.stdout.splitlines() if line.strip()]
if len(ref)!=len(cases):raise ValueError('Official oracle did not return each case')
records=[]
# Official code uses a truncated degree/radian constant; 0.05 nT is far
# below 0.1 nT coefficient rounding, and is fixed before viewing results.
limit_nt=0.05
for case,truth in zip(cases,ref):
    value=module.field(rows,*case);error=max(abs(a-b) for a,b in zip(value,truth))
    if not math.isfinite(error):raise ValueError('Nonfinite field residual')
    records.append(dict(input=case,max_component_error_nt=error))
rejected=0
for c in [(2031,6671.2,0,0),(2026,0,0,0),(2026,6671.2,90,0),(float('nan'),6671.2,0,0)]:
    try: module.field(rows,*c)
    except ValueError:rejected+=1
payload=dict(schema='gemini.qualification.igrf14.1',cases=len(records),rejected_inputs=rejected,
             acceptance_max_component_error_nt=limit_nt,max_component_error_nt=max(x['max_component_error_nt'] for x in records),
             production_grid_qualified=False,scope='geocentric field oracle, not coordinated grid/basis/driver transforms',records=records)
payload['passed']=payload['max_component_error_nt']<=limit_nt and rejected==4
a.output.write_text(json.dumps(payload,indent=2)+'\n')
print(json.dumps({k:v for k,v in payload.items() if k!='records'},indent=2))
raise SystemExit(0 if payload['passed'] else 1)
