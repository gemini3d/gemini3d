"""Native ETD checked against an 80-digit independent Decimal analytic oracle."""
import argparse
import decimal
import json
import math
from pathlib import Path
import subprocess


def main():
    p=argparse.ArgumentParser();p.add_argument('--exe',required=True)
    p.add_argument('--output',type=Path);a=p.parse_args()
    cases=[]
    for z in [0,1e-16,1e-12,1e-10,1.01e-10,1e-8,1e-6,.001,.009999,.01,.1,1,50]:
        for sign in [-1,1]:
            for f,source,dt in [(2,3,1),(1e12,1e4,.01),(-2,3,10),(1,-.1,1)]:
                cases.append((f,source,sign*z/dt,dt))
    # Equilibria across the branch threshold and tiny but dynamically important loss.
    cases += [(1e20,0,5e-11,1),(1e8,1,1e-8,1),(1,1,1,1)]
    proc=subprocess.run([a.exe],input=''.join(' '.join(map(str,c))+'\n' for c in cases),
                        text=True,capture_output=True,check=True)
    values=[float(x) for x in proc.stdout.split()]
    if len(values)!=len(cases):raise AssertionError('Missing native results')
    errors=[];old_failures=0
    with decimal.localcontext() as ctx:
        ctx.prec=80
        for case,value in zip(cases,values):
            f,p,l,dt=map(lambda x:decimal.Decimal(str(x)),case)
            exact=f+p*dt if l==0 else f*(-l*dt).exp()+p/l*(1-(-l*dt).exp())
            scale=max(abs(float(exact)),abs(float(f)),abs(float(p*dt)),1e-300)
            error=abs(value-float(exact))/scale
            if not math.isfinite(value) or error>5e-14:raise AssertionError((case,value,str(exact),error))
            z=float(l*dt)
            old=float(f)*math.exp(-z)+float(p/l)*(1-math.exp(-z)) if z>1e-10 else float(f+p*dt)
            old_failures+=abs(old-float(exact))/scale>5e-14
            errors.append(error)
    result=dict(schema='gemini.etd.oracle.1',passed=True,cases=len(cases),precision_digits=80,
                max_scaled_error=max(errors),budget=5e-14,old_formula_failing_cases=old_failures)
    if not old_failures:raise AssertionError('Test must detect the previous defect')
    if a.output:a.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))


if __name__=='__main__':main()
