"""Unsupported executables must reject exact-restart qualification before work."""
import argparse
import json
import os
import subprocess

p=argparse.ArgumentParser();p.add_argument('--alternate',required=True);p.add_argument('--density',required=True);a=p.parse_args()
results=[]
for exe,model in [(a.alternate,'alternate energy'),(a.density,'density and potential only')]:
    proc=subprocess.run([exe],env=dict(os.environ,GEMINI_EXACT_RESTART='1'),text=True,capture_output=True,timeout=20)
    passed=proc.returncode!=0 and ('Exact restart research profile excludes model: '+model) in proc.stdout+proc.stderr
    results.append(dict(model=model,returncode=proc.returncode,passed=passed))
print(json.dumps(results))
raise SystemExit(0 if all(r['passed'] for r in results) else 1)
