"""Independent coupled ODE oracle; acceptance does not validate rate amplitudes."""
import argparse
import csv
import hashlib
import json
import subprocess
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

p = argparse.ArgumentParser()
p.add_argument('--exe', type=Path, required=True)
p.add_argument('--output', type=Path, required=True)
a = p.parse_args()
a.output.parent.mkdir(parents=True, exist_ok=True)
csv_path = a.output.with_suffix('.csv')
subprocess.run([str(a.exe), str(csv_path)], check=True)
with csv_path.open() as f:
    rows = [{k: float(v) for k, v in row.items()} for row in csv.DictReader(f)]
assert len(rows) == 381 and all(np.isfinite(list(r.values())).all() for r in rows)
results = []
for r in rows:
    ca, cb, ra, rb, drag, k, dt, fraction = (r[x] for x in
        ('Ca', 'Cb', 'rhoa', 'rhob', 'drag', 'K', 'dt', 'fraction_a'))
    initial = np.array([r[x] for x in ('va0x','va0y','va0z','vb0x','vb0y','vb0z','Ta0','Tb0')])
    observed = np.array([r[x] for x in ('vax','vay','vaz','vbx','vby','vbz','Ta','Tb')])
    def rhs(t, y):
        force = drag*(y[3:6]-y[:3])
        heat = drag*np.dot(y[:3]-y[3:6], y[:3]-y[3:6])
        transfer = k*(y[7]-y[6])
        return dt*np.r_[force/ra, -force/rb,
            (transfer+fraction*heat)/ca, (-transfer+(1-fraction)*heat)/cb]
    if dt:
        oracle = solve_ivp(rhs, (0.,1.), initial, method='Radau', rtol=2e-12, atol=2e-11)
        control = solve_ivp(rhs, (0.,1.), initial, method='Radau', rtol=2e-13, atol=2e-12)
        assert oracle.success and control.success
        expected = control.y[:,-1]
        control_error = np.max(abs(oracle.y[:,-1]-expected)/np.maximum(abs(expected), 1.))
    else:
        expected=initial;control_error=0.
        # The native probe checks U and v bitwise at dt=0. CSV temperatures
        # include a rounded U/Cv conversion even when the operator is a no-op.
    error = np.max(abs(observed-expected)/np.maximum(abs(expected),1.))
    e0=ca*initial[6]+cb*initial[7]+ra*np.dot(initial[:3],initial[:3])/2+rb*np.dot(initial[3:6],initial[3:6])/2
    e1=ca*observed[6]+cb*observed[7]+ra*np.dot(observed[:3],observed[:3])/2+rb*np.dot(observed[3:6],observed[3:6])/2
    energy=abs(e1-e0)/abs(e0)
    momentum=np.max(abs(ra*(observed[:3]-initial[:3])+rb*(observed[3:6]-initial[3:6])))/max(
        ra*np.linalg.norm(initial[:3])+rb*np.linalg.norm(initial[3:6]),np.finfo(float).tiny)
    identity=max(r[x] for x in ('instantaneous_momentum','instantaneous_energy','instantaneous_thermal'))
    passed=error<=1e-10 and control_error<=1e-10 and energy<=1e-11 and momentum<=1e-11 and identity<=1e-11
    results.append({'species':[int(r['a']),int(r['b'])], 'level':int(r['level']),
        'oracle_relative_error':float(error), 'oracle_control_error':float(control_error),
        'energy_relative_error':float(energy), 'momentum_relative_error':float(momentum),
        'instantaneous_identity_error':float(identity), 'passed':bool(passed)})
result={'schema':'gemini.pair.exchange.validation.1', 'scope':'closed constant-rate pair operator; application source coupling unexecuted',
    'rate_magnitudes_independently_validated':False, 'case_count':len(results),
    'passed':all(r['passed'] for r in results), 'cases':results,
    'executable_sha256':hashlib.sha256(a.exe.read_bytes()).hexdigest(),
    'csv_sha256':hashlib.sha256(csv_path.read_bytes()).hexdigest()}
a.output.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='cases'}))
raise SystemExit(0 if result['passed'] else 1)
