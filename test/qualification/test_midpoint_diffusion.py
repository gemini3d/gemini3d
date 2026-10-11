"""Run native operator, verify observed stages, then compare independent laws."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import tempfile

p=argparse.ArgumentParser()
p.add_argument('--source',type=Path,required=True)
p.add_argument('--build',type=Path,required=True)
p.add_argument('--work',type=Path,required=True)
a=p.parse_args()
# Keep each invocation isolated, including concurrent or repeated CTest calls.
# No existing evidence or user directory is removed.
a.work.mkdir(parents=True,exist_ok=True)
invocation=Path(tempfile.mkdtemp(prefix='run-',dir=a.work))
native=invocation/'native';trace=invocation/'trace'
scripts=a.source/'scripts/qualification'
subprocess.run([sys.executable,str(scripts/'run_physics_native.py'),'--build',str(a.build),
    '--source',str(a.source/'test/qualification/midpoint_diffusion_probe.f90'),'--work',str(native)],check=True)
subprocess.run([sys.executable,str(scripts/'trace_diffusion_stages.py'),
    '--source',str(a.source/'src/numerical/diffusion/PDEparabolic.f90'),
    '--build',str(a.build),'--probe',str(native/'physics_native.o'),'--work',str(trace)],check=True)
try:
    revision=subprocess.check_output(['git','-C',str(a.source),'rev-parse','HEAD'],text=True).strip()
except subprocess.CalledProcessError:
    revision='source-export; see PACKAGE_PROVENANCE.json'
subprocess.run([sys.executable,str(a.source/'test/qualification/analyze_midpoint_diffusion.py'),
    '--work',str(native),'--trace',str(trace/'tr-stages.csv'),'--candidate',revision],check=True)
result=json.loads((native/'physics-results.json').read_text())
a.work.joinpath('latest-result.json').write_text(json.dumps({'invocation':str(invocation),
    'result':result},indent=2)+'\n')
assert result['passed'] and all(r['native_stage_trace_used'] and r['fixed_boundary_max_error_K']<=1e-11
    and r['floor_cell_events']==0 for r in result['PC03']['cases'])
