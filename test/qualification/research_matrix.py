"""Execute the frozen scope: driver preflight, restart trajectories, and transport.

All outputs are fresh. Retains diagnostics for expected rejections. No budget is
learned from these results. Requires the pinned native reference inputs.
"""
import argparse
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time
import h5py
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts/qualification'))
from run_research import input_signature,preflight
from check_transport import check,check_continuity
from check_sources import check as check_sources
from check_energy_operators import check as check_energy


def set_duration(case,duration,*,extend=False):
    cfg=case/'inputs/config.nml';text=cfg.read_text()
    pattern=re.compile(r'(?im)^([ \t]*tdur[ \t]*=)([^!\r\n]*)')
    matches=list(pattern.finditer(text))
    if len(matches)!=1:raise ValueError('Exactly one tdur assignment required')
    try:
        value=matches[0][2].strip().removesuffix(',').strip()
        previous=float(value.replace('D','e').replace('d','e'))
        duration=float(duration)
    except (TypeError,ValueError) as e:
        raise ValueError('tdur must be a positive finite duration') from e
    if not all(math.isfinite(value) and value>0 for value in (previous,duration)):
        raise ValueError('tdur must be a positive finite duration')
    if extend and duration<=previous:raise ValueError('Restart tdur must extend the first run')
    cfg.write_text(pattern.sub(lambda m:m[1]+f' {duration:.17g} ',text))


def restart_advanced(case,previous):
    frames={p.name for p in case.glob('????????_?????.??????.h5')}
    return bool(previous) and previous<frames and max(frames)>max(previous)


def checkpoint_shapes(case,ranks,layout):
    """Check identity, finite payloads and the C/Fortran allocation contract."""
    frames=sorted(case.glob('????????_?????.??????.h5'))
    if not frames:return False
    for frame in frames:
        with h5py.File(frame) as f:
            root=f['restart_runtime']
            prefix=root['prefix'][()].decode().strip()
            checkpoint=root['checkpoint'][()]
            input_sha=root['input_sha256'][()]
            executable_sha=root['executable_sha256'][()]
            if root['schema'][()]!=2 or tuple(root['layout'][...])!=tuple(layout):return False
            if checkpoint.decode().strip()!=frame.name or not prefix.startswith(frame.name+'.partial.'):return False
        for rank in range(ranks):
            with h5py.File(case/(prefix+f'.r{rank:08d}.h5')) as f:
                if f['schema'][()]!=2 or f['complete'][()]!=1 or f['rank'][()]!=rank:return False
                if tuple(f['layout'][...])!=tuple(layout) or f['checkpoint'][()]!=checkpoint:return False
                if f['generation'][()].decode().strip()!=prefix:return False
                if f['input_sha256'][()]!=input_sha or f['executable_sha256'][()]!=executable_sha:return False
                if f['fluid'].shape[0]!=35 or f['electro'].shape[0]!=7:return False
                if f['fluid'].shape[1:]!=f['electro'].shape[1:]:return False
                for field in ('fluid','electro','vi1','vi2','vi3'):
                    if f[field].ndim!=4 or f[field].dtype!=np.dtype('float64'):return False
                    if not np.isfinite(f[field][...]).all():return False
    return True


def compare(continuous,split,budget):
    rows=[]
    frames=sorted(continuous.glob('????????_?????.??????.h5'))
    if not frames or {f.name for f in frames}!={f.name for f in split.glob('????????_?????.??????.h5')}:
        return dict(passed=False,fields_checked=0,max_budget_fraction=None,error='Missing or unexpected output frame')
    for fa in frames:
        fb=split/fa.name
        with h5py.File(fa) as a,h5py.File(fb) as b:
            # Require the complete output sequence, including all continuation frames.
            for key in ['ns','Ts','vs1','Phi']:
                av=a['restart_core/'+key][...].astype(float);bv=b['restart_core/'+key][...].astype(float)
                if av.shape!=bv.shape or not np.isfinite(av).all():
                    return dict(passed=False,fields_checked=len(rows),max_budget_fraction=None,error='Invalid reference or output shape')
                for species in range(7) if key!='Phi' else [None]:
                    x=av[species] if species is not None else av;y=bv[species] if species is not None else bv
                    error=float(np.sqrt(np.mean((x-y)**2)));scale=float(np.sqrt(np.mean(x*x)))
                    bound=budget['absolute_floor'][key]+budget['trajectory_relative_l2']*scale
                    rows.append(dict(frame=fa.name,field=key,species=species,rms_error=error,bound=bound,
                                     max_abs=float(np.max(np.abs(x-y))),finite=bool(np.isfinite(y).all()),
                                     passed=bool(np.isfinite(y).all() and error<=bound)))
    return dict(passed=bool(rows) and all(r['passed'] for r in rows),fields_checked=len(rows),
                max_budget_fraction=max((r['rms_error']/r['bound'] for r in rows),default=None),metrics=rows)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--build',type=Path,required=True)
    p.add_argument('--inputs-root',type=Path,required=True);p.add_argument('--work',type=Path,required=True)
    p.add_argument('--mpiexec',default='mpiexec');p.add_argument('--quick',action='store_true');p.add_argument('--suffix',default='')
    p.add_argument('--extended-layouts',action='store_true',
                   help='also test four-rank 2x2 decomposition outside the frozen 1/2-rank profile')
    a=p.parse_args()
    if a.quick and a.extended_layouts:p.error('--extended-layouts requires the full matrix, not --quick')
    a.work.mkdir(parents=True,exist_ok=False)
    root=Path(__file__).resolve().parents[2]
    budget=json.loads((root/'docs/qualification/ACCEPTANCE_BUDGETS.json').read_text())['budgets']['checkpoint']
    exe=(a.build/('gemini.bin'+a.suffix)).resolve();exe_sha=hashlib.sha256(exe.read_bytes()).hexdigest()
    runs=[];comparisons=[];transport=[];continuity=[];sources=[];energy=[]
    cases=[('mini2dns_fang',1,[1,1]),('mini2dns_fang',2,[1,2]),('mini3d_fang',1,[1,1]),
           ('mini3d_fang',2,[1,2]),('mini3d_fang',2,[2,1])]
    if a.extended_layouts:cases.append(('mini3d_fang',4,[2,2]))
    if a.quick:cases=cases[:2]
    def create(name,case,duration):
        path=a.work/name;shutil.copytree(a.inputs_root/case/'inputs',path/'inputs')
        set_duration(path,duration)
        return path
    def run(path,label,ranks,layout,expected=None,executable=exe,restart=False):
        valid=preflight(path);(path/'forcing-preflight.json').write_text(json.dumps(valid,indent=2)+'\n')
        sig=input_signature(path);(path/'verified-inputs.json').write_text(json.dumps(sig,indent=2)+'\n')
        env=dict(os.environ,GEMINI_EXACT_RESTART='1',GEMINI_INPUT_SHA256=sig['sha256'],
                 GEMINI_EXECUTABLE_SHA256=hashlib.sha256(executable.read_bytes()).hexdigest(),GEMINI_NUMERICAL_AUDIT='1')
        cmd=[a.mpiexec,'-n',str(ranks),str(executable),str(path.resolve()),'-manual_grid',*[str(i) for i in layout]]
        previous={p.name for p in path.glob('????????_?????.??????.h5')} if restart else set()
        t=time.monotonic()
        proc=subprocess.run(cmd,env=env,text=True,capture_output=True,timeout=600)
        elapsed=time.monotonic()-t
        (a.work/(label+'.log')).write_text(proc.stdout+proc.stderr)
        passed=(proc.returncode==0 and checkpoint_shapes(path,ranks,layout)) if expected is None else proc.returncode!=0 and expected in proc.stdout+proc.stderr
        if restart:
            advanced=restart_advanced(path,previous)
            passed=passed and advanced
        row=dict(label=label,returncode=proc.returncode,passed=passed,expected_rejection=expected,wall_seconds=elapsed)
        if restart:row['restart_advanced']=advanced
        runs.append(row);print(json.dumps(row),flush=True)
        return passed
    for case,ranks,layout in cases:
        tag=case+'_'+str(layout[0])+'x'+str(layout[1])
        continuous=create(tag+'_continuous',case,300);split=create(tag+'_split',case,120)
        okay=run(continuous,tag+'_continuous',ranks,layout)
        okay=run(split,tag+'_first',ranks,layout) and okay
        # Preserve a restart seed for the adversarial probes before continuing.
        seed=a.work/(tag+'_seed');shutil.copytree(split,seed)
        set_duration(split,300,extend=True)
        okay=run(split,tag+'_restart',ranks,layout,restart=True) and okay
        if okay:
            comparisons.append(dict(case=tag,**compare(continuous,split,budget)))
            transport.append(dict(case=tag,**check(continuous,layout=layout)))
            continuity.append(dict(case=tag,**check_continuity(continuous)))
            sources.append(dict(case=tag,**check_sources(continuous,ranks)))
            energy.append(dict(case=tag,**check_energy(continuous,ranks)))
        if case=='mini2dns_fang' and ranks==2:
            for fault in ['inputs','layout','missing_state','incomplete_state','bad_dtype','bad_clock','bad_binary',
                          'swap_rank','mismatched_generation','nonfinite_fluid','nonfinite_electro',
                          'nonfinite_vi1','nonfinite_vi2','nonfinite_vi3']:
                bad=a.work/('reject_'+fault);shutil.copytree(seed,bad)
                set_duration(bad,300,extend=True)
                frame=sorted(bad.glob('????????_?????.??????.h5'))[-1]
                with h5py.File(frame,'r+') as f:
                    if fault=='layout':f['restart_runtime/layout'][...]=[2,1]
                    if fault=='bad_binary':f['restart_runtime/executable_sha256'][...]=b'0'*64
                    prefix=f['restart_runtime/prefix'][()].decode().strip()
                state=bad/(prefix+'.r00000000.h5')
                if fault=='inputs':
                    forcing=next((bad/'inputs/Efield').glob('????????_*.h5'))
                    with h5py.File(forcing,'r+') as f:f['Exit'][...]+=0.00001
                if fault=='missing_state':state.unlink()
                if fault=='swap_rank':
                    other=bad/(prefix+'.r00000001.h5')
                    first_bytes,other_bytes=state.read_bytes(),other.read_bytes()
                    state.write_bytes(other_bytes);other.write_bytes(first_bytes)
                if fault in ['incomplete_state','bad_dtype','bad_clock','mismatched_generation'] or fault.startswith('nonfinite_'):
                    with h5py.File(state,'r+') as f:
                        if fault=='incomplete_state':f['complete'][...]=0
                        if fault=='bad_clock':f['dt'][...]=np.nan
                        if fault=='mismatched_generation':f['generation'][...]=b'wrong-generation'
                        if fault.startswith('nonfinite_'):
                            field=fault.removeprefix('nonfinite_')
                            f[field][tuple(n-1 for n in f[field].shape)]=np.nan if field=='fluid' else np.inf
                        if fault=='bad_dtype':
                            data=f['fluid'][...].astype('float32');del f['fluid'];f['fluid']=data
                expected={'inputs':'Changed restart inputs','layout':'Incompatible runtime checkpoint layout',
                          'missing_state':'Missing runtime state file','incomplete_state':'Incomplete runtime checkpoint',
                          'bad_dtype':'requires float64','bad_clock':'Invalid runtime checkpoint clock',
                          'bad_binary':'Changed restart executable','swap_rank':'Runtime state rank mismatch',
                          'mismatched_generation':'Runtime state generation mismatch'}
                expected='Nonfinite runtime checkpoint payload: '+fault.removeprefix('nonfinite_') if fault.startswith('nonfinite_') else expected[fault]
                run(bad,'reject_'+fault,ranks,layout,expected)
            # A C++ run must produce the same physics as Fortran; each self-restart remains binary-bound.
            cpp_exe=(a.build/('gemini_c.bin'+a.suffix)).resolve()
            if cpp_exe.is_file():
                cpp=create('cpp_continuous',case,300)
                cpp_split=create('cpp_split',case,120)
                cpp_okay=run(cpp,'cpp_continuous',ranks,layout,executable=cpp_exe)
                cpp_okay=run(cpp_split,'cpp_first',ranks,layout,executable=cpp_exe) and cpp_okay
                set_duration(cpp_split,300,extend=True)
                cpp_okay=run(cpp_split,'cpp_restart',ranks,layout,executable=cpp_exe,restart=True) and cpp_okay
                if cpp_okay:
                    comparisons.append(dict(case='cpp_vs_fortran',**compare(continuous,cpp,budget)))
                    comparisons.append(dict(case='cpp_restart',**compare(cpp,cpp_split,budget)))
                    energy.append(dict(case='cpp_continuous',**check_energy(cpp,ranks)))
            else:
                runs.append(dict(label='cpp_required',passed=False,error='Required C++ executable is missing'))
    result=dict(schema='gemini.research.matrix.1',executable_sha256=exe_sha,runs=runs,
                 extended_layouts=a.extended_layouts,
                comparisons=comparisons,transport=transport,continuity=continuity,sources=sources,energy=energy,
                passed=all(r['passed'] for r in runs) and len(comparisons)==len(cases)+2 and
                       len(energy)==len(cases)+1 and all(r['passed'] for r in comparisons+transport+continuity+sources+energy))
    (a.work/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['comparisons','transport']},indent=2))
    return 0 if result['passed'] else 1
if __name__=='__main__':raise SystemExit(main())
