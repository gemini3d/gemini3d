"""Measure native cross-midnight restart error; never rewrite reference data.

Needs numpy and h5py. Reuses a downloaded mini2dns_fang input case. All runs
are isolated below --work. Exit 0 requires successful runs and invalid restart
rejection; reported numerical errors still require a science acceptance bound.
"""
import argparse
from pathlib import Path
import datetime as dt
import json
import os
import re
import shutil
import subprocess
import h5py
import numpy as np


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--build',type=Path,required=True)
    p.add_argument('--input-case',type=Path,required=True)
    p.add_argument('--work',type=Path,required=True)
    p.add_argument('--mpiexec',default='mpiexec')
    args=p.parse_args()
    work=args.work.resolve();work.mkdir(parents=True,exist_ok=True)
    origin=dt.datetime(2013,2,20,5)
    new_origin=dt.datetime(2013,2,20,23,59,30)
    def create(name,duration):
        case=work/name
        if case.exists(): raise RuntimeError('Use a new work directory: '+str(case))
        shutil.copytree(args.input_case/'inputs',case/'inputs')
        for f in sorted((case/'inputs').rglob('*.h5')):
            if not re.fullmatch(r'\d{8}_\d{5}\.\d{6}\.h5',f.name):continue
            d,s=f.stem.split('_'); instant=dt.datetime.strptime(d,'%Y%m%d')+dt.timedelta(seconds=float(s))
            shifted=new_origin+(instant-origin)
            seconds=(shifted-shifted.replace(hour=0,minute=0,second=0,microsecond=0)).total_seconds()
            target=f.with_name(shifted.strftime('%Y%m%d_')+f'{seconds:012.6f}.h5')
            f.rename(target)
            with h5py.File(target,'r+') as h:
                if 'time/ymd' in h:h['time/ymd'][...]=[shifted.year,shifted.month,shifted.day]
                if 'time/UTsec' in h:h['time/UTsec'][...]=seconds
                if 'time/UThour' in h:h['time/UThour'][...]=seconds/3600
        cfg=case/'inputs/config.nml'
        text=cfg.read_text()
        text=re.sub(r'(?im)^(UTsec0\s*=)[^!\n]*',r'\g<1> 86370.0 ',text)
        text=re.sub(r'(?im)^(tdur\s*=)[^!\n]*',r'\g<1> '+str(duration)+' ',text)
        cfg.write_text(text)
        return case
    runs=[]
    def run(case,label,expected=0):
        proc=subprocess.run([args.mpiexec,'-n','2',str(args.build.resolve()/'gemini.bin'),str(case)],
                            capture_output=True,text=True,timeout=180)
        (work/(label+'.log')).write_text(proc.stdout+proc.stderr)
        runs.append(dict(label=label,returncode=proc.returncode,expected_success=expected==0,
                         passed=(proc.returncode==0)==(expected==0)))
        return proc
    continuous=create('continuous',300)
    split=create('split',120)
    run(continuous,'continuous');run(split,'split_first')
    cfg=split/'inputs/config.nml';cfg.write_text(re.sub(r'(?im)^(tdur\s*=)[^!\n]*',r'\g<1> 300 ',cfg.read_text()))
    run(split,'split_restart')
    end='20130221_00270.000000.h5'
    metrics=[]
    if (continuous/end).exists() and (split/end).exists():
        with h5py.File(continuous/end) as a,h5py.File(split/end) as b:
            for key in ['nsall','Tsall','vs1all','Phiall','J1all','J2all','J3all','v2avgall','v3avgall']:
                av=a[key][...].astype(float);bv=b[key][...].astype(float)
                metrics.append(dict(field=key,finite=bool(np.isfinite(bv).all()),
                                    relative_l2=float(np.linalg.norm(bv-av)/max(np.linalg.norm(av),1e-30)),
                                    max_absolute=float(np.max(np.abs(bv-av))),bitwise_equal=bool(np.array_equal(av,bv))))
    for label in ['incomplete','timestamp','field_resolved']:
        case=work/label;shutil.copytree(split,case)
        cfg=case/'inputs/config.nml'
        cfg.write_text(re.sub(r'(?im)^(tdur\s*=)[^!\n]*',r'\g<1> 360 ',cfg.read_text()))
        frame=case/end
        if label=='incomplete':
            with h5py.File(frame,'r+') as f:del f['Tsall']
        elif label=='timestamp':
            with h5py.File(frame,'r+') as f:f['time/UThour'][...]=10
        else:
            cfg.write_text(cfg.read_text().replace('potsolve = 1','potsolve = 3'))
            # Historical slab-only records remain ineligible for field-resolved restart.
            for historical in case.glob('*.h5'):
                with h5py.File(historical,'r+') as f:
                    if 'restart_core' in f: del f['restart_core']
        proc=run(case,label,1)
        expected={'incomplete':'Incomplete restart frame','timestamp':'Restart timestamp mismatch',
                  'field_resolved':'Field-resolved restart requires'}[label]
        runs[-1]['passed'] &= expected in (proc.stdout+proc.stderr)
    payload=dict(schema='gemini.audit.restart.1',runs=runs,metrics=metrics,
                 software_checks_passed=all(x['passed'] for x in runs) and bool(metrics),
                 scientific_equivalence='not certified: choose application-specific bounds and complete checkpoint state')
    (work/'results.json').write_text(json.dumps(payload,indent=2)+'\n')
    print(json.dumps(payload,indent=2))
    return 0 if payload['software_checks_passed'] else 1

if __name__=='__main__':raise SystemExit(main())
