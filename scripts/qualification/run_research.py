"""Launch the admitted exact-restart research profile with immutable-input hashes.

Only tdur may change when extending a run. All files below inputs, including
forcing and grid data, bind to the checkpoint. This does not authenticate users.
"""
import argparse
import hashlib
import datetime as dt
from validate_driver import validate,utc
import json
import os
from pathlib import Path
import re
import subprocess


    root=case/'inputs'
    if root.is_symlink():raise ValueError('Input symlinks are unsupported in the verified profile')
    root=root.resolve()
    records={}
        if not path.is_file():continue
        data=path.read_bytes()
        if path.name=='config.nml':
            text=data.decode()
            if len(re.findall(r'(?im)^\s*tdur\s*=',text))!=1:raise ValueError('Exactly one tdur entry required')
            text=re.sub(r'(?im)^(\s*tdur\s*=)[^!\n]*',r'\1 <EXTENDABLE> ',text)
            data=text.encode()
        records[str(path.relative_to(root))]=hashlib.sha256(data).hexdigest()
    if 'config.nml' not in records:raise ValueError('Missing inputs/config.nml')
    digest=hashlib.sha256(json.dumps(records,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    return {'schema':'gemini.research.inputs.1','sha256':digest,'files':records,'mutable':['base.tdur']}


def preflight(case):
    # This launcher deliberately accepts the simple one-assignment-per-line
    # namelist convention of the pinned research inputs, not arbitrary Fortran.
    text=(case/'inputs/config.nml').read_text()
    values={}
    for line in text.splitlines():
        line=line.split('!',1)[0].strip()
        if not line or line.startswith(('&','/')):continue
        if '=' not in line:raise ValueError('Research namelist requires one complete assignment per line')
        key,value=line.split('=',1);key=key.strip().lower()
        if key in values:raise ValueError('Duplicate research configuration key: '+key)
        values[key]=value.strip().rstrip(',')
    date=[int(x.strip()) for x in values['ymd'].split(',')]
    start=dt.datetime(*date,tzinfo=dt.timezone.utc)+dt.timedelta(seconds=float(values['utsec0']))
    stop=start+dt.timedelta(seconds=float(values['tdur']))
    if date[0]>2025:raise ValueError('Research profile does not qualify production geometry beyond 2025')
    if values.get('allow_missing_spatial','.false.').lower() not in ('.false.','f'):raise ValueError('Research profile requires strict spatial coverage')
    for key in ['indat_size','indat_grid','indat_file']:
        target=(case/values[key].strip('"\'')).resolve()
        if not target.is_relative_to((case/'inputs').resolve()) or not target.is_file():
            raise ValueError('Initial state/grid paths must remain inside verified inputs')
    result={}
    for kind,folder_key,cadence_key in [('efield','e0_dir','dte0'),('precip','prec_dir','dtprec')]:
        folder=values[folder_key].strip('\"\'')
        source=(case/folder).resolve()
        if not source.is_relative_to((case/'inputs').resolve()):raise ValueError('Driver paths must remain inside verified inputs')
        result[kind]=validate(source,start,stop,float(values[cadence_key]),kind)
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--exe',required=True,type=Path);p.add_argument('--case',required=True,type=Path)
    p.add_argument('--ranks',type=int,choices=[1,2],default=2);p.add_argument('--mpiexec',default='mpiexec')
    p.add_argument('--layout',type=int,nargs=2);a=p.parse_args()
    forcing=preflight(a.case)
    (a.case/'forcing-preflight.json').write_text(json.dumps(forcing,indent=2)+'\n')
    record=input_signature(a.case)
    (a.case/'verified-inputs.json').write_text(json.dumps(record,indent=2)+'\n')
    env=os.environ.copy();env.update(HDF5_PLUGIN_PRELOAD='::',GEMINI_EXACT_RESTART='1',GEMINI_INPUT_SHA256=record['sha256'],
        GEMINI_EXECUTABLE_SHA256=hashlib.sha256(a.exe.read_bytes()).hexdigest())
    command=[a.mpiexec,'-n',str(a.ranks),str(a.exe.resolve()),str(a.case.resolve())]
    if a.layout:command+=['-manual_grid',*[str(x) for x in a.layout]]
    return subprocess.run(command,env=env).returncode
if __name__=='__main__':raise SystemExit(main())
