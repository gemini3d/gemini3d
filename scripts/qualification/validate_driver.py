"""Preflight a uniform-cadence HDF5 forcing stream; gaps are errors.

Run separately for each driver. It does not infer a cadence, fill missing
values, or resample nonuniform data. Times are UTC (naive ISO values mean UTC).
"""
import argparse
import hashlib
import datetime as dt
import json
from pathlib import Path
import re
import h5py
import numpy as np
from hdf_policy import inspect_dataset


def parse_frame(path):
    match=re.fullmatch(r'(\d{8})_(\d{5}\.\d{6})\.h5',path.name)
    if not match:return None
    sec=float(match[2])
    if not 0<=sec<86400:raise ValueError('Invalid UTC seconds in '+path.name)
    return dt.datetime.strptime(match[1],'%Y%m%d').replace(tzinfo=dt.timezone.utc)+dt.timedelta(seconds=sec)


def numeric_dataset(group, name):
    """Reject indirection before dereferencing; forcing is a self-contained file."""
    if name not in group or not isinstance(group.get(name, getlink=True), h5py.HardLink):
        raise ValueError('Missing required field or unsupported link: '+name)
    obj=group[name]
    if not isinstance(obj,h5py.Dataset) or obj.dtype.kind not in 'iuf':
        raise ValueError('Expected real numeric forcing dataset: '+name)
    inspect_dataset(obj)
    value=obj[...]
    if not np.isfinite(value).all():raise ValueError('Nonfinite forcing dataset: '+name)
    return value


def driver_grid(directory):
    sizes=[]
    with h5py.File(directory/'simsize.h5') as f:
        for aliases in [('llat','Nlat','lx2'),('llon','Nlon','lx3')]:
            key=next((k for k in aliases if k in f),aliases[0])
            value=numeric_dataset(f,key)
            if value.shape!=() or value<=0 or value!=int(value):
                raise ValueError('Driver size must be a positive integer scalar: '+key)
            sizes.append(int(value))
    with h5py.File(directory/'simgrid.h5') as f:
        for name,length in zip(['mlat','mlon'],sizes):
            value=numeric_dataset(f,name)
            if value.shape!=(length,) or (length>1 and not np.all(np.diff(value)>0)):
                raise ValueError('Driver coordinate shape or increasing order: '+name)
            if name=='mlat' and np.any(np.abs(value)>90):
                raise ValueError('Driver latitude outside degrees [-90,90]')
    return tuple(sizes)


def required_schema(f,kind,shape):
    lat,lon=shape
    fields=({'Qp':shape,'E0p':shape} if kind=='precip' else
            dict(Exit=shape,Eyit=shape,Vminx1it=shape,Vmaxx1it=shape,
                 Vminx2ist=(lat,),Vmaxx2ist=(lat,),Vminx3ist=(lon,),Vmaxx3ist=(lon,),flagdirich=()))
    for name,expected in fields.items():
        value=numeric_dataset(f,name)
        if value.shape!=expected:
            raise ValueError(f'Required field shape {name}: {value.shape}, expected {expected} (disk lat,lon order)')
        if name=='flagdirich' and float(value) not in (0.,1.,2.):
            raise ValueError('Invalid categorical flagdirich; expected 0, 1 or 2')
        if name=='Qp' and np.any(value<0):raise ValueError('Negative precipitation energy flux Qp')
        if name=='E0p' and np.any(value<=0):raise ValueError('Precipitation characteristic energy E0p must be positive')
    if 'time' not in f or not isinstance(f.get('time',getlink=True),h5py.HardLink) or not isinstance(f['time'],h5py.Group):
        raise ValueError('Required driver time group missing or linked')
    numeric_dataset(f['time'],'ymd')
    if not any(k in f['time'] for k in ['UTsec','UThour']):raise ValueError('Required UTC time metadata missing')


def validate(directory,start,stop,cadence,kind='generic'):
    if kind not in ('generic','efield','precip'):raise ValueError('Unsupported driver kind')
    shape=driver_grid(directory) if kind!='generic' else None
    if not np.isfinite(cadence) or cadence<=0 or stop<start:raise ValueError('Invalid requested time span or cadence')
    frames=sorted((instant,p) for p in directory.glob('*.h5') if (instant:=parse_frame(p)) is not None)
    if len(frames)<2:raise ValueError('At least two bracketing driver frames are required')
    if frames[0][0]>start or frames[-1][0]<stop:raise ValueError('Driver does not cover requested interval')
    for (ta,pa),(tb,pb) in zip(frames,frames[1:]):
        if abs((tb-ta).total_seconds()-cadence)>1e-6:
            raise ValueError(f'Missing or nonuniform cadence: {pa.name} -> {pb.name}')
    signature=None;count=0
    for instant,path in frames:
        fields={}
        with h5py.File(path) as f:
            visited={h5py.h5o.get_info(f.id).addr}
            def walk(group,prefix=''):
                nonlocal count
                for key in group:
                    name=prefix+'/'+key
                    if not isinstance(group.get(key,getlink=True),h5py.HardLink):
                        raise ValueError('External/soft driver links are unsupported: '+name)
                    obj=group[key]
                    if isinstance(obj,h5py.Group):
                        address=h5py.h5o.get_info(obj.id).addr
                        if address in visited:raise ValueError('Repeated/cyclic forcing group: '+name)
                        visited.add(address);walk(obj,name);continue
                    if obj.dtype.kind not in 'iuf':raise ValueError('Nonnumeric or non-real forcing dataset: '+name)
                    inspect_dataset(obj)
                    value=obj[...]
                    if not np.isfinite(value).all():raise ValueError('Nonfinite forcing dataset: '+str(path)+name)
                    fields[name]=list(obj.shape);count+=1
            walk(f)
            if kind!='generic':required_schema(f,kind,shape)
            if 'time/ymd' in f:
                date=np.asarray(f['time/ymd'][...])
                if date.shape!=(3,) or not np.array_equal(date,[instant.year,instant.month,instant.day]):
                    raise ValueError('Driver date metadata differs from filename: '+str(path))
            for key,multiplier in [('time/UTsec',1),('time/UThour',3600)]:
                if key in f:
                    stamp=np.asarray(f[key][...])
                    expected=(instant-instant.replace(hour=0,minute=0,second=0,microsecond=0)).total_seconds()
                    if stamp.size!=1 or not np.isfinite(stamp).all() or abs(float(stamp.item())*multiplier-expected)>1e-6:
                        raise ValueError('Driver UTC metadata differs from filename: '+str(path))
        if not fields:raise ValueError('Empty forcing frame: '+str(path))
        if signature is None:signature=fields
        elif fields!=signature:raise ValueError('Dataset/shape changes between driver frames: '+str(path))
    hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for _,p in frames}
    if shape is not None:
        for name in ['simsize.h5','simgrid.h5']:hashes[name]=hashlib.sha256((directory/name).read_bytes()).hexdigest()
    return dict(passed=True,kind=kind,required_schema_checked=kind!='generic',sha256=hashes,
                frames=len(frames),datasets_checked=count,cadence_seconds=cadence,
                first=frames[0][0].isoformat(),last=frames[-1][0].isoformat(),fields=signature)


def utc(value):
    value=dt.datetime.fromisoformat(value)
    return value.replace(tzinfo=dt.timezone.utc) if value.tzinfo is None else value.astimezone(dt.timezone.utc)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('directory',type=Path);p.add_argument('--start',required=True,type=utc)
    p.add_argument('--stop',required=True,type=utc);p.add_argument('--cadence',type=float,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--kind',choices=['generic','efield','precip'],default='generic',help='A concrete kind is required for qualification')
    a=p.parse_args()
    try:result=validate(a.directory,a.start,a.stop,a.cadence,a.kind)
    except (ValueError,OSError) as e:result=dict(passed=False,error=str(e))
    result['schema']='gemini.qualification.driver.1';a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
    return 0 if result['passed'] else 1
if __name__=='__main__':raise SystemExit(main())
