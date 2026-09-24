"""Lossless GEMINI full-output interchange; no invented CVTWIN MHD state.

Exports the stored float32 full-output fields unchanged, the grid's native
geomagnetic Cartesian basis, validity mask, and UTC metadata. Does not promote
output precision, interpolate, or provide a two-way physical coupling.
"""
import argparse
from datetime import datetime,timedelta,timezone
import hashlib
import json
import os
from pathlib import Path
import tempfile
import h5py
import numpy as np
from hdf_policy import numeric,dataset
from validate_driver import parse_frame

SPECIES=['O+','NO+','N2+','O2+','N+','H+','e-']
FIELDS={'nsall':('number_density','m^-3'), 'Tsall':('temperature','K'),
        'vs1all':('parallel_velocity','m s^-1'), 'J1all':('current_native_1','A m^-2'),
        'J2all':('current_native_2','A m^-2'), 'J3all':('current_native_3','A m^-2'),
        'v2avgall':('density_averaged_velocity_2','m s^-1'),
        'v3avgall':('density_averaged_velocity_3','m s^-1'), 'Phiall':('potential','V')}
GEOMETRY_UNITS=dict(r='m',alt='m',theta='rad',phi='rad',glat='degree_north',glon='degree_east',
                    basis_geomagnetic_ecef='1',valid='1',x1='m',x2='m',x3='m',signed_background_B1_T='T')


def sha(path):
    with open(path,'rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()


def finite(value,name):
    if not np.isfinite(value).all():raise ValueError('Nonfinite '+name)
    return value


def native_time(date,hour,filename_time):
    """Check the actual 10 ms filename contract, retaining native time separately.

    timeutils:utsec2filestem rounds positive seconds to centiseconds; the six
    displayed fractional digits are not a promise of microsecond precision.
    A rounded filename can cross midnight while its native date has not.
    """
    if date.shape!=(3,) or date.dtype.kind not in 'iu':raise ValueError('Invalid UTC date')
    if hour.shape!=() or hour.dtype.kind!='f' or not np.isfinite(hour):raise ValueError('Invalid UTC hour')
    if not 0<=float(hour)<24:raise ValueError('Invalid UTC hour range')
    try:day=datetime(*(int(v) for v in date),tzinfo=timezone.utc)
    except ValueError as error:raise ValueError('Invalid UTC date') from error
    seconds=float(hour)*3600
    # Division by 3600 at write time and multiplication here can move an
    # exact 5 ms tie by an ulp. Check the rounding interval, not a second
    # independently rounded value; retain the 10 ms filename lattice.
    offset=(filename_time-day).total_seconds()-seconds
    roundoff=4*np.spacing(max(abs(seconds),1.))
    if filename_time.microsecond%10000 or abs(offset)>.005+roundoff:
        raise ValueError('UTC time mismatch with native 10 ms filename rounding')
    return dict(utc=(day+timedelta(seconds=seconds)).isoformat(),
                filename_utc=filename_time.isoformat(),filename_resolution_seconds=.01,
                native_ymd=[int(v) for v in date],native_ut_hour=float(hour),
                native_ut_hour_hex=float(hour).hex(),native_ut_hour_dtype=str(hour.dtype))


def read_native(frame,grid):
    instant=parse_frame(frame)
    if instant is None:raise ValueError('GEMINI timestamp filename required')
    with h5py.File(grid) as g,h5py.File(frame) as f:
        if numeric(f,'flagoutput').shape!=() or int(numeric(f,'flagoutput'))!=1:
            raise ValueError('Only full output flagoutput=1 is supported')
        mask=numeric(g,'nullpts')
        if mask.ndim!=3 or not np.isin(mask,[0,1]).all():raise ValueError('Invalid null mask')
        shape=mask.shape;valid=(mask==0)
        if not valid.any():raise ValueError('Empty physical domain')
        for key in ['h1','h2','h3']:
            metric=numeric(g,key,tuple(n+4 for n in shape))[2:-2,2:-2,2:-2]
            if not np.allclose(metric,1,rtol=0,atol=1e-12):raise ValueError('Only Cartesian unit-metric interchange is supported')
        date=numeric(f,'time/ymd',(3,));hour=numeric(f,'time/UThour')
        clock=native_time(date,hour,instant)
        data={}
        for original,(name,units) in FIELDS.items():
            expected=(7,*shape) if original in ('nsall','Tsall','vs1all') else (shape[:2] if original=='Phiall' else shape)
            data[name]=finite(numeric(f,original,expected),original)
            if data[name].dtype.kind!='f':raise ValueError('Floating-point plasma fields required')
        basis=np.stack([finite(numeric(g,f'e{i}',(3,*shape)),f'e{i}') for i in (1,2,3)],axis=-1)
        basis=np.moveaxis(basis,0,-2)  # x3,x2,x1,Cartesian component,native component
        gram=np.einsum('...ji,...jk->...ik',basis.astype(float),basis.astype(float))
        if np.max(np.abs(gram-np.eye(3)))>2e-6 or np.max(np.abs(np.linalg.det(basis)-1))>2e-6:
            raise ValueError('Basis must be right-handed orthonormal')
        geometry={key:finite(numeric(g,key,shape),key) for key in ['r','theta','phi','glat','glon','alt']}
        geometry['basis_geomagnetic_ecef']=basis
        geometry['valid']=valid
        for key,n in zip(['x1','x2','x3'],shape[::-1]):
            centers=finite(numeric(g,key,(n+4,)),key)[2:-2]
            if n>1 and not np.all(np.diff(centers)>0):raise ValueError('Grid centers must increase')
            geometry[key]=centers
        b=numeric(g,'Bmag',shape);geometry['signed_background_B1_T']=finite(b,'Bmag')
    return data,geometry,dict(**clock,native_frame=frame.name,
        source_sha256=sha(frame),grid_sha256=sha(grid),species=SPECIES,
        axis_order=['x3','x2','x1'],basis='native orthonormal components expressed in geomagnetic Cartesian coordinates',
        interpretation='background_B1 is signed imposed field; no evolved magnetic perturbation is supplied',
        representation='stored native full-output precision; not the float64 restart state')


def export(frame,grid,output):
    data,geometry,meta=read_native(frame,grid)
    output.parent.mkdir(parents=True,exist_ok=True)
    if output.exists():raise FileExistsError(output)
    fd,tmp=tempfile.mkstemp(prefix=output.name+'.',suffix='.partial',dir=output.parent);os.close(fd)
    try:
        with h5py.File(tmp,'w') as f:
            f.attrs['schema']='gemini.exchange.1';f.attrs['metadata']=json.dumps(meta,sort_keys=True)
            for original,(name,units) in FIELDS.items():
                ds=f.create_dataset('state/'+name,data=data[name]);ds.attrs['units']=units;ds.attrs['native_name']=original
            for name,value in geometry.items():
                ds=f.create_dataset('geometry/'+name,data=value);ds.attrs['units']=GEOMETRY_UNITS[name]
        os.replace(tmp,output)
    finally:
        if os.path.exists(tmp):os.unlink(tmp)
    return dict(passed=True,output_sha256=sha(output),**meta)


def read_exchange(path):
    with h5py.File(path) as f:
        if f.attrs.get('schema')!='gemini.exchange.1':raise ValueError('Unsupported exchange schema')
        meta=json.loads(f.attrs['metadata'])
        if meta['species']!=SPECIES or meta['axis_order']!=['x3','x2','x1']:raise ValueError('Unsupported channel/axis order')
        for name,units in GEOMETRY_UNITS.items():
            if dataset(f,'geometry/'+name).attrs.get('units')!=units:raise ValueError('Geometry unit mismatch')
        valid=numeric(f,'geometry/valid')
        if valid.ndim!=3 or not np.isin(valid,[0,1]).all() or not valid.any():raise ValueError('Invalid exchange mask')
        data={}
        for original,(name,units) in FIELDS.items():
            ds=dataset(f,'state/'+name)
            if ds.attrs.get('units')!=units or ds.attrs.get('native_name')!=original:raise ValueError('Unit or quantity mismatch')
            expected=(7,*valid.shape) if original in ('nsall','Tsall','vs1all') else (valid.shape[:2] if original=='Phiall' else valid.shape)
            data[name]=finite(numeric(f,'state/'+name,expected),name)
        basis=finite(numeric(f,'geometry/basis_geomagnetic_ecef',(*valid.shape,3,3)),'basis')
        gram=np.einsum('...ji,...jk->...ik',basis.astype(float),basis.astype(float))
        if np.max(np.abs(gram-np.eye(3)))>2e-6 or np.max(np.abs(np.linalg.det(basis)-1))>2e-6:
            raise ValueError('Invalid exchange basis')
    return data,valid.astype(bool),basis,meta


def restore_output(exchange,output):
    """Restore the supported science fields, not a runnable restart checkpoint."""
    data,valid,basis,meta=read_exchange(exchange)
    if output.exists():raise FileExistsError(output)
    with h5py.File(output,'x') as f:
        f.attrs['purpose']='interchange round trip; not a simulation checkpoint'
        for original,(name,_) in FIELDS.items():f[original]=data[name]
    return output


def cvtwin_mhd6(_exchange):
    raise ValueError('No justified MHD6 mapping: GEMINI lacks evolved magnetic perturbations and species-resolved perpendicular velocities; CVTWIN is nondimensional incompressible MHD. A closure, scales, boundaries and exchange law are required.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--frame',type=Path,required=True)
    p.add_argument('--grid',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    print(json.dumps(export(a.frame,a.grid,a.output),indent=2))
