"""Screen the existing vertical Cartesian field approximation against IGRF-14.

Full-grid diagnostic, not a new production geometry. The old Cartesian model
imposes -50 uT along e1; varying only its magnitude would leave the field-line
direction and transport metric assumptions unresolved.
"""
import argparse
import json
from pathlib import Path
import h5py
import numpy as np
from hdf_policy import numeric
from igrf14 import coefficients,field
from state_exchange import sha


def compare(grid,coeff,years,magnitude_budget,direction_budget):
    if not years or not (0<magnitude_budget<1 and 0<direction_budget<90):raise ValueError('Predeclared geometry budgets required')
    with h5py.File(grid) as f:
        r,lat,lon,b,mask=[numeric(f,key) for key in ('r','glat','glon','Bmag','nullpts')]
        e1,er=numeric(f,'e1'),numeric(f,'er')
        if any(x.shape!=r.shape for x in (lat,lon,b,mask)) or not all(np.isfinite(x).all() for x in (r,lat,lon,b)):
            raise ValueError('Invalid geometry arrays')
        if not np.allclose(e1,er,atol=2e-6,rtol=0) or not np.allclose(b,-50e-6,atol=1e-10,rtol=0):
            raise ValueError('This screen supports only the current vertical Cartesian -50 uT model')
        if not np.isin(mask,[0,1]).all():raise ValueError('Invalid geometry mask')
    valid=mask==0
    if not valid.any():raise ValueError('No physical grid points')
    rows=coefficients(coeff);results=[]
    for year in years:
        magnitude=[];direction=[]
        for radius,latitude,longitude,old in zip(r[valid],lat[valid],lon[valid],b[valid]):
            north,east,down=field(rows,year,float(radius)/1000,float(latitude),float(longitude))
            norm=np.linalg.norm([north,east,down])
            magnitude.append(float(abs(abs(float(old))*1e9-norm)/norm))
            direction.append(float(np.degrees(np.arccos(np.clip(down/norm,-1,1)))))
        results.append(dict(year=year,points=len(magnitude),max_relative_magnitude_error=max(magnitude),
                            max_direction_error_degrees=max(direction),
                            passed=bool(max(magnitude)<=magnitude_budget and max(direction)<=direction_budget)))
    return dict(schema='gemini.geometry.screen.1',passed=all(x['passed'] for x in results),
                grid_sha256=sha(grid),coefficients_sha256=sha(coeff),years=results,
                magnitude_budget=magnitude_budget,direction_budget_degrees=direction_budget,
                scope='Full physical grid screen of existing vertical Cartesian field against full IGRF-14; production guard remains')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--grid',type=Path,required=True)
    p.add_argument('--coefficients',type=Path,required=True);p.add_argument('--years',type=float,nargs='+',required=True)
    p.add_argument('--magnitude-budget',type=float,required=True);p.add_argument('--direction-budget-degrees',type=float,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    result=compare(a.grid,a.coefficients,a.years,a.magnitude_budget,a.direction_budget_degrees)
    a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
    raise SystemExit(0 if result['passed'] else 2)
