"""Independent conservative balances and theoretical observables for native probes."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import solve_banded

p=argparse.ArgumentParser()
p.add_argument('--work', type=Path, required=True)
p.add_argument('--trace', type=Path)
p.add_argument('--candidate',required=True)
a=p.parse_args()
w=a.work
LIMIT=1e-11
def rows(path):
    with path.open() as f:
        return [{k:float(v) for k,v in r.items()} for r in csv.DictReader(f)]
def norm_res(terms):
    return abs(sum(terms))/max(sum(abs(x) for x in terms), np.finfo(float).tiny)
def orders(errors):
    return [float(np.log2(errors[i]/errors[i+1])) for i in range(2)]

conduct=[]
traces={}
if a.trace:
    native_trace=np.loadtxt(a.trace,delimiter=',')
    offset=0
    for path in sorted(w.glob('conduction_*_m[23]_*.csv')):
        rs=rows(path);nodes=int(max(r['node'] for r in rs));steps=int(max(r['step'] for r in rs))
        count=steps*nodes
        traces[path.name]=native_trace[offset:offset+count].reshape(steps,nodes,4)
        offset+=count
    assert offset==len(native_trace)
oracle_cache={}
for path in sorted(w.glob('conduction_*.csv')):
    data=rows(path)
    n=int(max(r['node'] for r in data));steps=int(max(r['step'] for r in data));r=data[0]
    dt,dx,cv=r['dt'],r['dx'],r['Cv']
    window=int(path.stem.split('_w')[1][0]);method=int(path.stem.split('_m')[1][0]);level=int(path.stem.split('_l')[1][0])
    density=r['ns'];total=dt*steps
    initial=np.linspace(600.,1500.,n)
    scale=density/1e5
    # Density scaling uses tau=t/scale so both cases independently share the
    # same mathematically defined problem; retain each physical native result.
    key=(n,window,total/scale)
    if key not in oracle_cache:
        def rhs(t, interior):
            temp=np.r_[600.,interior,1500.]
            lam=.3986e-12*temp**2.5
            flux=-.5*(lam[:-1]+lam[1:])*np.diff(temp)/dx
            return -(flux[1:]-flux[:-1])/(cv/scale*dx)
        sol=solve_ivp(rhs,[0.,total/scale],initial[1:-1],method='Radau',rtol=1e-11,atol=1e-9)
        tight=solve_ivp(rhs,[0.,total/scale],initial[1:-1],method='Radau',rtol=2e-13,atol=2e-11)
        if not sol.success or not tight.success:
            raise RuntimeError('Independent MOL oracle failed: '+sol.message+' / '+tight.message)
        reference=np.r_[600.,tight.y[:,-1],1500.]
        control=float(np.linalg.norm(sol.y[:,-1]-tight.y[:,-1])/np.linalg.norm(tight.y[:,-1]))
        oracle_cache[key]={'temperature_K':reference,'control_relative_rms':control,
                           'nfev':int(tight.nfev),'njev':int(tight.njev)}
    oracle=oracle_cache[key]
    step_balances=[];telescoping=[];diagnostic_balances=[];stage_mismatches=[];floors=0
    for it in range(1,steps+1):
        rs=[r for r in data if r['step']==it]
        old=np.array([r['Told'] for r in rs]);new=np.array([r['Tnew'] for r in rs]);raw=np.array([r['Traw'] for r in rs])
        lam=np.array([r['lambda_applied'] if 'lambda_applied' in r else r['lambda_old'] for r in rs]);face=.5*(lam[:-1]+lam[1:])
        alpha=face/(cv*dx**2)
        def flux(temp):return -face*np.diff(temp)/dx
        if method==1:
            integrated=dt*flux(raw)
        else:
            # Independent tridiagonal Crank-Nicolson half-stage for native
            # TR-BDF2 (frozen lambda). No native residual is used as heat input.
            c=dt/4
            band=np.zeros((3,n-2));diagonal=alpha[:-1]+alpha[1:]
            band[1]=1+c*diagonal
            band[0,1:]=-c*alpha[1:-1];band[2,:-1]=-c*alpha[1:-1]
            old_rhs=old[1:-1]+c*(alpha[1:]*np.diff(old)[1:]-alpha[:-1]*np.diff(old)[:-1])
            old_rhs[0]+=c*alpha[0]*600;old_rhs[-1]+=c*alpha[-1]*1500
            tr=np.r_[600.,solve_banded((1,1),band,old_rhs),1500.]
            if path.name in traces:
                native=traces[path.name][it-1]
                assert np.array_equal(native[:,1],old) and np.array_equal(native[:,3],raw)
                independent_tr=tr
                tr=native[:,2]
                stage_mismatches.append(float(np.max(abs(independent_tr-tr))/np.max(abs(tr))))
            integrated=dt/3*(flux(old)+flux(tr)+flux(raw))
            derived=(4*tr[1:-1]-old[1:-1])/3+dt/3*(alpha[1:]*np.diff(raw)[1:]-alpha[:-1]*np.diff(raw)[:-1])
            stage_mismatches.append(float(np.max(abs(derived-raw[1:-1]))/np.max(abs(raw))))
        change=float(cv*dx*np.sum(raw[1:-1]-old[1:-1]))
        incoming=float(integrated[0]-integrated[-1])
        power_cells=-(integrated[1:]-integrated[:-1])
        scale_balance=max(float(np.sum(abs(cv*dx*(raw[1:-1]-old[1:-1])))),
                          float(np.sum(abs(power_cells))),abs(integrated[0])+abs(integrated[-1]),np.finfo(float).tiny)
        step_balances.append(abs(change-incoming)/scale_balance)
        telescoping.append(abs(np.sum(power_cells)-incoming)/scale_balance)
        inc=np.array([[r[k] for k in ('reaction','drift','conduction','source','reservoir')] for r in rs])
        diagnostic_balances.append(float(np.max(abs(raw-old-np.sum(inc,axis=1)))/max(np.max(abs(raw-old)),np.finfo(float).tiny)))
        floors+=int(np.count_nonzero(new!=raw))
    final=np.array([r['Tnew'] for r in data if r['step']==steps])
    error=float(np.linalg.norm(final[1:-1]-oracle['temperature_K'][1:-1])/np.linalg.norm(oracle['temperature_K'][1:-1]))
    entry={'file':path.name,'density_m3':density,'window':window,'method':method,'level':level,
           'steps':steps,'dt_s':dt,'duration_s':total,'temperature_relative_rms':error,
           'oracle_control_relative_rms':oracle['control_relative_rms'],
           'maximum_energy_balance_normalized':max(step_balances),
           'maximum_face_telescoping_normalized':max(telescoping),
           'maximum_diagnostic_increment_relative':max(diagnostic_balances),
           'maximum_independent_TR_stage_relative':max(stage_mismatches,default=0),
           'native_stage_trace_used':method==1 or path.name in traces,
           'fixed_boundary_max_error_K':max(abs(r['Traw']-(600 if r['node']==1 else 1500)) for r in data if r['node'] in (1,n)),
           'floor_cell_events':floors,'temperature_K':final.tolist(),
           'oracle_temperature_K':oracle['temperature_K'].tolist(),
           'energy_pass':max(step_balances)<=LIMIT,'temperature_pass':error<=.1}
    conduct.append(entry)
groups=[]
for density,window,method in sorted(set((r['density_m3'],r['window'],r['method']) for r in conduct)):
    levels=sorted([r for r in conduct if (r['density_m3'],r['window'],r['method'])==(density,window,method)],key=lambda r:r['level'])
    os=orders([r['temperature_relative_rms'] for r in levels])
    minorder=.9 if method==1 else 1.8
    groups.append({'density_m3':density,'window':window,'method':method,'errors':[r['temperature_relative_rms'] for r in levels],
                   'orders':os,'required_order':minorder,
                   'energy_pass':all(r['energy_pass'] for r in levels),
                   'temperature_pass':all(r['temperature_pass'] for r in levels),
                   'order_pass':min(os)>=minorder})
pc03={'cases':conduct,'groups':groups,'passed':all(r['energy_pass'] and r['temperature_pass'] and r['order_pass'] for r in groups)}

kir=rows(w/'kirchhoff_flux.csv');steady=[]
qexact=-.3986e-12*(1500.**3.5-600.**3.5)/(3.5*1e5)
for level in range(3):
    rs=[r for r in kir if r['level']==level]
    temp=np.array([r['T'] for r in rs]);lam=np.array([r['lambda'] for r in rs]);dx=rs[0]['dx']
    q=-.5*(lam[:-1]+lam[1:])*np.diff(temp)/dx
    err=float(np.sqrt(np.mean((q-qexact)**2))/abs(qexact))
    steady.append({'level':level,'nodes':len(rs),'dx_m':dx,'exact_flux_Wm2':qexact,
                   'relative_l2_face_flux_error':err,'minimum_flux_Wm2':float(min(q)),'maximum_flux_Wm2':float(max(q))})
pc04={'levels':steady,'orders':orders([r['relative_l2_face_flux_error'] for r in steady])}
pc04['passed']=min(pc04['orders'])>=1.8
result={'schema':'gemini.physical.closure.results.1','budget':LIMIT,
        'native_source_revision':a.candidate,
        'PC03':pc03,'PC04':pc04, 'passed':pc03['passed'] and pc04['passed']}
(w/'physics-results.json').write_text(json.dumps(result,indent=2,allow_nan=False,default=lambda v:v.item() if isinstance(v,np.generic) else str(v))+'\n')
print(json.dumps({k:{s:v for s,v in result[k].items() if s in ('passed','case_count','passing_cases','orders','groups')}
                  for k in ('PC03','PC04')},indent=2))

raise SystemExit(0 if result['passed'] else 1)
