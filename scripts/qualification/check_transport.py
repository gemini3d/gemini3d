"""Check measured physical-cell transport balances; do not infer unknown sources.

Quantity: 1 species number, 2 parallel momentum (covariant weighting in curved
perpendicular sweeps), 3 internal energy. Null cells and ghosts are excluded.
min_face/max_face include MPI subdomain faces; they are not global boundaries.
"""
import argparse
import csv
import json
from pathlib import Path
import math


def check(directory,budget=1e-11,layout=None):
    files=sorted(directory.glob('transport-r*.csv'))
    if not files:raise ValueError('No native transport evidence')
    count=0;worst=0;coverage=set();worst_row=None;by_rank={}
    for path in files:
        rank=int(path.stem.split("-r")[1])
        if rank in by_rank:raise ValueError('Duplicate transport rank')
        by_rank[rank]={}
        with path.open() as f:
            for row in csv.DictReader(f):
                r={k:float(v) for k,v in row.items()}
                if not all(math.isfinite(v) for v in r.values()):raise ValueError('Nonfinite budget row')
                if r['dt']<=0 or r['scale']<=0:raise ValueError('Nonpositive time/scale')
                if any(r[k]!=int(r[k]) for k in ('quantity','species','axis')):raise ValueError('Noninteger stage identifier')
                scale=max(r['scale'],abs(r['outward']),2.2250738585072014e-308)
                # Recompute using independently recorded state change and actual face transport.
                residual=abs(r['delta']+r['outward'])/scale
                if residual>worst:worst=residual;worst_row={'file':path.name,'t':r['t'],'quantity':r['quantity'],'species':r['species'],'axis':r['axis']}
                key=(r['t'],r['dt'],int(r['quantity']),int(r['species']),int(r['axis']))
                if key in by_rank[rank]:raise ValueError('Duplicate native transport stage')
                by_rank[rank][key]=r
                coverage.add(key[2:])
                count+=1
    expected={(q,s,a) for q in (1,2,3) for s in range(1,8 if q==3 else 7) for a in (1,2,3)}
    complete=True
    for data in by_rank.values():
        clocks={key[:2] for key in data}
        complete=complete and bool(clocks) and all({key[2:] for key in data if key[:2]==clock}==expected for clock in clocks)
    paired=None;global_worst=None;interfaces=0
    if layout is not None:
        l2,l3=layout
        if l2*l3!=len(files) or set(by_rank)!=set(range(l2*l3)):
            raise ValueError('MPI rank inventory differs from declared layout')
        stages=set(by_rank[0])
        if any(set(rows)!=stages for rows in by_rank.values()):raise ValueError('MPI stages/times do not match')
        paired=0.;global_worst=0.
        for key in stages:
            axis=key[-1];rows=[by_rank[rank][key] for rank in range(l2*l3)]
            scale=max(sum(r['scale'] for r in rows),2.2250738585072014e-308)
            outer=0.;cutout=sum(r['outward']-r['min_face']-r['max_face'] for r in rows)
            for rank,r in enumerate(rows):
                i2=rank%l2;i3=rank//l2
                i,extent=(0,1) if axis==1 else ((i2,l2) if axis==2 else (i3,l3))
                if i==0:outer+=r['min_face']
                if i==extent-1:outer+=r['max_face']
                if i<extent-1:
                    neighbor=rank+(1 if axis==2 else l2)
                    paired=max(paired,abs(r['max_face']+rows[neighbor]['min_face'])/scale);interfaces+=1
            global_worst=max(global_worst,abs(sum(r['delta'] for r in rows)+outer+cutout)/scale)
    return dict(schema='gemini.native.transport.1',passed=count>0 and complete and coverage==expected and worst<=budget and
        (paired is None or max(paired,global_worst)<=budget),
        layout=layout,interface_pairs=interfaces,max_interface_mismatch=paired,max_global_residual=global_worst,
        rank_files=len(files),rows=count,per_step_coverage_complete=complete,max_normalized_residual=worst,budget=budget,worst=worst_row,
        coverage=[list(x) for x in sorted(coverage)],missing_coverage=[list(x) for x in sorted(expected-coverage)],
        scope='Native split-advection balance with physical/null boundaries; MPI interfaces paired when layout is supplied. Full coupled source and energy conservation is separate.')
def check_continuity(directory,budget=1e-11,charge_budget=1e-12):
    files=sorted(directory.glob('continuity-r*.csv'))
    if not files:raise ValueError('No native full-step continuity evidence')
    rows=0;worst=0;charge=0;cleanup={str(s):0. for s in range(1,7)};coverage=set();by_rank={}
    for path in files:
        rank=int(path.stem.split('-r')[1])
        if rank in by_rank:raise ValueError('Duplicate continuity rank')
        by_rank[rank]={}
        with path.open() as f:
            for row in csv.DictReader(f):
                r={k:float(v) for k,v in row.items()}
                if not all(math.isfinite(v) for v in r.values()):raise ValueError('Nonfinite continuity evidence')
                if r['scale']<=0:raise ValueError('Nonpositive continuity scale')
                if r['dt']<=0 or r['charge_fraction']<0:raise ValueError('Invalid continuity time/charge')
                species=int(r['species'])
                if species!=r['species'] or species not in range(1,7):raise ValueError('Invalid ion species')
                key=(r['t'],r['dt']);present=by_rank[rank].setdefault(key,set())
                if species in present:raise ValueError('Duplicate continuity stage')
                present.add(species)
                residual=abs(r['after']-r['before']-r['source']+r['outward']-r['cleanup'])/r['scale']
                worst=max(worst,residual);charge=max(charge,r['charge_fraction'])
                coverage.add(species);cleanup[str(species)]+=r['cleanup'];rows+=1
    complete=set(by_rank)==set(range(len(files))) and all(bool(clocks) and all(s==set(range(1,7)) for s in clocks.values()) for clocks in by_rank.values())
    if complete:complete=all(set(clocks)==set(by_rank[0]) for clocks in by_rank.values())
    return dict(schema='gemini.native.continuity.1',passed=rows>0 and complete and coverage==set(range(1,7)) and worst<=budget and charge<=charge_budget,
        rank_files=len(files),rows=rows,per_step_coverage_complete=complete,max_normalized_residual=worst,max_charge_fraction=charge,
        budget=budget,charge_budget=charge_budget,total_cleanup_number_by_species=cleanup,
        scope='Six native ion continuity equations with applied chemical/ionization source increments and explicit cleanup; electron quasi-neutrality. Reaction-rate accuracy, full momentum and energy balances are separate.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('directory',type=Path);p.add_argument('--output',type=Path,required=True);p.add_argument('--layout',type=int,nargs=2);a=p.parse_args()
    try:r=check(a.directory,layout=a.layout)
    except (ValueError,OSError) as e:r=dict(passed=False,error=str(e))
    a.output.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2));raise SystemExit(0 if r['passed'] else 1)
