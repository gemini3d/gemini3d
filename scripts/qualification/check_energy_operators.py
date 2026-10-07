"""Check every native compression and parabolic-energy stage and species.

The boundary channel is the change in algebraically prescribed endpoint cells.
It is NOT mislabeled as conductive heat flux. The interior conduction channel
uses flux divergences evaluated at the actual implicit solver stages. Whole-step
coupled momentum/energy and physical closure validation remain separate gates.
Four ranks support the optional 2x2 regression, not an expanded science profile.
"""
import argparse
import json
import math
from pathlib import Path
from check_sources import rows

EXPECTED={(stage,s) for stage in (1,2) for s in range(1,8)}


def check(directory,ranks,budget=1e-11):
    directory=Path(directory)
    if ranks not in (1,2,4) or not math.isfinite(budget) or budget<=0:
        raise ValueError('Admitted rank count and positive finite budget required')
    expected={f'energy-operators-r{r:08d}.csv' for r in range(ranks)}
    if {p.name for p in directory.glob('energy-operators-r*.csv')}!=expected:
        raise ValueError('Energy operator rank inventory mismatch')
    clocks=None;count=0;worst=0.;totals={};row_residual=0.
    for rank in range(ranks):
        stages={}
        for r in rows(directory/f'energy-operators-r{rank:08d}.csv'):
            key=(int(r['stage']),int(r['species']));clock=(r['t'],r['dt'])
            if key not in EXPECTED or key in stages.setdefault(clock,set()):
                raise ValueError('Unknown or duplicate energy operator stage')
            stages[clock].add(key)
            terms=[r[f'term{i}_J'] for i in range(1,6)]
            scale=r['scale_J'];delta=r['delta_J']
            if scale<=0 or scale<(1-1e-12)*max(abs(r['before_J']),abs(r['after_J']),sum(map(abs,terms))):
                raise ValueError('Invalid energy scale')
            if key[0]==1 and any(terms[i]!=0 for i in (2,3,4)):
                raise ValueError('Unexpected compression channels')
            residual=delta-math.fsum(terms)
            worst=max(worst,abs(residual)/scale,abs(r['after_J']-r['before_J']-delta)/scale)
            row_residual=max(row_residual,abs(r['residual_J']-residual)/scale)
            for i,value in enumerate(terms,1):
                label=f'{key[0]}:{key[1]}:{i}';totals[label]=totals.get(label,0.)+value
            count+=1
        if not stages or any(v!=EXPECTED for v in stages.values()):
            raise ValueError('Missing energy stage/species at a time step')
        these=set(stages)
        if clocks is not None and these!=clocks:raise ValueError('MPI energy clocks differ')
        clocks=these
        continuity={(r['t'],r['dt']) for r in rows(directory/f'continuity-r{rank:08d}.csv')}
        if these!=continuity:raise ValueError('Energy/continuity step inventory differs')
    return dict(schema='gemini.energy.operators.1',passed=max(worst,row_residual)<=budget,
        ranks=ranks,rows=count,steps=len(clocks),budget=budget,max_normalized_residual=worst,
        max_reported_residual_inconsistency=row_residual,integrated_channels_J=totals,
        scope='Compression pressure/viscosity and implicit parabolic operator terms; algebraic endpoint reservoir changes explicit; not whole-step coupled energy conservation')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('directory',type=Path)
    p.add_argument('--ranks',type=int,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    try:result=check(a.directory,a.ranks)
    except (ValueError,KeyError,TypeError,OSError) as e:result=dict(passed=False,error=str(e))
    a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
    raise SystemExit(0 if result['passed'] else 1)
