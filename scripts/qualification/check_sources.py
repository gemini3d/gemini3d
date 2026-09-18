"""Validate native aggregate ETD production/loss and temperature-floor ledgers.

No inferred residual is used as a source. This does not certify the underlying
reaction/closure physics or the still-incomplete whole-step energy balance.
"""
import argparse
import csv
import json
import math
from pathlib import Path

EXPECTED={(q,s) for q in (1,2,3) for s in range(1,8 if q==3 else 7)}


def rows(path):
    with path.open() as stream:
        for raw in csv.DictReader(stream):
            row={k:float(v) for k,v in raw.items()}
            if not all(math.isfinite(v) for v in row.values()):raise ValueError('Nonfinite ledger')
            for key in ('quantity','species','stage'):
                if key in row and row[key]!=int(row[key]):raise ValueError('Noninteger ledger identifier')
            if row['dt']<=0:raise ValueError('Nonpositive time step')
            yield row


def check(directory,ranks,budget=1e-11):
    if ranks not in (1,2):raise ValueError('Unsupported research layout')
    expected_files={f'sources-r{r:08d}.csv' for r in range(ranks)}
    if {p.name for p in directory.glob('sources-r*.csv')}!=expected_files:
        raise ValueError('Incomplete or unexpected source rank inventory')
    count=0;worst=0.;worst_state=0.;common=None;floors={};floor_rows=0
    for rank in range(ranks):
        stages={}
        for row in rows(directory/f'sources-r{rank:08d}.csv'):
            clock=(row['t'],row['dt']);key=(int(row['quantity']),int(row['species']))
            if key in stages.setdefault(clock,set()):raise ValueError('Duplicate source stage')
            stages[clock].add(key)
            # The native scale integrates absolute state. It must at least
            # dominate the absolute signed integrals it claims to normalize.
            scale=row['scale']
            lower=max(abs(row[k]) for k in ('before','after','production','integrated_loss'))
            if scale<=0 or scale<(1-1e-12)*lower:raise ValueError('Invalid source scale')
            worst=max(worst,abs(row['delta']-row['production']+row['integrated_loss'])/scale)
            worst_state=max(worst_state,abs(row['after']-row['before']-row['delta'])/scale)
            count+=1
        if not stages or any(s!=EXPECTED for s in stages.values()):
            raise ValueError('Missing quantity/species at one or more source steps')
        clocks=set(stages)
        if common is not None and clocks!=common:raise ValueError('MPI source clocks do not match')
        common=clocks
        continuity={(r['t'],r['dt']) for r in rows(directory/f'continuity-r{rank:08d}.csv')}
        if clocks!=continuity:raise ValueError('Source/continuity steps do not match')
        seen=set()
        for row in rows(directory/f'temperature-floor-r{rank:08d}.csv'):
            clock=(row['t'],row['dt']);key=(*clock,int(row['stage']),int(row['species']))
            if key in seen or clock not in clocks:raise ValueError('Duplicate or unknown temperature floor step')
            if row['energy_added_J']<0:raise ValueError('Negative energy attributed to a lower floor')
            seen.add(key);label=f"{int(row['stage'])}:{int(row['species'])}"
            floors[label]=floors.get(label,0.)+row['energy_added_J'];floor_rows+=1
        expected={(*clock,stage,s) for clock in clocks for stage in (1,2) for s in range(1,8)}
        if seen!=expected:raise ValueError('Missing standard-energy temperature-floor stages')
    return dict(schema='gemini.native.sources.1',passed=max(worst,worst_state)<=budget,
                ranks=ranks,rows=count,steps=len(common),budget=budget,max_normalized_residual=worst,
                max_state_delta_inconsistency=worst_state,temperature_floor_rows=floor_rows,
                temperature_floor_energy_J=floors,
                scope='Frozen-coefficient aggregate production/loss ODEs and explicit standard-energy floors; not full coupled energy conservation or reaction-rate validation')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('directory',type=Path)
    p.add_argument('--ranks',type=int,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    try:result=check(a.directory,a.ranks)
    except (ValueError,KeyError,OSError,TypeError) as e:result=dict(passed=False,error=str(e))
    a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
    raise SystemExit(0 if result['passed'] else 1)
