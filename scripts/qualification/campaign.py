"""Event lineage, matched rollout evaluation and conservative UQ admission.

These components are software-tested. A synthetic test or a valid manifest is
not evidence of physical fidelity, population coverage or a trained surrogate.
"""
import hashlib
import math
from pathlib import Path
import re
import numpy as np
from evaluate_gates import evaluate

SPLITS={'train','validation','calibration','test'}


def artifact(root,record):
    root=Path(root).resolve();rel=Path(record['path']);p=(root/rel).resolve()
    if rel.is_absolute() or not p.is_relative_to(root) or not p.is_file():raise ValueError('Invalid corpus artifact path')
    with p.open('rb') as f:digest=hashlib.file_digest(f,'sha256').hexdigest()
    if digest!=record['sha256']:raise ValueError('Corpus artifact hash mismatch')
    return p


def validate_manifest(manifest,root):
    if manifest.get('schema')!='gemini.corpus.1':raise ValueError('Unsupported corpus schema')
    if not re.fullmatch('[0-9a-f]{40}',manifest.get('candidate_commit','')):raise ValueError('Exact solver commit required')
    if not manifest.get('target_population'):raise ValueError('Target population required')
    seen=set();groups={};splits=set();counts={s:0 for s in SPLITS}
    for event in manifest['events']:
        eid=event['id'];split=event['split']
        if not eid or eid in seen or split not in SPLITS:raise ValueError('Duplicate event or unknown split')
        seen.add(eid);splits.add(split);counts[split]+=1
        if event.get('solver_commit')!=manifest['candidate_commit']:raise ValueError('Mixed solver commits')
        # Parent event/forcing lineage survives crops, renamed replicas and re-runs.
        for kind in ['physical_event','lineage','forcing_sha256','initial_state_sha256','data_sha256']:
            value=event.get(kind)
            if not isinstance(value,str) or not value:raise ValueError('Missing event identity: '+kind)
            if kind.endswith('sha256') and not re.fullmatch('[0-9a-f]{64}',value):raise ValueError('Invalid SHA256')
            key=(kind,value)
            if key in groups and groups[key]!=split:raise ValueError('Event/forcing/content leakage across splits')
            groups[key]=split
        if event['artifact']['sha256']!=event['data_sha256']:raise ValueError('Inconsistent content identity')
        artifact(root,event['artifact'])
        for kind in ['inputs','grid','forcing','initial_state']:artifact(root,event[kind])
        if event['forcing']['sha256']!=event['forcing_sha256']:raise ValueError('Forcing identity mismatch')
        if event['initial_state']['sha256']!=event['initial_state_sha256']:raise ValueError('Initial-state identity mismatch')
    if splits!=SPLITS:raise ValueError('All four nonempty event splits are required')
    training={e['id'] for e in manifest['events'] if e['split']=='train'}
    if set(manifest.get('normalization_fit_events',[]))!=training:raise ValueError('Normalization must use training events only')
    return dict(schema='gemini.corpus.check.1',passed=True,counts=counts,events=len(seen),
                scope='Content and split integrity; physical validation is a separate prerequisite')


def admit_corpus(manifest,root,register,package_root):
    result=validate_manifest(manifest,root)
    if manifest['candidate_commit']!=register['candidate_commit']:raise ValueError('Corpus candidate differs from qualification')
    evaluation=evaluate(register,Path(package_root))
    blockers=[g for g in evaluation['open_gates'] if g['id'] in ('R02','R04','R07')]
    if blockers:raise ValueError('Scientific corpus admission blocked by '+','.join(g['id'] for g in blockers))
    return dict(**result,prerequisite_records_passed=True,
                limitation='Release owner must authenticate the independent approval records; hashes do not authenticate a reviewer')


def rollout_metrics(truth,prediction,scale):
    y=np.asarray(truth,float);p=np.asarray(prediction,float);s=np.asarray(scale,float)
    if y.shape!=p.shape or y.ndim<2 or not y.size or s.shape!=(y.shape[-1],):raise ValueError('Matched nonempty full rollout/channel shapes required')
    if not all(np.isfinite(x).all() for x in (y,p,s)) or np.any(s<=0):raise ValueError('Finite predictions and positive frozen scales required')
    error=(p-y)/s
    return dict(rms_by_channel=np.sqrt(np.mean(error**2,axis=tuple(range(y.ndim-1)))).tolist(),
                trajectory_score=float(np.max(np.abs(error))),frames=y.shape[0])


def matched_benchmark(truth,initial,candidates,scale,budget):
    """Evaluate complete externally generated rollouts on identical targets.

    Each candidate includes prediction, end-to-end_seconds and model/data hashes;
    deployment latency must be measured on the eventual target, not inferred here.
    """
    y=np.asarray(truth,float);x=np.asarray(initial,float)
    if x.shape!=y.shape[1:]:raise ValueError('Initial state and rollout shape mismatch')
    if not candidates:raise ValueError('At least one actual candidate is required')
    baseline=rollout_metrics(y,np.broadcast_to(x,y.shape),scale)
    results={}
    for name,item in candidates.items():
        for key in ['model_sha256','data_sha256']:
            if not re.fullmatch('[0-9a-f]{64}',item.get(key,'')):raise ValueError('Candidate provenance required')
        seconds=item['end_to_end_seconds']
        if not math.isfinite(seconds) or seconds<=0:raise ValueError('Measured runtime required')
        m=rollout_metrics(y,item['prediction'],scale)
        for key in ('normalized_rms_max','wall_seconds_max'):
            value=budget.get(key)
            if not isinstance(value,(float,int)) or not math.isfinite(value) or value<=0:
                raise ValueError('Positive finite benchmark budgets required')
        results[name]=dict(**m,end_to_end_seconds=seconds,
            passed=max(m['rms_by_channel'])<=budget['normalized_rms_max'] and seconds<=budget['wall_seconds_max'],
            persistence_rms_ratio=max(m['rms_by_channel'])/max(max(baseline['rms_by_channel']),1e-300))
    return dict(persistence=baseline,candidates=results,passed=all(r['passed'] for r in results.values()))


def conformal(scores,event_ids,alpha):
    scores=np.asarray(scores,float)
    if scores.ndim!=1 or len(scores)==0 or len(event_ids)!=len(scores) or len(set(event_ids))!=len(scores):
        raise ValueError('One independent trajectory score per unique event is required')
    if not np.isfinite(scores).all() or np.any(scores<0) or not 0<alpha<1:raise ValueError('Invalid calibration scores or alpha')
    k=math.ceil((len(scores)+1)*(1-alpha))
    radius=float(np.sort(scores)[k-1]) if k<=len(scores) else math.inf
    return dict(radius=radius,events=len(scores),order_statistic=k,alpha=alpha,
                assumption='Exchangeable independent events in the declared population; not adjacent time windows')


def uq_decision(calibration,features,training_min,training_max,width_cap,population,calibrated_population):
    x,lo,hi=[np.asarray(a,float) for a in (features,training_min,training_max)]
    if x.shape!=lo.shape or x.shape!=hi.shape or not all(np.isfinite(a).all() for a in (x,lo,hi)) or np.any(lo>hi):
        raise ValueError('Finite matched shift-detection features required')
    if not math.isfinite(width_cap) or width_cap<=0:raise ValueError('Positive predeclared width cap required')
    reasons=[];radius=calibration['radius']
    if not math.isfinite(radius):reasons.append('insufficient_calibration')
    elif radius<0:raise ValueError('Negative interval radius')
    elif 2*radius>width_cap:reasons.append('interval_too_wide')
    if np.any((x<lo)|(x>hi)):reasons.append('outside_training_envelope')
    if not population or population!=calibrated_population:reasons.append('population_mismatch')
    return dict(status='abstain' if reasons else 'eligible_for_research',reasons=reasons,control_enabled=False,
                limitation='Envelope test detects only declared feature excursions; it is not a proof of no distribution shift')


def coverage(scores,event_ids,calibration,calibration_event_ids):
    if set(event_ids)&set(calibration_event_ids):raise ValueError('Calibration/test event leakage')
    x=np.asarray(scores,float)
    if x.ndim!=1 or not len(x) or len(event_ids)!=len(x) or len(set(event_ids))!=len(x) or not np.isfinite(x).all() or np.any(x<0):
        raise ValueError('Independent finite nonnegative test-event scores required')
    if not math.isfinite(calibration['radius']) or calibration['radius']<0:
        raise ValueError('Finite nonnegative interval radius required')
    if len(set(calibration_event_ids))!=len(calibration_event_ids) or len(calibration_event_ids)!=calibration['events']:
        raise ValueError('Calibration event inventory mismatch')
    return dict(events=len(x),covered=int(np.sum(x<=calibration['radius'])),
                empirical_coverage=float(np.mean(x<=calibration['radius'])),width=2*calibration['radius'],
                scope='Empirical held-out event coverage; no population confidence claim')


def qualify_coverage(scores,event_ids,calibration,calibration_event_ids,*,coverage_target,width_cap,confidence):
    """One-sided exact binomial lower bound for independent held-out events.

    Targets and confidence must be fixed before evaluating test events. This
    does not establish event independence or exchangeability from arrays.
    """
    from scipy.stats import beta
    if not (0<coverage_target<1 and 0<confidence<1 and math.isfinite(width_cap) and width_cap>0):
        raise ValueError('Predeclared coverage, width and confidence criteria required')
    result=coverage(scores,event_ids,calibration,calibration_event_ids)
    k,n=result['covered'],result['events']
    lower=0. if k==0 else float(beta.ppf(1-confidence,k,n-k+1))
    return dict(**result,confidence=confidence,coverage_target=coverage_target,width_cap=width_cap,
                population_coverage_lower_bound=lower,method='One-sided exact Clopper-Pearson',
                passed=bool(lower>=coverage_target and result['width']<=width_cap),
                assumptions='Independent held-out events from the fixed target population; no test-dependent tuning')
