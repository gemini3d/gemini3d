"""Check authenticated GitHub Actions API records against the exact candidate.

Read-only collector. Does not push a branch, start jobs or post a review. The
raw API responses are retained. This establishes hosted job status, not physics
approval or leak freedom; those require their own measured artifacts.
"""
import argparse,json,re,subprocess
from pathlib import Path

REQUIRED={'native (Debug)','native (Release)','current-hdf5','sanitizers'}
WORKFLOW='.github/workflows/preintegration.yml'
REQUIRED_STEPS={
    'native (Debug)':{'Run native suite and reference comparisons','Verify restart record precision and rejection behavior'},
    'native (Release)':{'Run native suite and reference comparisons','Verify restart record precision and rejection behavior',
                        'Qualify exact restart and native numerical budgets'},
    'current-hdf5':{'Build and test current-library profile','Qualify current-library restart and numerical budgets',
                    'Measure complete workload with kernel memory accounting'},
    'sanitizers':{'Check audit probes and native C++ consumer','Verify LeakSanitizer detects the intentional leak control',
                   'Check leaks in standalone ABI and interpolation probes',
                   'Check instrumented full-step Fortran and C++ checkpoints','Check leaks in the actual native application'}}


def verify(run,jobs,repository,commit,queried_attempt):
    if not re.fullmatch(r'[0-9a-f]{40}',commit):raise ValueError('Exact candidate commit required')
    if not re.fullmatch(r'[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+',repository):raise ValueError('Explicit owner/repository required')
    issues=[]
    if run.get('repository',{}).get('full_name','').lower()!=repository.lower():issues.append('repository_mismatch')
    if run.get('head_sha')!=commit:issues.append('commit_mismatch')
    if run.get('path','').split('@')[0]!=WORKFLOW:issues.append('workflow_mismatch')
    if run.get('status')!='completed' or run.get('conclusion')!='success':issues.append('workflow_not_successful')
    if not isinstance(run.get('run_attempt'),int) or run['run_attempt']<1:issues.append('invalid_run_attempt')
    if queried_attempt!=run.get('run_attempt'):issues.append('query_attempt_mismatch')
    seen={}
    for job in jobs:
        name=job.get('name');seen.setdefault(name,[]).append(job)
        if job.get('run_id')!=run.get('id') or job.get('head_sha')!=commit:issues.append('job_identity_mismatch')
        # The API's attempt-specific endpoint supplies the binding even when
        # older job response schemas omit an individual run_attempt field.
        if 'run_attempt' in job and job['run_attempt']!=run.get('run_attempt'):issues.append('job_attempt_mismatch')
        if name in REQUIRED:
            if job.get('status')!='completed' or job.get('conclusion')!='success':issues.append('job_not_successful:'+name)
            steps=job.get('steps',[])
            if not steps or any(s.get('conclusion') in ('failure','cancelled','timed_out') for s in steps):
                issues.append('failed_or_missing_steps:'+name)
            for required in REQUIRED_STEPS[name]:
                matches=[s for s in steps if s.get('name')==required]
                if len(matches)!=1 or matches[0].get('conclusion')!='success':
                    issues.append('required_step_not_successful:'+name+':'+required)
    for name in sorted(REQUIRED):
        if len(seen.get(name,[]))!=1:issues.append('missing_or_duplicate_job:'+name)
    return dict(schema='gemini.hosted.evidence.1',passed=not issues,issues=issues,repository=repository,
        candidate_commit=commit,run_url=run.get('html_url'),run_id=run.get('id'),run_attempt=run.get('run_attempt'),
        required_jobs=sorted(REQUIRED),limitation='Hosted execution record only; result artifacts and independent science gates remain required')


def api(endpoint):
    p=subprocess.run(['gh','api',endpoint],capture_output=True,text=True,check=True,timeout=60)
    return json.loads(p.stdout)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--repository',required=True)
    p.add_argument('--commit',required=True);p.add_argument('--run-id',type=int,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if not re.fullmatch(r'[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+',a.repository) or a.run_id<=0:
        p.error('Explicit owner/repository and positive run ID required')
    a.output.mkdir(parents=True,exist_ok=False)
    run=api(f'repos/{a.repository}/actions/runs/{a.run_id}')
    (a.output/'run.json').write_text(json.dumps(run,indent=2)+'\n')
    jobs=[];page=1
    while True:
        payload=api(f'repos/{a.repository}/actions/runs/{a.run_id}/attempts/{run["run_attempt"]}/jobs?per_page=100&page={page}')
        (a.output/f'jobs-{page}.json').write_text(json.dumps(payload,indent=2)+'\n')
        jobs.extend(payload['jobs'])
        if len(jobs)>=payload['total_count']:break
        if not payload['jobs']:raise ValueError('Incomplete job pagination')
        page+=1
    result=verify(run,jobs,a.repository,a.commit,run['run_attempt'])
    (a.output/'result.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
    return 0 if result['passed'] else 2
if __name__=='__main__':raise SystemExit(main())
