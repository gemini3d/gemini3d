"""Measure a complete command tree in a fresh delegated Linux cgroup v2.

No sampling fallback and no mount/privilege changes. cgroup memory.peak includes
charged file cache and kernel memory; it is not the sum of process RSS samples.
An unavailable delegation is a blocked measurement (exit 2), never a pass.
"""
import argparse
import datetime
import json
import math
import os
from pathlib import Path
import subprocess
import time
import uuid


def measure(command,parent,output,wall_budget,memory_budget_mib,timeout):
    if not command or not all(math.isfinite(x) and x>0 for x in (wall_budget,memory_budget_mib,timeout)):
        raise ValueError('Positive finite budgets and command required')
    group=parent/('gemini-qualification-'+uuid.uuid4().hex)
    result=dict(schema='gemini.kernel.memory.1',passed=False,status='blocked',command=command,
                host_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                memory_budget_mib=memory_budget_mib,wall_budget_seconds=wall_budget)
    proc=None
    try:
        group.mkdir()
        required=['cgroup.procs','cgroup.events','memory.peak','memory.events','memory.max']
        if any(not (group/k).is_file() for k in required):raise OSError('Delegated cgroup lacks memory controller')
        if (group/'memory.peak').read_text().strip()!='0':raise OSError('Fresh cgroup has preexisting memory usage')
        (group/'memory.max').write_text(str(int(memory_budget_mib*2**20)))
        if (group/'memory.swap.max').exists():(group/'memory.swap.max').write_text('0')
        # Move the new child itself before exec, before any command work starts.
        def enter():
            with (group/'cgroup.procs').open('w') as stream:stream.write(str(os.getpid()))
        start=time.monotonic()
        with output.with_suffix('.log').open('w') as log:
            proc=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT,preexec_fn=enter,start_new_session=True)
            try:code=proc.wait(timeout=timeout)
            except subprocess.TimeoutExpired:
                if (group/'cgroup.kill').exists():(group/'cgroup.kill').write_text('1')
                else:os.killpg(proc.pid,9)
                proc.wait();raise RuntimeError('Workload exceeded measurement timeout')
        elapsed=time.monotonic()-start
        peak=int((group/'memory.peak').read_text())
        events=dict(line.split() for line in (group/'memory.events').read_text().splitlines())
        populated=dict(line.split() for line in (group/'cgroup.events').read_text().splitlines()).get('populated')
        okay=code==0 and 0<peak<=memory_budget_mib*2**20 and elapsed<=wall_budget and populated=='0' and int(events['oom_kill'])==0
        result.update(status='passed' if okay else 'failed',passed=okay,returncode=code,wall_seconds=elapsed,
                      peak_charged_memory_bytes=peak,events=events,populated_after_exit=populated)
    except (OSError,RuntimeError,subprocess.SubprocessError) as e:
        result.update(error=str(e),status='blocked' if proc is None else 'failed')
    finally:
        if group.exists():
            if (group/'cgroup.events').exists() and 'populated 1' in (group/'cgroup.events').read_text():
                if (group/'cgroup.kill').exists():(group/'cgroup.kill').write_text('1')
            try:group.rmdir()
            except OSError as e:result.update(passed=False,status='failed',cleanup_error=str(e))
    output.write_text(json.dumps(result,indent=2)+'\n')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--parent',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--wall-budget',type=float,default=600)
    p.add_argument('--memory-budget-mib',type=float,default=2048);p.add_argument('--timeout',type=float,default=660)
    p.add_argument('command',nargs=argparse.REMAINDER);a=p.parse_args()
    command=a.command[1:] if a.command[:1]==['--'] else a.command
    result=measure(command,a.parent,a.output,a.wall_budget,a.memory_budget_mib,a.timeout)
    print(json.dumps(result,indent=2));raise SystemExit(0 if result['passed'] else 2)
