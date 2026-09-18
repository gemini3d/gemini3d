"""Probe required host capabilities without changing privileges or controls.

Leak qualification requires a positive control: detecting an intentional leak
and accepting a leak-free program. A successful compiler exit is insufficient.
"""
import argparse,hashlib,json,os,platform,subprocess,tempfile
from pathlib import Path
from kernel_memory import measure


def leak_probe(cc,work):
    work=Path(work);work.mkdir(parents=True,exist_ok=False)
    code='''#include <stdlib.h>
__attribute__((noinline)) static void allocate(int leak) {
    volatile char *p = malloc(1237);
    if (!p) abort();
    p[0] = 1;
    if (!leak) free((void *)p);
}
int main(int argc, char **argv) { (void)argv; allocate(argc > 1); return 0; }
'''
    source=work/'leak_control.c';source.write_text(code);exe=work/'leak_control'
    command=[cc,'-O0','-g','-fsanitize=address','-fno-omit-frame-pointer',str(source),'-o',str(exe)]
    built=subprocess.run(command,capture_output=True,text=True,timeout=60)
    (work/'compile.log').write_text(built.stdout+built.stderr)
    if built.returncode:return dict(passed=False,status='blocked',reason='Probe compilation failed',command=command)
    env=dict(os.environ,ASAN_OPTIONS='detect_leaks=1:halt_on_error=1',LSAN_OPTIONS='exitcode=23')
    records=[]
    for label,args in [('clean',[]),('intentional_leak',['leak'])]:
        try:
            run=subprocess.run([str(exe),*args],env=env,capture_output=True,text=True,timeout=30)
            output=run.stdout+run.stderr;(work/(label+'.log')).write_text(output)
            okay=(run.returncode==0 and 'LeakSanitizer' not in output) if not args else (
                run.returncode==23 and 'detected memory leaks' in output and '1237 byte(s)' in output)
            records.append(dict(case=label,returncode=run.returncode,passed=okay))
        except (OSError,subprocess.SubprocessError) as e:records.append(dict(case=label,passed=False,error=str(e)))
    return dict(schema='gemini.leak.capability.1',passed=all(r['passed'] for r in records),
                status='capable' if all(r['passed'] for r in records) else 'blocked',runs=records,
                executable_sha256=hashlib.sha256(exe.read_bytes()).hexdigest(),
                limitation='Runtime controls only; actual instrumented GEMINI application runs are separately required')


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--cc',default='gcc')
    p.add_argument('--work',type=Path,required=True);p.add_argument('--only',choices=['all','leak','memory'],default='all');p.add_argument('--cgroup-parent',type=Path,default=Path('/sys/fs/cgroup'))
    a=p.parse_args();a.work.mkdir(parents=True,exist_ok=False)
    leaks=leak_probe(a.cc,a.work/'leak') if a.only in ('all','leak') else None
    memory=measure(['/bin/true'],a.cgroup_parent,a.work/'memory.json',600,2048,30) if a.only in ('all','memory') else None
    result=dict(schema='gemini.host.capabilities.1',platform=platform.platform(),leak=leaks,memory=memory,
                requested=a.only,passed=all(r['passed'] for r in (leaks,memory) if r is not None))
    (a.work/'result.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
    return 0 if result['passed'] else 2
if __name__=='__main__':raise SystemExit(main())
