"""Build an observational TR-stage trace copy without changing the native algorithm."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shlex
import subprocess

p=argparse.ArgumentParser()
p.add_argument('--source',type=Path,required=True)
p.add_argument('--build',type=Path,required=True)
p.add_argument('--probe',type=Path,required=True)
p.add_argument('--work',type=Path,required=True)
a=p.parse_args();a.build=a.build.absolute();a.work=a.work.absolute();a.work.mkdir(parents=True,exist_ok=False)
original=a.source.read_text()
needle="call check_solver_info(info,'TRBDF21D BDF2 stage')"
assert original.count(needle)==1
modified=original.replace(needle,needle+'\ncall trace_stage(Ts,TR,TRBDF21D)')
trace='''
subroutine trace_stage(old,tr,new)
real(wp),intent(in)::old(:),tr(:),new(:)
integer::u,i,stat
character(4096)::path
call get_environment_variable('GEMINI_PHYSICAL_TR_TRACE',path,status=stat)
if(stat/=0.or.len_trim(path)==0)return
open(newunit=u,file=trim(path),status='unknown',position='append')
do i=1,size(old)
  write(u,'(i0,3(",",es25.17))')i,old(i),tr(i),new(i)
enddo
close(u)
end subroutine
'''
assert modified.count('end module PDEparabolic')==1
modified=modified.replace('end module PDEparabolic',trace+'\nend module PDEparabolic')
source=a.work/'PDEparabolic_trace.f90';source.write_text(modified)
obj=a.work/'PDEparabolic_trace.o';exe=a.work/'physics_trace'
# Use the compiler's actual native flags, including Debug bounds checking.
commands=subprocess.check_output(['ninja','-C',str(a.build),'-t','commands','PDEparabolic'],text=True)
compile_tokens=None
for line in commands.splitlines():
    tok=shlex.split(line)
    if '-c' in tok and '-o' in tok and Path(tok[tok.index('-o')+1]).name=='PDEparabolic.f90.o':
        compile_tokens=tok
assert compile_tokens is not None
cc=compile_tokens
cc[cc.index('-c')+1]=str(source)
cc[cc.index('-o')+1]=str(obj)
intel=Path(cc[0]).name in ('ifx','ifort')
flags=[]
skip=False
for t in cc:
    if skip:
        skip=False
        continue
    if t in ('-J','-module'):
        skip=True
    elif t!='-fpreprocessed' and not t.startswith(('-J','-module=')):
        flags.append(t)
cc=flags
cc.append('-I'+str(a.build/'include'))
if intel:
    cc.extend(['-module',str(a.work)])
else:
    cc.extend(['-ffree-line-length-none','-J'+str(a.work)])
with (a.work/'compile.log').open('w') as f:subprocess.run(cc,cwd=a.build,stdout=f,stderr=subprocess.STDOUT,check=True)
commands=subprocess.check_output(['ninja','-C',str(a.build),'-t','commands','qualification_continuity_floor'],text=True)
tokens=shlex.split(commands.splitlines()[-1]);start=tokens.index('&&')+1;end=tokens.index('&&',start)
cmd=tokens[start:end]
libraries=[t for t in cmd if Path(t).name in ('libgemini3d.a','libgemini3d.debug.a')]
assert len(libraries)==1
library=libraries[0]
cmd=[str(a.probe.absolute()) if t.endswith('continuity_floor.f90.o') else t for t in cmd]
cmd.insert(cmd.index(library),str(obj))
cmd[cmd.index('-o')+1]=str(exe);cmd=[t for t in cmd if not t.startswith('-Wl,--dependency-file=')]
with (a.work/'link.log').open('w') as f:subprocess.run(cmd,cwd=a.build,stdout=f,stderr=subprocess.STDOUT,check=True)
tracepath=a.work/'tr-stages.csv'
r=subprocess.run([str(exe),str(a.work)],capture_output=True,text=True,env=dict(os.environ,GEMINI_PHYSICAL_TR_TRACE=str(tracepath)))
(a.work/'native.log').write_text(r.stdout+r.stderr)
baseline=a.probe.parent
comparisons={}
for q in sorted(baseline.glob('*.csv')):
    copy=a.work/q.name
    comparisons[q.name]={'bitwise_equal':q.read_bytes()==copy.read_bytes(),
                         'sha256':hashlib.sha256(copy.read_bytes()).hexdigest()}
result={'returncode':r.returncode,'observational_copy':True,
        'original_source_sha256':hashlib.sha256(original.encode()).hexdigest(),
        'traced_source_sha256':hashlib.sha256(modified.encode()).hexdigest(),
        'compile_command':cc,'link_command':cmd,'comparisons':comparisons,
        'all_unmodified_outputs_bitwise_equal':all(v['bitwise_equal'] for v in comparisons.values())}
(a.work/'trace-receipt.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k in ('returncode','all_unmodified_outputs_bitwise_equal')}))
if r.returncode or not result['all_unmodified_outputs_bitwise_equal']:raise SystemExit(1)
