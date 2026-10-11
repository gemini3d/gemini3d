"""Compile probes against unmodified GEMINI objects, preserving execution receipts."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shlex
import subprocess

p = argparse.ArgumentParser()
p.add_argument('--build', type=Path, required=True)
p.add_argument('--source', type=Path, required=True)
p.add_argument('--work', type=Path, required=True)
a = p.parse_args()
a.build = a.build.absolute()
a.source = a.source.absolute()
a.work = a.work.absolute()
a.work.mkdir(parents=True, exist_ok=False)
obj = a.work / 'physics_native.o'
exe = a.work / 'physics_native'
compile_command = [os.environ['FC'], '-g', '-fcheck=all', '-ffree-line-length-none',
                   '-I'+str(a.build/'include'), '-c', str(a.source), '-o', str(obj)]
with (a.work/'compile.log').open('w') as f:
    r = subprocess.run(compile_command, stdout=f, stderr=subprocess.STDOUT)
if r.returncode:
    raise SystemExit(f'Compile failed; see {a.work}/compile.log')
commands = subprocess.check_output(['ninja', '-C', str(a.build), '-t', 'commands',
                                    'qualification_continuity_floor'], text=True)
tokens = shlex.split(commands.splitlines()[-1])
start = tokens.index('&&')+1
end = tokens.index('&&', start)
command = tokens[start:end]
assert sum(t.endswith('continuity_floor.f90.o') for t in command) == 1
command = [str(obj) if t.endswith('continuity_floor.f90.o') else t for t in command]
command[command.index('-o')+1] = str(exe)
command = [t for t in command if not t.startswith('-Wl,--dependency-file=')]
with (a.work/'link.log').open('w') as f:
    r = subprocess.run(command, cwd=a.build, stdout=f, stderr=subprocess.STDOUT)
if r.returncode:
    raise SystemExit(f'Link failed; see {a.work}/link.log')
r = subprocess.run([str(exe), str(a.work)], env=dict(os.environ, GEMINI_NUMERICAL_AUDIT='0'),
                   capture_output=True, text=True)
(a.work/'native.log').write_text(r.stdout+r.stderr)
receipt = {
    'schema': 'gemini.physics.native.receipt.1',
    'native_returncode': r.returncode,
    'probe_source_sha256': hashlib.sha256(a.source.read_bytes()).hexdigest(),
    'probe_executable_sha256': hashlib.sha256(exe.read_bytes()).hexdigest(),
    'compile_command': compile_command,
    'link_command': command,
    'objects': [],
    'files': []
}
for t in command:
    q = a.build/t
    if t.endswith('.o') or t.endswith('.a') or '.so' in t:
        if q.is_file():
            receipt['objects'].append({'file': str(q), 'sha256': hashlib.sha256(q.read_bytes()).hexdigest()})
for q in sorted(a.work.glob('*.csv')):
    receipt['files'].append({'file':q.name, 'bytes':q.stat().st_size,
                            'sha256':hashlib.sha256(q.read_bytes()).hexdigest()})
(a.work/'execution.json').write_text(json.dumps(receipt, indent=2)+'\n')
print(json.dumps({'native_returncode':r.returncode, 'output_files':len(receipt['files']),
                  'work': str(a.work)}))
raise SystemExit(r.returncode)
