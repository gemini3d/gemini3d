"""Native synthetic checkpoint identity/finite-state tests, not a science run."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import h5py
import numpy as np

from research_matrix import checkpoint_shapes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exe', type=Path, required=True)
    parser.add_argument('--mpiexec', required=True)
    parser.add_argument('--numproc-flag', default='-n')
    parser.add_argument('--mpi-pref', action='append', default=[])
    parser.add_argument('--mpi-post', action='append', default=[])
    parser.add_argument('--layout', choices=('1x2', '2x2'), required=True)
    parser.add_argument('--work', type=Path, required=True)
    args = parser.parse_args()
    layout = tuple(map(int, args.layout.split('x')))
    ranks = int(np.prod(layout))
    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    rejected_work = work / 'invalid-matrix-options'
    invalid = subprocess.run(
        [sys.executable, str(Path(__file__).with_name('research_matrix.py')),
         '--build', str(args.exe.parent), '--inputs-root', str(work / 'absent-inputs'),
         '--work', str(rejected_work), '--quick', '--extended-layouts'],
        text=True, capture_output=True, timeout=30)
    assert invalid.returncode == 2 and '--extended-layouts requires the full matrix' in invalid.stderr
    assert not rejected_work.exists()
    env = dict(os.environ, GEMINI_EXACT_RESTART='1', GEMINI_INPUT_SHA256='1' * 64,
               GEMINI_EXECUTABLE_SHA256=hashlib.sha256(args.exe.read_bytes()).hexdigest())
    rows = []

    def run(mode, case, expected=None):
        proc = subprocess.run(
            [args.mpiexec, args.numproc_flag, str(ranks), *args.mpi_pref,
             str(args.exe.resolve()), *args.mpi_post, mode, str(case),
             *map(str, layout)], env=env, text=True, capture_output=True, timeout=30)
        output = proc.stdout + proc.stderr
        (work / f'{case.name}-{mode}.log').write_text(output)
        passed = proc.returncode == 0 if expected is None else proc.returncode != 0 and expected in output
        rows.append(dict(case=case.name, mode=mode, passed=passed,
                         expected_rejection=expected, returncode=proc.returncode))
        if not passed:
            raise AssertionError(f'{case.name} {mode}: expected {expected!r}\n{output}')

    def fresh(name, source=None):
        path = work / name
        if path.exists():
            shutil.rmtree(path)
        if source is None:
            path.mkdir()
        else:
            shutil.copytree(source, path)
        return path

    def frame(case):
        return next(case.glob('????????_?????.??????.h5'))

    def sidecar(case, rank=0):
        with h5py.File(frame(case)) as f:
            prefix = f['restart_runtime/prefix'][()].decode().strip()
        return case / f'{prefix}.r{rank:08d}.h5'

    seed = fresh('seed')
    run('write', seed)
    run('read', seed)
    assert checkpoint_shapes(seed, ranks, layout)
    copied = fresh('copied', seed)
    run('read', copied)
    with h5py.File(frame(seed)) as root:
        prefix = root['restart_runtime/prefix'][()].decode().strip()
        assert root['restart_runtime/schema'][()] == 2
        assert root['restart_runtime/checkpoint'][()].decode().strip() == frame(seed).name
        for rank in range(ranks):
            with h5py.File(sidecar(seed, rank)) as state:
                assert state['schema'][()] == 2
                assert state['rank'][()] == rank
                assert tuple(state['layout'][...]) == layout
                assert state['generation'][()].decode().strip() == prefix
                assert state['checkpoint'][()] == root['restart_runtime/checkpoint'][()]
                assert state['input_sha256'][()] == root['restart_runtime/input_sha256'][()]
                assert state['executable_sha256'][()] == root['restart_runtime/executable_sha256'][()]

    bad = fresh('swap_rank', seed)
    first, second = sidecar(bad), sidecar(bad, ranks - 1)
    first_bytes, second_bytes = first.read_bytes(), second.read_bytes()
    first.write_bytes(second_bytes)
    second.write_bytes(first_bytes)
    run('read', bad, 'Runtime state rank mismatch')

    # Same root timestamp, same rank, same inputs/executable, different publication.
    other = fresh('other_generation')
    run('write', other)
    bad = fresh('mismatched_generation', seed)
    shutil.copyfile(sidecar(other), sidecar(bad))
    run('read', bad, 'Runtime state generation mismatch')

    mutations = [
        ('root_schema', 'restart_runtime/schema', 1, 'Unsupported runtime checkpoint schema', True),
        ('root_identity', 'restart_runtime/checkpoint', b'other.h5', 'root identity mismatch', True),
        ('root_layout', 'restart_runtime/layout', [2, 1], 'Incompatible runtime checkpoint layout', True),
        ('root_executable', 'restart_runtime/executable_sha256', b'0' * 64, 'Changed restart executable', True),
        ('state_schema', 'schema', 1, 'Unsupported runtime state schema', False),
        ('state_root', 'checkpoint', b'other.h5', 'Runtime state root identity mismatch', False),
        ('state_layout', 'layout', [2, 1], 'Runtime state layout mismatch', False),
        ('state_executable', 'executable_sha256', b'0' * 64, 'Changed restart executable', False),
        ('state_input', 'input_sha256', b'0' * 64, 'Changed restart inputs', False),
        ('incomplete', 'complete', 0, 'Incomplete runtime checkpoint', False),
        ('clock', 'dt', np.nan, 'Invalid runtime checkpoint clock', False),
    ]
    for name, key, value, expected, root in mutations:
        bad = fresh(name, seed)
        with h5py.File(frame(bad) if root else sidecar(bad, ranks - 1), 'r+') as f:
            f[key][...] = value
        run('read', bad, expected)

    bad = fresh('missing_state', seed)
    sidecar(bad, ranks - 1).unlink()
    run('read', bad, 'Missing runtime state file')
    bad = fresh('bad_dtype', seed)
    with h5py.File(sidecar(bad, ranks - 1), 'r+') as f:
        values = f['fluid'][...].astype('float32')
        del f['fluid']
        f['fluid'] = values
    run('read', bad, 'requires float64 rank-four state')
    for field in ('fluid', 'electro', 'vi1', 'vi2', 'vi3'):
        for label, value in (('nan', np.nan), ('inf', np.inf), ('neginf', -np.inf)):
            bad = fresh(f'{field}_{label}', seed)
            with h5py.File(sidecar(bad, ranks - 1), 'r+') as f:
                f[field][tuple(n - 1 for n in f[field].shape)] = value
            assert not checkpoint_shapes(bad, ranks, layout)
            run('read', bad, f'Nonfinite runtime checkpoint payload: {field}')
    bad = fresh('save_nonfinite', seed)
    original = frame(bad).read_bytes()
    run('write_nan', bad, 'Nonfinite runtime checkpoint payload: vi3')
    assert frame(bad).read_bytes() == original
    run('read', bad)
    (work / 'results.json').write_text(json.dumps(dict(layout=layout, ranks=ranks, runs=rows), indent=2) + '\n')
    print(f'{len(rows)} native checkpoint checks passed on {args.layout}; synthetic state only.')


if __name__ == '__main__':
    main()
