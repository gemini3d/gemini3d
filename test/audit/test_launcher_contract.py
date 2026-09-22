"""Exercise the real launcher without allowing the child to modify simulation output."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

import h5py

from run_native_probes import BASE


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--exe", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    args = parser.parse_args()
    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    sim = work / "simulation's $data with spaces"
    (sim / "inputs").mkdir(parents=True, exist_ok=True)
    (sim / "inputs/config.nml").write_text(BASE)
    with h5py.File(sim / "inputs/simsize.h5", "w") as handle:
        for name, size in (("lx1", 8), ("lx2", 9), ("lx3", 6)):
            handle[name] = size
    for name in ("20130220_18000.000000.h5", "20130220_18060.000000.h5", "restart.h5"):
        (sim / name).write_bytes(("pre-existing output " + name).encode())
    child = work / "child's $executable with spaces"
    record = work / "argv.json"
    child.write_text(f"#!{sys.executable}\nimport json,os,sys\n"
                     "from pathlib import Path\n"
                     "Path(os.environ['ARGV_RECORD']).write_text(json.dumps(sys.argv[1:]))\n"
                     "sys.exit(int(os.environ.get('CHILD_EXIT', '0')))\n")
    child.chmod(0o755)
    mpi = work / "MPI's launcher with spaces"
    mpi.write_text(f"#!{sys.executable}\nimport subprocess,sys\n"
                   "assert sys.argv[1:3] == ['-n', '6'], sys.argv\n"
                   "sys.exit(subprocess.call(sys.argv[3:]))\n")
    mpi.chmod(0o755)

    def hashes():
        return {str(path.relative_to(sim)): hashlib.sha256(path.read_bytes()).hexdigest()
                for path in sim.rglob("*") if path.is_file()}

    before = hashes()
    for cpus in ("1", "6"):
        env = dict(os.environ, GEMINI_CPU=cpus, ARGV_RECORD=str(record))
        for flags in ([], ["-dryrun"], ["-plan"], ["-start_time", "2013", "2", "20", "18060"]):
            for child_exit in ("0", "7"):
                if record.exists():
                    record.unlink()
                result = subprocess.run([str(args.exe), str(sim), "-exe", str(child), "-mpiexec", str(mpi), *flags],
                                        env=dict(env, CHILD_EXIT=child_exit), capture_output=True, text=True, timeout=10)
                assert hashes() == before, "launcher modified pre-existing output"
                assert (result.returncode == 0) == ("-plan" in flags or child_exit == "0"), result.stderr
                assert f'"Ncpu": {cpus}' in result.stdout, result.stdout
                if "-plan" in flags:
                    assert not record.exists(), "planning launched the simulation"
                else:
                    assert json.loads(record.read_text()) == [str(sim), *flags]
    print("16 launcher preservation, quoting, failure, planning and single-CPU contracts passed")


if __name__ == "__main__":
    main()
