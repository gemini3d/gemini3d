import typing as T
import os
import subprocess
import shutil
from pathlib import Path
from .utils import get_mpi_count

Pathlike = T.Union[str, Path]
cwd = os.getcwd()


def runner(mpiexec: Pathlike, gemexe: Pathlike, config_file: Pathlike, out_dir: Path):

    if not mpiexec:
        mpiexec = "mpiexec"
    mpi_root = os.getenv("MPI_ROOT")
    if mpi_root:
        mpi_root += "/bin"
    mpiexec = shutil.which(mpiexec, path=mpi_root)
    if not mpiexec:
        raise FileNotFoundError("Need mpiexec to run simulations")
    print("mpiexec:", mpiexec)

    if gemexe:
        gemexe = Path(gemexe).expanduser()
        if not gemexe.is_file():
            raise FileNotFoundError(gemexe)
    else:
        build_dir = Path(__file__).resolve().parents[1] / 'build'
        if not build_dir.is_dir():
            raise NotADirectoryError(build_dir)
        gemexe = shutil.which("gemini.bin", path=str(build_dir))
        if not gemexe:
            raise FileNotFoundError("Cannot find gemini.bin")
    print("gemini executable:", gemexe)

    out_dir = Path(out_dir).expanduser()
    if out_dir.is_file():
        raise NotADirectoryError(out_dir)
    if not out_dir.is_dir():
        out_dir.mkdir(parents=True, exist_ok=True)

    Nmpi = get_mpi_count(config_file, cwd=cwd)

    cmd = [str(mpiexec), "-n", str(Nmpi), str(gemexe), str(config_file), str(out_dir)]
    print(" ".join(cmd))
    ret = subprocess.run(cmd)

    raise SystemExit(ret.returncode)
