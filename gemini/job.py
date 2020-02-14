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
        msg = "Need mpiexec to run simulations"
        if os.name == "nt":
            msg += "\n\nTypically Windows users will use any one of:"
            msg += "\na) Windows Subsystem for Linux (WSL) <-- recommended \nb) Cygwin"
            msg += "\nc) Intel Parallel Studio for Windows (or WSL)"
        raise FileNotFoundError(msg)
    print("Using mpiexec:", mpiexec)

    if gemexe:
        gemexe = Path(gemexe).expanduser()
        if not gemexe.is_file():
            raise FileNotFoundError(gemexe)
    else:
        build_dir = Path(__file__).resolve().parents[1] / "build"
        if not build_dir.is_dir():
            raise NotADirectoryError(build_dir)
        gemexe = shutil.which("gemini.bin", path=str(build_dir))
        if not gemexe:
            raise FileNotFoundError("Cannot find gemini.bin")
    gemexe = str(gemexe)

    ret = subprocess.run(gemexe, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, universal_newlines=True)
    if ret.returncode != 77:
        raise RuntimeError(
            f"{gemexe} was not runnable on your platform. Try recompiling on this computer type."
            "E.g. different HPC nodes may not have the CPU feature sets."
            f"{ret.stderr}"
        )

    print("using gemini executable", gemexe)

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
