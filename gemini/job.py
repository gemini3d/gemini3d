import typing as T
import os
import logging
import subprocess
import shutil
from pathlib import Path

from .utils import get_mpi_count
from .config import read_config
from .hpc import hpc_job, hpc_batch_create

Pathlike = T.Union[str, Path]
cwd = os.getcwd()


def runner(mpiexec: Pathlike, gemexe: Pathlike, config_file: Pathlike, out_dir: Path) -> int:

    config_file = Path(config_file).resolve()
    # load configuration to know what directories to check
    p = read_config(config_file)
    for k in ("indat_size", "indat_grid", "indat_file"):
        if not p[k].is_file():
            ok = initialize_simulation(config_file, p)
            if not ok or not p[k].is_file():
                raise RuntimeError("could not initialize simulation. Try doing this manually.")
            break

    if p.get("flagE0file") == 1:
        E0dir = p["E0dir"].resolve()
        if not E0dir.is_dir():
            raise FileNotFoundError(E0dir)
    if p.get("flagprecfile") == 1:
        precdir = p["precdir"].resolve()
        if not precdir.is_dir():
            raise FileNotFoundError(precdir)

    # check MPIexec
    mpiexec = check_mpiexec(mpiexec)
    logging.info(f"Using mpiexec: {mpiexec}")

    gemexe = check_gemini_exe(gemexe)
    logging.info(f"using gemini executable: {gemexe}")

    out_dir = check_outdir(out_dir)

    Nmpi = get_mpi_count(config_file, cwd=cwd)

    cmd = [str(mpiexec), "-n", str(Nmpi), str(gemexe), str(config_file), str(out_dir)]
    print(" ".join(cmd), "\n")

    batcher = hpc_job()
    if batcher:
        job_file = hpc_batch_create(batcher, out_dir, cmd)  # noqa: F841
        # hpc_submit_job(job_file)
        print("Please examine batch file", job_file, "and when ready submit the job as usual.")
        ret = 0
    else:
        ret = subprocess.run(cmd).returncode

    return ret


def initialize_simulation(config_file: Path, p: T.Dict[str, T.Any], matlab: Pathlike = None) -> bool:
    """
    TODO: these functions will be in Python

    GNU Octave doesn't have the HDF5 API needed for Gemini
    """

    env = os.environ
    if not os.environ.get("GEMINI_ROOT"):
        env["GEMINI_ROOT"] = Path(__file__).parents[1].as_posix()

    matlab = shutil.which(matlab) if matlab else shutil.which("matlab")
    if not matlab:
        return False

    cmd = [matlab, "-batch", "config()"]
    print("Initializing simulation: ", cmd)

    ret = subprocess.run(cmd, cwd=config_file.parent, env=env)

    return ret.returncode == 0


def check_mpiexec(mpiexec: Pathlike) -> str:
    """ check that specified mpiexec exists on this system """

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
        raise EnvironmentError(msg)

    return mpiexec


def check_gemini_exe(gemexe: Pathlike) -> str:
    """
    check that Gemini exectuable can run on this system

    If not given a specific full path to gemini.bin, looks for gemini.bin under:

        build
        build / Release
        build / Debug
    """

    if gemexe:
        gemexe = Path(gemexe).expanduser()
        if not gemexe.is_file():
            raise EnvironmentError(f"Cannot find gemini.bin in {gemexe}")
    else:
        build_dir = Path(__file__).resolve().parents[1] / "build"
        if not build_dir.is_dir():
            raise EnvironmentError(f"Build directory missing: {build_dir}")

        for d in (build_dir, build_dir / "Release", build_dir / "Debug"):
            gemexe = shutil.which("gemini.bin", path=str(d))
            if gemexe:
                break
        if not gemexe:
            raise EnvironmentError(f"Cannot find gemini.bin under {build_dir}")

    gemexe = str(Path(gemexe).resolve())

    ret = subprocess.run(gemexe, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, universal_newlines=True)
    if ret.returncode != 77:
        raise RuntimeError(
            f"{gemexe} was not runnable on your platform. Try recompiling on this computer type."
            "E.g. different HPC nodes may not have the CPU feature sets."
            f"{ret.stderr}"
        )

    return gemexe


def check_outdir(out_dir: Pathlike) -> Path:

    out_dir = Path(out_dir).expanduser().resolve()
    if out_dir.is_file():
        raise NotADirectoryError(f"please specify output DIRECTORY, you specified {out_dir}")
    if not out_dir.is_dir():
        out_dir.mkdir(parents=True, exist_ok=True)

    return out_dir
