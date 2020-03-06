import typing as T
import os
import logging
import subprocess
import shutil
from pathlib import Path
import functools

from .utils import get_mpi_count
from .config import read_config, get_config_filename
from .hpc import hpc_batch_detect, hpc_batch_create
from .model_setup import model_setup

Pathlike = T.Union[str, Path]


def runner(pr: T.Dict[str, T.Any]) -> int:

    config_file = get_config_filename(pr["config_file"])
    # load configuration to know what directories to check
    p = read_config(config_file)
    for k in ("indat_size", "indat_grid", "indat_file"):
        f = p[k].expanduser().resolve()
        if pr["force"] or not f.is_file():
            ok = initialize_simulation(p, pr)
            if not ok:
                raise RuntimeError("\ncould not initialize simulation parameters. Is there a problem with config.nml?")
            if not f.is_file():
                raise FileNotFoundError(f"\ntried to initialize simulation but missing expected output file {f}")
            break

    if p.get("flagE0file") == 1:
        E0dir = p["E0dir"].expanduser().resolve()
        if not E0dir.is_dir():
            ok = initialize_simulation(p, pr)
            if not ok or not E0dir.is_dir():
                raise FileNotFoundError(E0dir)

    if p.get("flagprecfile") == 1:
        precdir = p["precdir"].expanduser().resolve()
        if not precdir.is_dir():
            ok = initialize_simulation(p, pr)
            if not ok or not precdir.is_dir():
                raise FileNotFoundError(precdir)

    # build checks
    mpiexec = check_mpiexec(pr["mpiexec"])
    logging.info(f"Detected mpiexec: {' '.join(mpiexec)}")

    gemexe = check_gemini_exe(pr["gemexe"])
    logging.info(f"using gemini executable: {gemexe}")

    out_dir = check_outdir(pr["out_dir"])

    Nmpi = get_mpi_count(config_file)

    cmd = mpiexec + ["-n", str(Nmpi), str(gemexe), str(config_file), str(out_dir)]
    print(" ".join(cmd), "\n")

    batcher = hpc_batch_detect()
    if batcher:
        job_file = hpc_batch_create(batcher, out_dir, cmd)  # noqa: F841
        # hpc_submit_job(job_file)
        print("Please examine batch file", job_file, "and when ready submit the job as usual.")
        ret = 0
    else:
        ret = subprocess.run(cmd).returncode

    return ret


@functools.lru_cache()
def wsl_available() -> bool:
    """
    heuristic to detect if Windows Subsystem for Linux is available.

    Uses presence of /etc/os-release in the WSL image to say Linux is there.
    This is a de facto file standard across Linux distros.
    """

    has_wsl = False
    if os.name == "nt" and shutil.which("wsl"):
        has_wsl = wsl_file_exist("/etc/os-release")

    return has_wsl


def wsl_file_exist(file: Pathlike) -> bool:
    """
    path is specified as if in WSL
    NOT //wsl$/Ubuntu/etc/os-release
    but /etc/os-release
    """
    if os.name != "nt":
        return False

    try:
        return subprocess.run(["wsl", "test", "-f", str(file)], timeout=10).returncode == 0
    except subprocess.TimeoutExpired:
        return False


def check_compiler():

    fc = os.environ.get("FC")
    fc = shutil.which(fc) if fc else shutil.which("gfortran")
    if not fc:
        raise EnvironmentError("Cannot find Fortran compiler e.g. Gfortran")


def initialize_simulation(p: T.Dict[str, T.Any], pr: T.Dict[str, T.Any]) -> bool:
    """
    TODO: these functions will be in Python

    GNU Octave doesn't have the HDF5 API needed for Gemini
    """

    env = os.environ
    if not os.environ.get("GEMINI_ROOT"):
        env["GEMINI_ROOT"] = Path(__file__).parents[1].as_posix()
    matlab_script_dir = Path(env["GEMINI_ROOT"], "setup")

    if pr["matlab"]:
        if not shutil.which("matlab"):
            raise EnvironmentError("Matlab not found")

        check_compiler()

        cmd = ["matlab", "-batch", f"model_setup('{p['nml']}')"]
        print("Initializing simulation: \n", " ".join(cmd))

        ret = subprocess.run(cmd, cwd=matlab_script_dir, env=env)
        return ret.returncode == 0
    else:
        model_setup(p["nml"], pr["out_dir"])
        return True


def check_mpiexec(mpiexec: Pathlike) -> T.List[str]:
    """ check that specified mpiexec exists on this system """

    if not mpiexec:
        mpiexec = "mpiexec"
    mpi_root = os.getenv("MPI_ROOT")
    if mpi_root:
        mpi_root += "/bin"
    mpiexec = shutil.which(mpiexec, path=mpi_root)

    if mpiexec:
        mpi_exec = [mpiexec]
    else:
        msg = "Need mpiexec to run simulations"
        if os.name == "nt":
            msg += "\n\n Typically Windows users use Windows Subsystem for Linux (WSL)"
            if wsl_available():
                msg += "\n 😊  WSL appears to be already installed on your PC, look in the Start menu for Ubuntu or see:"
                mpi_exec = ["wsl", "mpiexec"]
                ret = subprocess.run(mpi_exec + ["--version"])
                if ret.returncode != 0:
                    mpi_exec = []
            msg += "\n 📖  WSL install guide: https://docs.microsoft.com/en-us/windows/wsl/install-win10"

    if not mpi_exec:
        raise EnvironmentError(msg)

    return mpi_exec


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
            raise EnvironmentError(f"\nCannot find gemini.bin under {build_dir}")

    gemexe = str(Path(gemexe).resolve())

    ret = subprocess.run(gemexe, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, universal_newlines=True)
    if ret.returncode != 77:
        raise RuntimeError(
            f"\n{gemexe} was not runnable on your platform. Try recompiling on this computer type."
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
