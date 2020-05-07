import typing as T
import os
import logging
import subprocess
import shutil
from pathlib import Path

from .mpi import get_mpi_count
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
            model_setup(p["nml"], pr["out_dir"])
            if not f.is_file():
                raise FileNotFoundError(f"\ntried to initialize simulation but missing expected output file {f}")
            break

    if p.get("flagE0file") == 1:
        E0dir = p["E0dir"].expanduser().resolve()
        if not E0dir.is_dir():
            model_setup(p["nml"], pr["out_dir"])
            if not E0dir.is_dir():
                raise FileNotFoundError(E0dir)

    if p.get("flagprecfile") == 1:
        precdir = p["precdir"].expanduser().resolve()
        if not precdir.is_dir():
            model_setup(p["nml"], pr["out_dir"])
            if not precdir.is_dir():
                raise FileNotFoundError(precdir)

    # build checks
    mpiexec = check_mpiexec(pr["mpiexec"])
    logging.info(f"Detected mpiexec: {' '.join(mpiexec)}")

    gemexe = check_gemini_exe(pr["gemexe"])
    logging.info(f"using gemini executable: {gemexe}")

    out_dir = check_outdir(pr["out_dir"])

    Nmpi = get_mpi_count(config_file, pr["cpu_count"])

    cmd = mpiexec + ["-n", str(Nmpi), str(gemexe), str(config_file), str(out_dir)]
    if pr["out_format"]:
        cmd += ["-out_format", pr["out_format"]]

    # %% attempt dry run, but don't fail in case intended for HPC
    logging.info('Attempting Gemini dry run of first time step')
    dryrun_ok = subprocess.run(cmd + ["-dryrun"], stdout=subprocess.DEVNULL).returncode == 0
    if dryrun_ok:
        logging.info("OK: Gemini dry run")
    else:
        logging.error("Gemini dry run failed.")

    batcher = hpc_batch_detect()
    if batcher:
        job_file = hpc_batch_create(batcher, out_dir, cmd)  # noqa: F841
        # hpc_submit_job(job_file)
        print("Please examine batch file", job_file, "and when ready submit the job as usual.")
        ret = 0
    elif dryrun_ok:
        print("\nBEGIN Gemini run with command:")
        print(" ".join(cmd), "\n")
        ret = subprocess.run(cmd).returncode
    else:
        ret = -1
        logging.error("Did not run Gemini due to dry run failure")

    return ret


def check_compiler():

    fc = os.environ.get("FC")
    fc = shutil.which(fc) if fc else shutil.which("gfortran")
    if not fc:
        raise EnvironmentError("Cannot find Fortran compiler e.g. Gfortran")


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
        raise EnvironmentError("Need mpiexec to run simulations")

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
        build_dir = Path(__file__).resolve().parents[2] / "build"
        if not build_dir.is_dir():
            raise EnvironmentError(f"GEMINI build directory missing: {build_dir}")

        for d in (build_dir, build_dir / "Release", build_dir / "Debug"):
            gemexe = shutil.which("gemini.bin", path=str(d))
            if gemexe:
                break
        if not gemexe:
            raise EnvironmentError(f"\nCannot find gemini.bin under {build_dir}")

    gemexe = str(Path(gemexe).resolve())

    ret = subprocess.run(gemexe, stdout=subprocess.DEVNULL)
    if ret.returncode != 77:
        raise RuntimeError(
            f"\n{gemexe} was not runnable on your platform. Try recompiling on this computer type."
            "E.g. different HPC nodes may not have the CPU feature sets."
        )

    return gemexe


def check_outdir(out_dir: Pathlike) -> Path:

    out_dir = Path(out_dir).expanduser().resolve()
    if out_dir.is_file():
        raise NotADirectoryError(f"please specify output DIRECTORY, you specified {out_dir}")
    if not out_dir.is_dir():
        out_dir.mkdir(parents=True, exist_ok=True)

    return out_dir
