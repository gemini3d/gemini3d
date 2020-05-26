import os
import subprocess
import binascii
from pathlib import Path
import typing as T
import shutil

R = Path(__file__).parent


def hpc_submit_job(batcher: str, job_file: Path):
    """
    submits batch job
    """

    if batcher == "qsub":
        subprocess.run(["qsub", str(job_file)])
    else:
        raise LookupError(f"batcher {batcher} not yet listed. Please raise a Gemini Github Issue.")


def hpc_batch_create(batcher: str, out_dir: Path, cmd: T.Sequence[str]) -> Path:
    """
    creates HPC batch scripts for known systems

    assumes that user-specific parameters like account number are already set
    as environment variables
    or static configuration files not handled by this scripts.

    This function assumes a script template exists, and it merely appends lines
    to the end of that template.

    TODO:

    1. determine estimated wallclock time to request on HPC
    2. determine requested HPC RAM per node (limit is master node)
    3. format number of nodes request
    """

    template_dir = R / "templates"
    Nchar = 6  # arbitrary number of characters

    if batcher == "qsub":
        template_file = template_dir / "qsub_template.sh"
        template = template_file.read_text()
        job_file = out_dir / f"job_{binascii.b2a_hex(os.urandom(Nchar)).decode('ascii')}.sh"
        print("writing job file", job_file)
        text = template + "\n" + " ".join(cmd)
        job_file.write_text(text)
    else:
        raise LookupError(f"batcher {batcher} not yet listed. Please raise a Gemini Github Issue.")

    return job_file


def hpc_batch_detect() -> str:
    """
    Assuming a known job batching system, we will create a template for the user
    to verify and then the user will run.
    """

    batcher = None

    if shutil.which("qsub"):
        batcher = "qsub"

    return batcher
