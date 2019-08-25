import subprocess
import shutil

git = shutil.which("git")


def gitrev() -> str:
    if not git:
        return ""

    return subprocess.check_output(
        [git, "rev-parse", "--short", "HEAD"], universal_newlines=True
    ).strip()
