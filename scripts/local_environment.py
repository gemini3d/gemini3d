"""Install, verify, and run a local GEMINI environment without shell activation.

Native package installation is explicit; running a model never downloads packages
or elevates privileges. Windows users run this inside WSL.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import venv


SOURCE = Path(__file__).resolve().parents[1]
REQUIREMENTS = SOURCE / "scripts/requirements-local.txt"
MPI_ENVIRONMENT = ("LD_LIBRARY_PATH", "DYLD_LIBRARY_PATH", "LD_PRELOAD",
                   "DYLD_INSERT_LIBRARIES", "I_MPI_ROOT", "MPI_ROOT", "MPI_HOME")


def execute(command, *, env=None, capture=False, cwd=None):
    return subprocess.run(
        [str(arg) for arg in command], check=True, env=env, cwd=cwd,
        text=True, stdout=subprocess.PIPE if capture else None,
    ).stdout


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def paths(root):
    return root / "venv/bin", root / "build", root / "install"


def environment(root):
    bin_dir, _, prefix = paths(root)
    env = os.environ.copy()
    env["PATH"] = os.pathsep.join((str(bin_dir), str(prefix / "bin"), env.get("PATH", "")))
    env["VIRTUAL_ENV"] = str(root / "venv")
    env.pop("PYTHONHOME", None)
    env.pop("PYTHONPATH", None)
    env["PYTHONNOUSERSITE"] = "1"
    env["HDF5_PLUGIN_PRELOAD"] = "::"
    return env


def native_prerequisites(env):
    missing = [name for name in ("mpiexec", "mpifort", "cc", "c++")
               if not shutil.which(name, path=env["PATH"])]
    if missing:
        raise RuntimeError("Missing native tools: " + ", ".join(missing) +
                           ". Run scripts/install-local.sh --system-deps or load your HPC toolchain.")


def python_check(bin_dir, env):
    execute([bin_dir / "python", "-I", "-c",
             "import sys; assert sys.version_info >= (3, 11); "
             "import numpy, scipy, h5py; "
             "assert h5py.h5z.filter_avail(h5py.h5z.FILTER_DEFLATE)"], env=env)
    execute([bin_dir / "python", "-I", "-m", "pip", "check"], env=env)


def source_identity():
    try:
        revision = execute(["git", "-C", SOURCE, "rev-parse", "HEAD"], capture=True).strip()
        dirty = bool(execute(["git", "-C", SOURCE, "status", "--porcelain"], capture=True).strip())
        return {"revision": revision, "dirty": dirty}
    except (OSError, subprocess.CalledProcessError):
        return {"revision": None, "dirty": None}


def source_inventory():
    try:
        names = execute(["git", "-C", SOURCE, "ls-files", "-z", "--cached", "--others",
                         "--exclude-standard", "--", "src", "app", "include", "cmake",
                         "test", "scripts", "CMakeLists.txt", "CMakePresets.json", "options.cmake"],
                        capture=True)
        return {name: digest(SOURCE / name) for name in sorted(set(names.split("\0")))
                if name and (SOURCE / name).is_file() and not (SOURCE / name).is_symlink()}
    except (OSError, subprocess.CalledProcessError):
        return {}


def executable(prefix, build_type):
    if build_type not in ("Debug", "Release"):
        raise RuntimeError("Invalid recorded build type; rerun the installer.")
    return prefix / "bin" / ("gemini.bin.debug" if build_type == "Debug" else "gemini.bin")


def model_resources(cache):
    names = []
    if cache.get("gemini3d_msis2", "").upper() in ("ON", "TRUE", "1", "YES"):
        names.append("msis21.parm")
    if cache.get("gemini3d_hwm14", "").upper() in ("ON", "TRUE", "1", "YES"):
        names += ["hwm123114.bin", "dwm07b104i.dat", "gd2qd.dat"]
    return names


def install(args):
    root = args.root
    root.mkdir(parents=True, exist_ok=True)
    # Invalidate previous success before attempting an update.
    record = root / "environment.json"
    record.unlink(missing_ok=True)
    bin_dir, build, prefix = paths(root)
    if not (root / "venv/pyvenv.cfg").is_file():
        venv.EnvBuilder(with_pip=True).create(root / "venv")
    env = environment(root)
    native_prerequisites(env)
    execute([bin_dir / "python", "-I", "-m", "pip", "install", "-r", REQUIREMENTS], env=env)
    python_check(bin_dir, env)
    cmake = bin_dir / "cmake"
    configure = [cmake]
    if args.source_cache:
        configure += ["-C", args.source_cache.resolve()]
    execute([*configure, "-S", SOURCE, "-B", build, "-G", "Ninja",
             f"-DCMAKE_BUILD_TYPE={args.build_type}", f"-DCMAKE_INSTALL_PREFIX={prefix}",
             "-DCMAKE_INSTALL_LIBDIR=lib", f"-DPython_EXECUTABLE={bin_dir / 'python'}",
             "-Dgemini3d_require_qualification=ON", "-DBUILD_TESTING=ON",
             "-Dgemini3d_BUILD_TESTING=ON", "-Dgemini3d_test_simulations=ON",
             f"-DMPIEXEC_MAX_NUMPROCS={args.jobs}"], env=env)
    execute([cmake, "--build", build, "--parallel", args.jobs], env=env)
    ctest = bin_dir / "ctest"
    tests = [ctest, "--test-dir", build, "--output-on-failure", "--no-tests=error",
             "--parallel", args.jobs]
    if not args.reference_tests:
        execute([*tests, "-R", "^(HDF5_standalone_.*|GeminiMPIstandalone|GeminiMUMPSstandalone)$"], env=env)
        tests += ["-L", "unit"]
    execute(tests, env=env)
    execute([cmake, "--install", build], env=env)
    exe = executable(prefix, args.build_type)
    # Run from the installed resource directory, not the source/build directory.
    execute([exe, "-h"], env=env, cwd=prefix / "bin")
    inventory = json.loads(execute([ctest, "--test-dir", build, "--show-only=json-v1"],
                                  env=env, capture=True))
    cache = {}
    for line in (build / "CMakeCache.txt").read_text().splitlines():
        if "=" in line and not line.startswith(("#", "//")):
            key, value = line.split("=", 1)
            name = key.split(":", 1)[0]
            if name.startswith(("CMAKE_C_COMPILER", "CMAKE_CXX_COMPILER", "CMAKE_Fortran_COMPILER",
                                "MPI_", "MPIEXEC_", "HDF5_", "LAPACK_", "SCALAPACK_", "gemini3d_")):
                cache[name] = value
    launcher = Path(cache["MPIEXEC_EXECUTABLE"])
    manifest = {
        "schema": "gemini.local-environment.1", "source": source_identity(),
        "source_inventory_sha256": source_inventory(),
        "platform": platform.platform(), "python": sys.version,
        "build_type": args.build_type, "requirements_sha256": digest(REQUIREMENTS),
        "native_libraries": json.loads((SOURCE / "cmake/libraries.json").read_text()),
        "packages": json.loads(execute([bin_dir / "python", "-I", "-m", "pip", "list", "--format=json"],
                                      env=env, capture=True)),
        "native_configuration": cache,
        "model_resources_sha256": {name: digest(prefix / "bin" / name) for name in model_resources(cache)},
        "mpi_launcher": str(launcher), "mpi_launcher_sha256": digest(launcher),
        "mpi_environment": {name: env.get(name) for name in MPI_ENVIRONMENT},
        "tests": [test["name"] for test in inventory["tests"]],
        "verification": "all-registered-tests" if args.reference_tests else "unit-tests-only",
        "executable_sha256": digest(exe),
    }
    record.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Installed and checked: {prefix}\nEvidence: {record}")
    if not args.reference_tests:
        print("Unit tests passed; full reference/scientific qualification remains separate.")


def check(root):
    record = root / "environment.json"
    if not record.is_file():
        raise RuntimeError("No completed installation. Run scripts/install-local.sh first.")
    saved = json.loads(record.read_text())
    if saved.get("schema") != "gemini.local-environment.1":
        raise RuntimeError("Unsupported local environment record; rerun the installer.")
    bin_dir, _, prefix = paths(root)
    exe = executable(prefix, saved.get("build_type"))
    if saved.get("requirements_sha256") != digest(REQUIREMENTS) or saved.get("executable_sha256") != digest(exe):
        raise RuntimeError("Installation changed since verification; rerun the installer.")
    resources = {name: digest(prefix / "bin" / name)
                 for name in model_resources(saved["native_configuration"])}
    if resources != saved.get("model_resources_sha256"):
        raise RuntimeError("Installed model resources changed; rerun the installer.")
    env = environment(root)
    launcher = saved.get("mpi_launcher")
    current = shutil.which("mpiexec", path=env["PATH"])
    if not launcher or not current or not Path(current).samefile(launcher):
        raise RuntimeError("MPI launcher changed; load the installation's MPI environment or reinstall.")
    if digest(Path(launcher)) != saved.get("mpi_launcher_sha256") or \
            {name: env.get(name) for name in MPI_ENVIRONMENT} != saved.get("mpi_environment"):
        raise RuntimeError("MPI runtime environment changed since verification; rerun the installer.")
    python_check(bin_dir, env)
    packages = json.loads(execute([bin_dir / "python", "-I", "-m", "pip", "list", "--format=json"],
                                 env=env, capture=True))
    if packages != saved.get("packages"):
        raise RuntimeError("Python packages changed since verification; rerun the installer.")
    execute([exe, "-h"], env=env, cwd=prefix / "bin")
    return exe, [launcher, saved["native_configuration"]["MPIEXEC_NUMPROC_FLAG"]], env


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("install", "check", "run"))
    parser.add_argument("--root", type=Path, default=SOURCE / "build/local")
    parser.add_argument("--jobs", type=int, default=2)
    parser.add_argument("--build-type", choices=("Debug", "Release"), default="Release")
    parser.add_argument("--reference-tests", action="store_true")
    parser.add_argument("--source-cache", type=Path,
                        help="sources.cmake produced by scripts/offline_libraries.cmake")
    parser.add_argument("--system-deps", action="store_true",
                        help="Handled only by install-local.sh before Python is started")
    parser.add_argument("--ranks", type=int, default=1)
    parser.add_argument("--case", type=Path)
    arguments = list(sys.argv[1:] if argv is None else argv)
    extra = []
    if "--" in arguments:
        index = arguments.index("--")
        arguments, extra = arguments[:index], arguments[index + 1:]
    args = parser.parse_args(arguments)
    args.root = args.root.expanduser().resolve()
    if args.jobs < 1 or args.ranks < 1:
        parser.error("--jobs and --ranks must be positive")
    if args.source_cache and not args.source_cache.is_file():
        parser.error("--source-cache must be an existing CMake cache initialization file")
    if args.action == "run" and (args.case is None or not args.case.is_dir()):
        parser.error("run requires an existing --case directory")
    if args.action != "run" and extra:
        parser.error("Solver arguments after -- are only supported for run")
    if sys.version_info < (3, 11):
        parser.error("Python 3.11+ is required")
    if os.name == "nt":
        parser.error("Use install-local.ps1 and run this script inside Windows WSL")
    try:
        if args.action == "install":
            install(args)
        else:
            exe, mpi, env = check(args.root)
            if args.action == "run":
                execute([*mpi, args.ranks, exe, args.case.resolve(), *extra],
                        env=env, cwd=exe.parent)
    except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError) as error:
        parser.exit(1, f"Local environment error: {error}\n")


if __name__ == "__main__":
    main()
