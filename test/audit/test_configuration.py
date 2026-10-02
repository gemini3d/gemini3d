"""Check isolation and environment composition in the generated simulation tests."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import tempfile


def properties(test):
    return {prop["name"]: prop["value"] for prop in test.get("properties", [])}


def check_configuration(ctest, build, mpi_tmpdir, windows=False):
    result = subprocess.run(
        [ctest, "--test-dir", str(build), "--show-only=json-v1"],
        check=True, capture_output=True, text=True, timeout=30,
    )
    tests = {test["name"]: test for test in json.loads(result.stdout)["tests"]}
    frontends = [name for name in tests if name.startswith("gemini_run:")]
    assert frontends, "No simulation frontend tests registered"
    for frontend_name in frontends:
        case = frontend_name.split(":")[1]
        frontend = tests[frontend_name]
        native = tests[f"gemini:{case}"]
        dryrun = tests[f"gemini:{case}:dryrun"]
        copy = tests[f"{case}:copy_data_frontend"]
        download = tests[f"{case}:download"]
        for test in (copy, download):
            command = test["command"]
            assert Path(command[command.index("-P") + 1]).is_file(), test["name"]
        copy_args = copy["command"]
        frontend_dir = next(arg.removeprefix("-Doutdir:PATH=") for arg in copy_args
                            if arg.startswith("-Doutdir:PATH="))
        native_dir = next(arg.removeprefix("-Doutdir:PATH=") for arg in download["command"]
                          if arg.startswith("-Doutdir:PATH="))
        assert Path(frontend_dir) != Path(native_dir), case
        assert frontend_dir in frontend["command"], case
        assert native_dir in native["command"] and native_dir in dryrun["command"], case
        assert properties(copy)["FIXTURES_REQUIRED"] == [f"{case}:download_fxt"], case
        assert properties(copy)["FIXTURES_SETUP"] == [f"{case}:frontend_copy_fxt"], case
        assert f"{case}:frontend_copy_fxt" in properties(frontend)["FIXTURES_REQUIRED"], case
        assert set(properties(dryrun)["FIXTURES_REQUIRED"]) == {
            f"{case}:download_fxt", "gemini_exe_fxt",
        }, case
        assert properties(native)["FIXTURES_REQUIRED"] == [f"{case}:dryrun"], case

        for test in (frontend, dryrun, native):
            props = properties(test)
            workdir = props["WORKING_DIRECTORY"]
            environment = props["ENVIRONMENT_MODIFICATION"]
            assert f"HWMPATH=set:{workdir}" in environment, test["name"]
            if mpi_tmpdir:
                assert f"TMPDIR=set:{mpi_tmpdir}" in environment, test["name"]
            if windows:
                assert any(value.startswith("PATH=path_list_prepend:")
                           for value in environment), test["name"]
            assert props["RESOURCE_LOCK"] == ["cpu_mpi"], test["name"]
        print(f"{case}: isolated frontend/native directories and preserved test environment")


def check_fixture(cmake, ctest):
    with tempfile.TemporaryDirectory(prefix="gemini test configuration ") as work:
        root = Path(work)
        source = Path(__file__).resolve().parent / "configuration"
        for windows in (False, True):
            build = root / ("windows" if windows else "native")
            mpi_tmpdir = root / "mpi tmp"
            subprocess.run(
                [cmake, "-S", str(source), "-B", str(build),
                 f"-Dmpi_tmpdir:PATH={mpi_tmpdir}",
                 f"-Dtest_windows:BOOL={'ON' if windows else 'OFF'}"],
                check=True, capture_output=True, text=True, timeout=30,
            )
            check_configuration(ctest, build, str(mpi_tmpdir), windows=windows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--ctest", default="ctest")
    parser.add_argument("--cmake", default="cmake")
    parser.add_argument("--build", type=Path)
    parser.add_argument("--mpi-tmpdir", default="")
    args = parser.parse_args()
    if args.build:
        check_configuration(args.ctest, args.build, args.mpi_tmpdir, windows=os.name == "nt")
    else:
        check_fixture(args.cmake, args.ctest)
