"""Local installer orchestration and fail-closed run checks (no package downloads)."""
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
from unittest.mock import patch

SPEC = importlib.util.spec_from_file_location(
    "local_environment", Path(__file__).resolve().parents[2] / "scripts/local_environment.py")
local = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(local)


class LocalEnvironment(unittest.TestCase):
    @unittest.skipUnless(shutil.which("ctest") and shutil.which("cmake"), "CTest is required")
    def test_real_ctest_execution_receipt_and_skipped_test_rejection(self):
        with tempfile.TemporaryDirectory(prefix="gemini ctest ") as directory:
            root = Path(directory)
            cmake = Path(shutil.which("cmake")).as_posix()
            tests = root / "CTestTestfile.cmake"
            tests.write_text(
                f'add_test(passes "{cmake}" -E true)\n'
                f'add_test(optional "{cmake}" -E false)\n'
                'set_tests_properties(optional PROPERTIES DISABLED TRUE)\n')
            result = local.verified_tests(shutil.which("ctest"), root, None, 1, "fixture",
                                          required={"passes"})
            self.assertEqual(result["passed"], ["passes"])
            self.assertEqual(result["disabled"], ["optional"])
            self.assertEqual(result["junit_sha256"], local.digest(root / result["junit"]))
            with self.assertRaisesRegex(RuntimeError, "required tests missing or disabled"):
                local.verified_tests(shutil.which("ctest"), root, None, 1, "fixture", required={"optional"})
            tests.write_text(f'add_test(skipped "{cmake}" -E true)\n'
                             'set_tests_properties(skipped PROPERTIES SKIP_RETURN_CODE 0)\n')
            with self.assertRaisesRegex(RuntimeError, "executed test results"):
                local.verified_tests(shutil.which("ctest"), root, None, 1, "fixture")
            tests.write_text(f'add_test(failed "{cmake}" -E false)\n')
            with self.assertRaises(subprocess.CalledProcessError):
                local.verified_tests(shutil.which("ctest"), root, None, 1, "fixture")

    def test_test_receipt_rejects_missing_disabled_duplicate_or_stale_evidence(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for tests in ([], [{"name": "test", "properties": [{"name": "DISABLED", "value": True}]}],
                          [{"name": "test"}, {"name": "test"}]):
                with self.subTest(tests=tests), patch.object(
                        local, "execute", return_value=json.dumps({"tests": tests})) as run:
                    with self.assertRaisesRegex(RuntimeError, "no enabled tests or duplicate"):
                        local.verified_tests("ctest", root, None, 1, "fixture")
                    self.assertEqual(run.call_count, 1)
            report = root / "local-fixture.xml"
            passed = '<testsuite><testcase name="test" status="run"/></testsuite>'
            reports = [None, '<testsuite/>', passed.replace('name="test"', 'name="other"'),
                       passed.replace('</testsuite>', '<testcase name="test" status="run"/></testsuite>'),
                       passed.replace('/>', '><skipped/></testcase>'),
                       passed.replace('/>', '><failure/></testcase>'),
                       passed.replace('/>', '><error/></testcase>'),
                       passed.replace('status="run"', 'status="notrun"')]
            for contents in reports:
                with self.subTest(contents=contents):
                    report.write_text(passed)
                    def execute(command, **kwargs):
                        if "--show-only=json-v1" in command:
                            return '{"tests":[{"name":"test"}]}'
                        self.assertFalse(report.exists(), "Old test results must be removed before CTest runs")
                        if contents is not None:
                            report.write_text(contents)
                    with patch.object(local, "execute", side_effect=execute):
                        with self.assertRaises((RuntimeError, FileNotFoundError)):
                            local.verified_tests("ctest", root, None, 1, "fixture")

    @unittest.skipUnless(shutil.which("pwsh"), "PowerShell is required for WSL wrapper tests")
    def test_wsl_install_options_and_exit_status(self):
        script = str(local.SOURCE / "scripts/install-local.ps1").replace("'", "''")
        command = """
function wsl.exe {
    $global:LASTEXITCODE = 0
    if ($args[3] -eq 'wslpath') {
        '/checkout with spaces/scripts/install-local.sh'
    } else {
        ConvertTo-Json -Compress -InputObject @($args)
        $global:LASTEXITCODE = INSTALL_EXIT
    }
}
& 'SCRIPT' -Distribution Ubuntu-24.04 -SystemDeps -Root '/local env' `
    -Jobs 3 -BuildType Debug -ReferenceTests -SourceCache '/cache dir/sources.cmake'
""".replace("SCRIPT", script)
        run = subprocess.run(["pwsh", "-NoProfile", "-NonInteractive", "-Command",
                              command.replace("INSTALL_EXIT", "0")],
                             capture_output=True, text=True)
        self.assertEqual(run.returncode, 0, run.stderr)
        self.assertEqual(json.loads(run.stdout), [
            "--distribution", "Ubuntu-24.04", "--exec", "bash",
            "/checkout with spaces/scripts/install-local.sh",
            "--jobs", "3", "--build-type", "Debug", "--system-deps", "--root", "/local env",
            "--reference-tests", "--source-cache", "/cache dir/sources.cmake",
        ])
        failed = subprocess.run(["pwsh", "-NoProfile", "-NonInteractive", "-Command",
                                 command.replace("INSTALL_EXIT", "17")],
                                capture_output=True, text=True)
        self.assertNotEqual(failed.returncode, 0)
        self.assertIn("Local installation failed (exit 17)", failed.stderr)
        for old, new in (("-Jobs 3", "-Jobs 0"), ("-BuildType Debug", "-BuildType Invalid")):
            invalid = subprocess.run(["pwsh", "-NoProfile", "-NonInteractive", "-Command",
                                      command.replace("INSTALL_EXIT", "0").replace(old, new)],
                                     capture_output=True, text=True)
            self.assertNotEqual(invalid.returncode, 0)
            self.assertEqual(invalid.stdout, "")

    def test_run_arguments_after_action_and_separator(self):
        with tempfile.TemporaryDirectory(prefix="gemini local ") as directory:
            root = Path(directory)
            exe = root / "install/bin/gemini.bin"
            with patch.object(local, "check", return_value=(exe, ["/mpi/mpiexec", "-n"], {})), \
                    patch.object(local, "execute") as run:
                local.main(["run", "--root", directory, "--case", directory, "--ranks", "2",
                            "--", "-dryrun", "-manual_grid", "1", "2"])
            self.assertEqual(run.call_args.args[0],
                             ["/mpi/mpiexec", "-n", 2, exe, root, "-dryrun", "-manual_grid", "1", "2"])
            self.assertEqual(run.call_args.kwargs["cwd"], exe.parent)

    def test_check_requires_success_record_and_never_installs(self):
        with tempfile.TemporaryDirectory() as directory, patch.object(local, "execute") as run:
            with self.assertRaisesRegex(RuntimeError, "No completed installation"):
                local.check(Path(directory))
            run.assert_not_called()

    def test_debug_executable(self):
        self.assertEqual(local.executable(Path("/local"), "Debug"), Path("/local/bin/gemini.bin.debug"))
        with self.assertRaises(RuntimeError):
            local.executable(Path("/local"), "../bad")

    def test_missing_native_prerequisites(self):
        with patch.object(local.shutil, "which", return_value=None):
            with self.assertRaisesRegex(RuntimeError, "Missing native tools"):
                local.native_prerequisites({"PATH": ""})

    def test_failed_update_invalidates_old_success(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "venv").mkdir()
            (root / "venv/pyvenv.cfg").write_text("fixture")
            record = root / "environment.json"
            record.write_text("{}")
            with patch.object(local, "native_prerequisites"), patch.object(
                    local, "execute", side_effect=subprocess.CalledProcessError(1, ["pip"])):
                with self.assertRaises(SystemExit):
                    local.main(["install", "--root", directory])
            self.assertFalse(record.exists())

    @unittest.skipUnless(shutil.which("bash"), "POSIX installer uses bash")
    def test_failed_system_update_invalidates_old_success(self):
        with tempfile.TemporaryDirectory(prefix="gemini setup ") as directory:
            root = Path(directory)
            shutil.copyfile(local.SOURCE / "scripts/install-local.sh", root / "install-local.sh")
            (root / "install-system-deps.sh").write_text("exit 17\n")
            record = root / "environment.json"
            for option in (["--root", directory], [f"--root={directory}"]):
                record.write_text("{}")
                run = subprocess.run(["bash", str(root / "install-local.sh"), "--system-deps", *option])
                self.assertEqual(run.returncode, 17)
                self.assertFalse(record.exists())

    def test_changed_mpi_launcher_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            exe = root / "install/bin/gemini.bin"
            exe.parent.mkdir(parents=True)
            exe.write_text("fixture")
            old, new = root / "old-mpi", root / "new-mpi"
            old.write_text("OpenMPI")
            new.write_text("MPICH")
            (root / "environment.json").write_text(json.dumps({
                "schema": "gemini.local-environment.1", "build_type": "Release",
                "requirements_sha256": local.digest(local.REQUIREMENTS),
                "executable_sha256": local.digest(exe), "mpi_launcher": str(old),
                "native_configuration": {}, "model_resources_sha256": {},
            }))
            with patch.object(local.shutil, "which", return_value=str(new)), \
                    patch.object(local, "execute") as run:
                with self.assertRaisesRegex(RuntimeError, "MPI launcher changed"):
                    local.check(root)
                run.assert_not_called()

    def test_changed_executable_rejected_before_running(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            exe = root / "install/bin/gemini.bin"
            exe.parent.mkdir(parents=True)
            exe.write_text("different binary")
            (root / "environment.json").write_text(json.dumps({
                "schema": "gemini.local-environment.1", "build_type": "Release",
                "requirements_sha256": local.digest(local.REQUIREMENTS),
                "executable_sha256": "0" * 64,
            }))
            with patch.object(local, "execute") as run:
                with self.assertRaisesRegex(RuntimeError, "changed since verification"):
                    local.check(root)
                run.assert_not_called()

    def test_success_checks_before_publication(self):
        self.check_successful_install(reference_tests=False)

    def test_successful_reference_install_records_complete_execution(self):
        self.check_successful_install(reference_tests=True)

    def check_successful_install(self, reference_tests):
        with tempfile.TemporaryDirectory(prefix="gemini install ") as directory:
            root = Path(directory)
            (root / "venv").mkdir()
            (root / "venv/pyvenv.cfg").write_text("fixture")
            (root / "build").mkdir()
            launcher = root / "mpiexec"
            launcher.write_text("fixture")
            (root / "build/CMakeCache.txt").write_text(
                f"gemini3d_realbits:STRING=64\nMPIEXEC_EXECUTABLE:FILEPATH={launcher}\n"
                "MPIEXEC_NUMPROC_FLAG:STRING=-n\ngemini3d_msis2:BOOL=ON\n")
            source_cache = root / "sources.cmake"
            source_cache.write_text("# hash-verified offline source locations\n")
            exe = root / "install/bin/gemini.bin"
            exe.parent.mkdir(parents=True)
            exe.write_text("fixture")
            resource = exe.parent / "msis21.parm"
            resource.write_text("model parameter fixture")
            commands = []
            libraries = sorted(local.LIBRARY_TESTS)
            all_tests = sorted(local.LIBRARY_TESTS | local.REFERENCE_TESTS | {"audit:example"})

            def execute(command, **kwargs):
                commands.append([str(item) for item in command])
                self.assertFalse((root / "environment.json").exists())
                names = libraries if "-R" in command else ["audit:example"] if "-L" in command else all_tests
                if "--show-only=json-v1" in command:
                    return json.dumps({"tests": [{"name": name} for name in names]})
                if "--output-junit" in command:
                    report = command[command.index("--output-junit") + 1]
                    report.write_text("<testsuite>" + "".join(
                        f'<testcase name="{name}" status="run"/>' for name in names) + "</testsuite>")
                if "--format=json" in command:
                    return "[]"

            with patch.object(local, "native_prerequisites"), patch.object(local, "python_check"), \
                    patch.object(local, "execute", side_effect=execute), \
                    patch.object(local, "source_identity", return_value={}), \
                    patch.object(local, "source_inventory", return_value={}):
                local.main(["install", "--root", directory, "--jobs", "2",
                            "--source-cache", str(source_cache)] +
                           (["--reference-tests"] if reference_tests else []))
            saved = json.loads((root / "environment.json").read_text())
            self.assertEqual(saved["tests"], all_tests)
            if reference_tests:
                self.assertEqual(saved["verification"], "enabled-registered-tests")
                self.assertEqual(saved["test_results"]["reference"]["passed"], all_tests)
                self.assertEqual(set(saved["test_results"]), {"reference"})
            else:
                self.assertEqual(saved["verification"], "unit-tests-only")
                self.assertEqual(saved["test_results"]["libraries"]["passed"], libraries)
                self.assertEqual(saved["test_results"]["unit"]["passed"], ["audit:example"])
                self.assertTrue(any("-R" in c and "HDF5_standalone" in c[c.index("-R") + 1] for c in commands))
            for result in saved["test_results"].values():
                self.assertEqual(result["junit_sha256"], local.digest(root / "build" / result["junit"]))
            self.assertEqual(saved["model_resources_sha256"], {"msis21.parm": local.digest(resource)})
            self.assertTrue(any("-Dgemini3d_require_qualification=ON" in c for c in commands))
            self.assertTrue(any("-C" in c and str(source_cache) in c for c in commands))
            test_indices = [i for i, c in enumerate(commands) if "--no-tests=error" in c]
            install_index = next(i for i, c in enumerate(commands) if "--install" in c)
            self.assertEqual(len(test_indices), 1 if reference_tests else 2)
            self.assertTrue(all(index < install_index for index in test_indices))
            resource.write_text("changed model parameters")
            with patch.object(local, "execute") as run:
                with self.assertRaisesRegex(RuntimeError, "model resources changed"):
                    local.check(root)
                run.assert_not_called()

    def test_reference_install_requires_all_nine_comparisons_before_publication(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "venv").mkdir()
            (root / "venv/pyvenv.cfg").write_text("fixture")
            names = sorted(local.LIBRARY_TESTS | local.REFERENCE_TESTS)
            self.assertEqual(len(local.REFERENCE_TESTS), 9)
            for disabled in (False, True):
                tests = [{"name": name} for name in names]
                if disabled:
                    tests[-1]["properties"] = [{"name": "DISABLED", "value": True}]
                else:
                    tests.pop()
                def execute(command, **kwargs):
                    self.assertNotIn("--install", command)
                    self.assertNotIn("--output-junit", command)
                    if "--show-only=json-v1" in command:
                        return json.dumps({"tests": tests})
                with self.subTest(disabled=disabled), patch.object(local, "native_prerequisites"), \
                        patch.object(local, "python_check"), patch.object(local, "execute", side_effect=execute):
                    with self.assertRaises(SystemExit):
                        local.main(["install", "--root", directory, "--reference-tests"])
                self.assertFalse((root / "environment.json").exists())


if __name__ == "__main__":
    unittest.main()
