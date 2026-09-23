"""Exercise optional/required discovery and stale-cache handling without native downloads."""
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]


def find_cmake():
    for candidate in (os.environ.get("CMAKE"), shutil.which("cmake"), "cmake"):
        if not candidate:
            continue
        resolved = shutil.which(candidate) or candidate
        if resolved == "cmake" or Path(resolved).is_file():
            return resolved
    return "cmake"


class QualificationDependencies(unittest.TestCase):
    def setUp(self):
        self.work = tempfile.TemporaryDirectory()
        self.addCleanup(self.work.cleanup)
        self.root = Path(self.work.name)
        self.source = self.root / "source"
        self.source.mkdir()
        (self.source / "cmake").mkdir()
        shutil.copy2(ROOT / "cmake/cpu_count.cmake", self.source / "cmake/cpu_count.cmake")
        (self.source / "CMakeLists.txt").write_text(
            "cmake_minimum_required(VERSION 3.25)\n"
            "project(qualification_dependency_probe LANGUAGES NONE)\n"
            "set(gemini3d_IS_TOP_LEVEL TRUE)\n"
            'set(BUILD_TESTING ON CACHE BOOL "Test gate")\n'
            f'include("{ROOT.as_posix()}/options.cmake")\n'
            f'include("{ROOT.as_posix()}/cmake/python.cmake")\n'
            'file(WRITE "${CMAKE_BINARY_DIR}/found.txt" '
            '"${Python_Interpreter_FOUND};${NUMPY_FOUND};${H5PY_FOUND};${SCIPY_FOUND}")\n'
        )
        # Deterministic import probes work even on an optional, stdlib-only host.
        self.modules = self.root / "modules"
        self.modules.mkdir()
        for module in ("numpy", "h5py", "scipy"):
            (self.modules / (module + ".py")).write_text("__version__ = 'test-fixture'\n")
        self.env = dict(os.environ, PYTHONPATH=str(self.modules))
        cmake = os.environ.get("CMAKE")
        if cmake and Path(cmake).is_file():
            self.env["PATH"] = str(Path(cmake).parent) + os.pathsep + self.env.get("PATH", "")

    def configure(self, *args):
        return subprocess.run(
            [find_cmake(), "-S", str(self.source), "-B", str(self.root / "build"),
             f"-DPython_EXECUTABLE={sys.executable}", *args],
            env=self.env, capture_output=True, text=True, timeout=60
        )

    def found(self):
        return (self.root / "build/found.txt").read_text().split(";")

    def test_successful_required_configuration(self):
        result = self.configure("-Dgemini3d_require_qualification=ON")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertTrue(all(value.upper() == "TRUE" for value in self.found()))

    def test_each_failed_import_refreshes_cache_and_fails_required_mode(self):
        for module in ("numpy", "h5py", "scipy"):
            with self.subTest(module=module):
                path = self.modules / (module + ".py")
                result = self.configure("-Dgemini3d_require_qualification=OFF",
                                        "-DPYGEMINI_FOUND=TRUE")
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                self.assertTrue(all(value.upper() == "TRUE" for value in self.found()))
                path.write_text("raise ImportError('missing qualification dependency fixture')\n")
                result = self.configure()
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                self.assertEqual(self.found()[("numpy", "h5py", "scipy").index(module) + 1], "FALSE")
                result = self.configure("-Dgemini3d_require_qualification=ON")
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("Required qualification dependencies unavailable", result.stderr)
                self.assertIn(module, result.stderr)
                path.write_text("__version__ = 'test-fixture-restored'\n")

    def test_missing_interpreter_clears_cached_imports(self):
        self.assertEqual(self.configure().returncode, 0)
        missing = f"-DPython_EXECUTABLE={self.root / 'missing-python'}"
        result = self.configure(missing)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(self.found()[1:], ["FALSE"] * 3)
        result = self.configure(missing, "-Dgemini3d_require_qualification=ON")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Python >= 3.11", result.stderr)

    @unittest.skipIf(os.name == "nt", "POSIX fake interpreter launcher")
    def test_python_before_file_digest_is_not_accepted(self):
        executable = self.root / "old-python"
        executable.write_text(
            f"#!{sys.executable}\n"
            "import os,sys\n"
            "if len(sys.argv)>2 and sys.argv[1]=='-c' and 'version_info' in sys.argv[2]:\n"
            "    sys.version_info = (3,10,99,'final',0)\n"
            "    exec(sys.argv[2])\n"
            "else:\n"
            f"    os.execv({sys.executable!r},[{sys.executable!r},*sys.argv[1:]])\n"
        )
        executable.chmod(0o755)
        old = f"-DPython_EXECUTABLE={executable}"
        result = self.configure(old)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(self.found()[1:], ["FALSE"] * 3)
        result = self.configure(old, "-Dgemini3d_require_qualification=ON")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Python >= 3.11", result.stderr)

    def test_required_qualification_cannot_disable_tests(self):
        for flag in ("BUILD_TESTING", "gemini3d_BUILD_TESTING"):
            with self.subTest(flag=flag):
                result = self.configure("-Dgemini3d_require_qualification=ON",
                                        "-DBUILD_TESTING=ON", "-Dgemini3d_BUILD_TESTING=ON", f"-D{flag}=OFF")
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("requires BUILD_TESTING=ON", result.stderr)

    def test_unusable_cmake_environment_falls_back_to_path(self):
        cmake = os.environ.get("CMAKE")
        self.addCleanup(lambda: os.environ.pop("CMAKE", None) if cmake is None else os.environ.__setitem__("CMAKE", cmake))
        os.environ["CMAKE"] = "D:/not-valid-under-wsl/cmake.exe"
        result = self.configure("-Dgemini3d_require_qualification=ON")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertTrue(all(value.upper() == "TRUE" for value in self.found()))


if __name__ == "__main__":
    unittest.main()
