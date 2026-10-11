"""Check trace command rewriting for native GNU and Intel Ninja commands."""
from pathlib import Path
import runpy
import shlex
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[2]


class TraceCompilation(unittest.TestCase):
    def test_native_compiler_flags(self):
        for compiler, flags, source_name in (
            ("gfortran", ["-fpreprocessed", "-Jinclude", "-fcheck=all"],
             "PDEparabolic.f90-pp.f90"),
            ("gfortran", ["-fpreprocessed", "-J", "include", "-fcheck=all"],
             "PDEparabolic.f90-pp.f90"),
            ("/opt/intel/bin/ifx", ["-module", "include", "-CB"],
             "PDEparabolic.f90-pp.f90"),
            ("ifort", ["-module=include", "-CB"], "PDEparabolic.f90"),
        ):
            with self.subTest(compiler=compiler, flags=flags), tempfile.TemporaryDirectory() as tmp:
                build = Path(tmp)
                source = build / "PDEparabolic.f90"
                source.write_text("call check_solver_info(info,'TRBDF21D BDF2 stage')\n"
                                  "end module PDEparabolic\n")
                probe = build / "native/physics_native.o"
                probe.parent.mkdir()
                work = build / "trace"
                compile_command = shlex.join(
                    [compiler, "-Iinclude", *flags, "-c", source_name, "-o", "original.o"])
                commands = (
                    "cmake -E cmake_ninja_depends --ddi=original.o.ddi\n"
                    + compile_command + "\n"
                    + shlex.join([compiler, "-c", "other.f90", "-o", "other.o"]) + "\n")
                link_command = (
                    ": && ifx continuity_floor.f90.o libgemini3d.a -o original && :\n")
                argv = ["trace_diffusion_stages.py", "--source", str(source),
                        "--build", str(build), "--probe", str(probe), "--work", str(work)]
                with patch.object(sys, "argv", argv), \
                        patch("subprocess.check_output", side_effect=[commands, link_command]), \
                        patch("subprocess.run", return_value=subprocess.CompletedProcess(
                            [], 0, stdout="", stderr="")) as run:
                    runpy.run_path(str(ROOT / "scripts/qualification/trace_diffusion_stages.py"),
                                   run_name="__main__")
                command = run.call_args_list[0].args[0]
                self.assertEqual(command[0], compiler)
                self.assertEqual(command[command.index("-c") + 1],
                                 str(work / "PDEparabolic_trace.f90"))
                self.assertEqual(command[command.index("-o") + 1],
                                 str(work / "PDEparabolic_trace.o"))
                self.assertIn("-Iinclude", command)
                self.assertNotIn("-fpreprocessed", command)
                self.assertNotIn("include", command)
                if compiler == "gfortran":
                    self.assertIn("-fcheck=all", command)
                    self.assertIn("-J" + str(work), command)
                    self.assertNotIn("-Jinclude", command)
                    self.assertIn("-ffree-line-length-none", command)
                else:
                    self.assertIn("-CB", command)
                    self.assertEqual(command[command.index("-module") + 1], str(work))
                    self.assertNotIn("-module=include", command)
                    self.assertNotIn("-ffree-line-length-none", command)
                    self.assertFalse(any(t.startswith("-J") for t in command))


if __name__ == "__main__":
    unittest.main()
