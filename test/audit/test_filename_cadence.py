"""Both configuration formats must reject filename cadences below 10 ms."""
import argparse
from pathlib import Path
import subprocess

from run_native_probes import BASE


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--exe", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    cases = 0
    for cadence, accepted in (("0.001", False), ("0.009999", False), ("0.01", True),
                               ("0.015", True), ("60", True)):
        groups = {
            "dtout": BASE.replace("dtout=60", f"dtout={cadence}"),
            "dtglowout": BASE + f"\n&glow dtglow=1,dtglowout={cadence} /\n",
            "dtprec": BASE + f"\n&precip dtprec={cadence},prec_dir='p' /\n",
            "dtE0": BASE + f"\n&efield dtE0={cadence},E0_dir='e' /\n",
            "dtneu": BASE + f"\n&neutral_perturb dtneu={cadence},source_dir='n' /\n",
            "dtsolflux": BASE + f"\n&solflux dtsolflux={cadence},solfluxdir='s' /\n",
            "dtneuBGfile": BASE + f"\n&neutralBG_file dtneuBGfile={cadence},neutralBGdir='n' /\n",
        }
        # Deprecated INI still reaches the same cadence guard.
        ini = ("20 2 2013\n18000\n300\n{dtout}\n100 100 4\n0.5\n1000\n1\n0\n1\n0\n"
               "simsize.h5\nsimgrid.h5\ninitial.h5\n0\n0\n0\n{glow}\n")
        groups["ini_dtout"] = ini.format(dtout=cadence, glow="0")
        groups["ini_dtglowout"] = ini.format(dtout="60", glow=f"1\n1\n{cadence}")
        for name, config in groups.items():
            path = args.work / (name + (".ini" if name.startswith("ini_") else ".nml"))
            path.write_text(config)
            result = subprocess.run([str(args.exe), str(path)], capture_output=True, text=True, timeout=10)
            assert (result.returncode == 0) == accepted, (name, cadence, result.stdout, result.stderr)
            if not accepted:
                assert "filename cadence must be at least 0.01" in result.stderr, result.stderr
            cases += 1
    print(f"{cases} namelist/INI filename cadence contracts passed")


if __name__ == "__main__":
    main()
