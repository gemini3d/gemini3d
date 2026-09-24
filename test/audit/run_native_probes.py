"""Targeted negative/edge probes of compiled GEMINI, not physical validation.

Usage: python audit/run_native_probes.py --build /path/to/cmake-build --output results.json
Run in the same compiler/MPI/HDF5 runtime environment used to build GEMINI.
"""
from pathlib import Path
import argparse
import json
import math
import subprocess
import tempfile

BASE = """&base
ymd=2013,2,20
UTsec0=18000
tdur=300
dtout=60
activ=100,100,4
tcfl=0.5
Teinf=1000
/
&flags
potsolve=1
flagoutput=1
/
&files
indat_size='inputs/simsize.h5'
indat_grid='inputs/simgrid.h5'
indat_file='inputs/initial_conditions.h5'
/
"""


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--build', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    bins = args.build.resolve() / 'test/audit'
    results = []
    with tempfile.TemporaryDirectory(prefix='gemini-audit-') as tmp:
        tmp = Path(tmp)

        def run(name, exe, cmdargs, succeeds, predicate=None):
            r = subprocess.run([str(bins / exe), *map(str, cmdargs)], cwd=tmp,
                               text=True, capture_output=True, timeout=25)
            ok = (r.returncode == 0) == succeeds
            if predicate is not None:
                ok = ok and predicate(r.stdout + r.stderr)
            results.append(dict(name=name, passed=ok, returncode=r.returncode,
                                expected_success=succeeds, stdout=r.stdout, stderr=r.stderr))

        run('null_partition:low_density_restart', 'audit_sanity', ['ns', 1, 'null_partition', tmp], True)
        for field in ('ns', 'v1', 'v2', 'v3', 'ts'):
            run(f'{field}:ghost_nan_ignored', 'audit_sanity', [field, 1, 'valid', tmp], True,
                lambda s: 'ACCEPTED' in s)
            for species in range(1, 8):
                run(f'{field}:nan:species{species}', 'audit_sanity', [field, species, 'nan', tmp], False,
                    lambda s: 'finite' in s.lower())
        for species in range(1, 8):
            run(f'ns:negative:species{species}', 'audit_sanity', ['ns', species, 'negative', tmp], False,
                lambda s: 'negative' in s.lower() or 'density' in s.lower())

        for species in range(1, 8):
            run(f'ts:negative:species{species}', 'audit_sanity', ['ts', species, 'negative', tmp], False,
                lambda s: 'negative temperature' in s.lower())

        def floor_ok(s):
            rows = [line.split() for line in s.splitlines() if 'AUDIT_CONFIG' in line]
            return bool(rows) and all(math.isclose(float(row[-1]), 1e-100, rel_tol=1e-12, abs_tol=0)
                                      for row in rows)

        good = {
            'minimal_defaults': BASE,
            'inline_indented_uppercase': BASE + '\n  &NEUTRAL_BG msis_version=21 /\n',
            'inline_comment': BASE + '\n&evibcool ! group comment\nflagevibcool=1\n/\n',
            'empty_optional': BASE + '\n&neutral_BG\n/\n',
        }
        for name, value in good.items():
            f = tmp / (name + '.nml'); f.write_text(value)
            pred = floor_ok
            if name == 'inline_indented_uppercase':
                pred = lambda s: floor_ok(s) and any('AUDIT_CONFIG' in x and x.split()[1] == '21' for x in s.splitlines())
            run('config:' + name, 'audit_config', [f], True, pred)

        first = tmp / 'saved_first.nml'
        second = tmp / 'saved_second.nml'
        first.write_text(BASE.replace('flagoutput=1', 'flagoutput=1\nflagperiodic=1') + '\n&neutral_BG msis_version=21 /\n')
        second.write_text(BASE + '\n&neutral_BG /\n')
        def saved_reset(s):
            rows = [line.split() for line in s.splitlines() if 'AUDIT_CONFIG' in line]
            return (len(rows) == 2 and rows[0][1] == '21' and rows[1][1] == '0'
                    and rows[0][-2] == '1' and rows[1][-2] == '0' and floor_ok(s))
        run('config:repeated_read_resets_defaults', 'audit_config', [first, second], True, saved_reset)

        bad = {
            'missing_activity': BASE.replace('activ=100,100,4\n', ''),
            'missing_file': BASE.replace("indat_grid='inputs/simgrid.h5'\n", ''),
            'zero_cfl': BASE.replace('tcfl=0.5', 'tcfl=0'),
            'large_cfl': BASE.replace('tcfl=0.5', 'tcfl=1.1'),
            'negative_duration': BASE.replace('tdur=300', 'tdur=-1'),
            'nan_duration': BASE.replace('tdur=300', 'tdur=NaN'),
            'inductive_not_implemented': BASE.replace('potsolve=1', 'potsolve=2'),
            'invalid_output': BASE.replace('flagoutput=1', 'flagoutput=9'),
            'zero_precip_cadence': BASE + '\n&precip dtprec=0, prec_dir="p" /\n',
            'zero_efield_cadence': BASE + '\n&efield dtE0=0, E0_dir="e" /\n',
            'zero_density_floor': BASE + '\n&mindens_user mindens_userval=0 /\n',
            'nan_density_floor': BASE + '\n&mindens_user mindens_userval=NaN /\n',
        }
        bad.update({
            'invalid_month': BASE.replace('ymd=2013,2,20','ymd=2013,13,20'),
            'invalid_leap_day': BASE.replace('ymd=2013,2,20','ymd=2013,2,29'),
            'unsupported_msis': BASE+'\n&neutral_BG msis_version=99 /\n',
            'invalid_periodic': BASE.replace('flagoutput=1','flagoutput=1\nflagperiodic=2'),
            'unsupported_fang': BASE+'\n&fang flag_fang=100 /\n',
            'invalid_diffusion': BASE+'\n&diffusion diffsolvetype=3 /\n',
            'invalid_capacitance': BASE+'\n&capacitance flagcap=3 /\n',
            'nan_capacitance': BASE+'\n&capacitance magcap=NaN /\n',
            'invalid_FBI': BASE+'\n&FBI flagFBI=3 /\n',
            'invalid_cooling': BASE+'\n&evibcool flagevibcool=2 /\n',
            'zero_milestone': BASE+'\n&milestone mcadence=0 /\n',
            'negative_background_flux': BASE+'\n&precip_BG PhiWBG=-1 /\n',
            'zero_background_energy': BASE+'\n&precip_BG W0BG=0 /\n',
            'invalid_spectrum': BASE+'\n&fang_pars diff_num_flux=99 /\n',
            'invalid_kappa': BASE+'\n&fang_pars kappa=2 /\n',
            'nan_kappa': BASE+'\n&fang_pars kappa=NaN /\n',
            'missing_precip_dir': BASE+'\n&precip dtprec=60 /\n',
            'missing_efield_dir': BASE+'\n&efield dtE0=60 /\n',
            'missing_solar_dir': BASE+'\n&solflux dtsolflux=60 /\n',
            'missing_neutral_dir': BASE+'\n&neutralBG_file dtneuBGfile=60 /\n',
            'oversized_cadence': BASE+'\n&precip dtprec=86401, prec_dir="p" /\n',
            'nonfinite_neutral_cadence': BASE+'\n&neutral_BG flagneuBG=.true., dtneuBG=NaN /\n',
            'invalid_neutral_interpolation': BASE+'\n&neutral_perturb dtneu=60, source_dir="n", interptype=2 /\n',
            'invalid_source_latitude': BASE+'\n&neutral_perturb dtneu=60, source_dir="n", sourcemlat=91 /\n',
        })
        for name, value in bad.items():
            f = tmp / (name + '.nml'); f.write_text(value)
            run('config:' + name, 'audit_config', [f], False, lambda s: 'config:' in s or 'impossible' in s)

    payload = dict(schema='gemini.audit.native_probes.1', passed=all(r['passed'] for r in results),
                   count=len(results), results=results)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + '\n')
    print(json.dumps({k:v for k,v in payload.items() if k != 'results'}))
    for r in results:
        if not r['passed']:
            print('FAILED', r['name'], r['returncode'], (r['stdout'] + r['stderr'])[-600:])
    return 0 if payload['passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
