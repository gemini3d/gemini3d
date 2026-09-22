# Remaining-gate implementation update

See [CONTINUATION.md](CONTINUATION.md) for the stable source integrator, independent source/floor ledgers, new exchange/data/observation/UQ components, stronger HDF5 input policy, optional-executable guards and kernel-memory workflow. [STAGE_B_PROTOCOL.json](STAGE_B_PROTOCOL.json) records the latest authorization and component contracts. The issued package contains actual execution evidence and the authoritative fourteen-gate status; component tests do not close independent physics or hosted-execution requirements.

---

# Qualification continuation with authorized scope

The current application/model card and predeclared budgets are `MODEL_CARD.md`, `APPROVED_PATH.md` and `ACCEPTANCE_BUDGETS.json`. The latest user message authorizes these choices; another scope/budget approval is not required. See `RESTART_CONTRACT.md` and `NUMERICAL_ACCOUNTING.md` for the implemented continuation. `remaining_work.json` preserves the full program obligations, including work beyond this research branch.

New functionality: driver-specific required-field and axis validation; corrected Nlon reader alias; atomic HDF5 publication; opt-in full-step, input/executable/layout-bound restart; actual native face fluxes and MPI interface cancellation; whole-step ion continuity with applied sources and explicit density cleanup. Hosted CI includes the research matrix. The accompanying exact-commit verification report determines which profiles were actually executed.

The material below documents the preceding engineering increment and its limitations **at that time**. In particular, the prior core-only/non-atomic restart description is superseded for the new opt-in profile, and the prior application-approval prerequisite is resolved by the user's authorization.

---

# GEMINI qualification continuation - 17 September 2026

This branch contains engineering corrections and reproducible qualification tools. It is ready for branch review when the accompanying verification manifest is green. Scientific, operational and coupled-model qualification remain separate gates. The issued package contains the candidate commit, test logs, evidence hashes and all 18 remaining-work decisions. No gate is closed by changing a label or suppressing a failed test.

## Corrections in this continuation

- Merge upstream `9320b8912343ebd9bfb10a68c1e987f2278ff147`. Its removal of the universal density-floor override agrees with the earlier correction. Preserve explicit case configuration for ESF stabilization.
- Reject missing spatial driver coverage by default. `inputdata%coverage` retains per-target masks and `spatial_missing_count`; only explicit `&input_coverage allow_missing_spatial=.true. /` permits warned, masked zero fill. Masks are in memory, not new per-output datasets. Strict coverage is the qualified default; permitting missing sites requires an application-specific downstream mask contract.
- Reject nonfinite interpolation data and nonuniform/missing driver frames. A separate HDF5 preflight validates the entire forcing stream, timestamps, fields, shapes and links before launch.
- Preserve density, temperature, parallel velocity and full 3D potential in `/restart_core/` at float64. Validate schema, completion marker and actual dataset precision; prefer these fields to legacy analysis arrays. Existing float32 analysis fields keep their layout. A completion marker detects common partial records; it does not make writes atomic or certify every auxiliary variable.
- Keep a frozen empirical neutral background at the configured start epoch during restart. Time-updated backgrounds still use the current epoch. This corrects a regression introduced by the newly merged upstream initialization change.
- Correct the x2 boundary-gradient metric index when h2 varies across x3. A manufactured field checks both the boundary and the interior.
- Extend TRBDF2 with optional midpoint boundary values for time-varying boundary conditions, retaining existing callers' stage-constant behavior. Manufactured source and spatial/time refinement checks exercise the operator.
- Reject capacitance combinations that do not enter the implemented dynamic solver: only a 3D Cartesian domain, `potsolve=1`, and potential-gradient boundary `flagdirich=0` are accepted. Prescribed-potential, current-boundary, field-resolved and 2D capacitance combinations are unsupported. This prevents silently incomplete equations and the observed nonfinite-current failure.

## Qualification tools

| Tool | What it establishes | Limit |
| --- | --- | --- |
| `test/audit/diffusion_convergence.f90` | Spatial and temporal convergence; manufactured source and boundary balance | One diffusion operator, not all coupled equations |
| `scripts/qualification/budgets.py` | Volume-weighted species/mass/charge/momentum/energy and clipping arithmetic | Full solver source and face-flux instrumentation still required |
| `scripts/qualification/igrf14.py` | Independent degree-13 geocentric IGRF-14 field oracle | Does not replace the legacy grid, metric or basis transformations |
| `test/qualification/restart_core_checks.py` | Native precision/preference checks and malformed-core rejection | Complete auxiliary state and trajectory budgets remain open |
| `test/qualification/mode_experiments.py` | Native optional-mode execution and C++/Fortran consistency | Synthetic 120-second runs are not independent physical validation |
| `scripts/qualification/inventory.py` | Actual build/Python/package/model-data inventory with hashes | Run again on the deployment target; license/advisory approval is separate |
| `scripts/qualification/evaluate_gates.py` | Candidate/evidence/hash/approval-record integrity | Release owner must authenticate reviewer identity and authority |

The gate checker requires exactly R01-R18, an exact candidate commit, existing hashed evidence, and an approval record when required. Missing gates, modified evidence, stale commits and evidence paths outside the package fail closed. A signed scope exclusion is needed to remove a gate from the proposed release; later integration items are retained instead of silently marked passed.

## Reproduce

Use a fresh build directory after derived-type or toolchain changes. For simultaneous or externally synchronized workspaces, set `-Dgemini3d_test_run_root=/absolute/fresh-output-root` to isolate native case outputs; the default is `build/test_runs`. This directory is reserved for CTest cases and its dated outputs are reset on each run. Install GNU C/C++/Fortran, MPI, LAPACK/BLAS/ScaLAPACK, HDF5 Fortran/HL with deflate, and Python **3.11 or newer** (qualification hashing uses `hashlib.file_digest`) with NumPy, h5py and SciPy. Hosted jobs install the pinned versions in `requirements.txt`.

`-Dgemini3d_require_qualification=ON` requires `BUILD_TESTING=ON`, `gemini3d_BUILD_TESTING=ON` and all these Python dependencies; missing dependencies fail configuration. Import checks run again on every configure, including when the interpreter changes. The default is OFF for optional native-only builds: unavailable Python-dependent tests are not registered, which is **not complete qualification**. The `qualification-debug` and `qualification-release` configure presets require these dependencies and retain all nine pinned native cases.

```sh
cmake -Doffline_dir=/absolute/cache -P scripts/offline_libraries.cmake
cmake -S . -B build -G Ninja -C /absolute/cache/sources.cmake \
  -DCMAKE_BUILD_TYPE=Release -Dgemini3d_realbits=64 \
  -Dgemini3d_test_simulations=ON -Dgemini3d_require_qualification=ON -DMPIEXEC_MAX_NUMPROCS=2
cmake --build build --parallel 2
ctest --test-dir build --output-on-failure --output-junit native.xml
python test/qualification/restart_core_checks.py --exe build/gemini.bin \
  --case build/test_runs/mini2dns_fang --work /absolute/new-restart-work
python test/qualification/mode_experiments.py --build build \
  --inputs build/test_runs/mini2dns_fang/inputs \
  --inputs3d build/test_runs/mini3d_fang/inputs --work /absolute/new-mode-work
```

For a current HDF5 profile, run the verified-source builder before GEMINI configuration, then supply `-DHDF5_ROOT=/absolute/hdf5-install` in a new GEMINI build directory. Its working directory must be fresh. The package includes the SHA256-verified source archive for offline use.

```sh
cmake -Dqualification_work=/absolute/new-hdf5-work \
  -Dqualification_prefix=/absolute/hdf5-install \
  -Dsource_archive=/absolute/cache/hdf5-2.2.0.tar.gz \
  -P scripts/build_qualification_hdf5.cmake
```

For GNU address/undefined-behavior instrumentation use a new Debug build with `-Dgemini3d_sanitizers=ON`. Set `ASAN_OPTIONS=detect_leaks=0:halt_on_error=1` only when leak checking is deliberately excluded; this is not leak qualification. The supplied hosted workflow has a separate leak job step. External MPI/HDF5/BLAS libraries are not automatically instrumented.

`.github/workflows/preintegration.yml` runs on pull requests, pushes and manual dispatches without excluding qualification/script changes. It runs Debug/Release native comparisons and the full research matrix, verified HDF5 2.2.0, sanitizer checks, a deliberate-leak positive control, and unsuppressed application leak checks. Artifacts bind the checked-out candidate commit, repository, run and attempt and retain hashed diagnostics even on failure. This workflow definition is not evidence of a hosted pass.

Kernel-memory measurements require an already delegated cgroup-v2 memory controller. `GEMINI_QUALIFICATION_CGROUP` can name that parent through a repository variable; the default is `/sys/fs/cgroup`. No privilege/mount workaround, sampling substitution or skipped measurement can produce a pass. Missing delegation, LeakSanitizer capability failure, or an MPI/dependency leak remains a failed/blocked qualification requiring a capable host and triage.

The Cartesian remapper accepts only Boolean or integer zero/one validity masks of the source-cell shape. Floating-point, string, object and nonbinary integer masks are rejected rather than coerced to truth; any invalid/cut cell still rejects this all-covered remapping contract.

Run each forcing stream through `validate_driver.py` with explicit start/stop/cadence before launching. Naive timestamps mean UTC. The validator rejects gaps; it does not infer cadence, impute data or approve extrapolation.

## Supported evidence and open scope

The core profile is Linux x86-64, real64, GNU Fortran 13.2 / GCC 13.3, OpenMPI 4.1.6, GLOW and MSIS2, two ranks on one node. HDF5 1.10.10 is retained as a compatibility baseline and HDF5 2.2.0 is the current-library profile. The nine pinned reference cases, archive hashes and comparison tolerances are unchanged. Full-output comparisons cover all seven species.

HWM14 has unit and native C++/Fortran smoke evidence. Its entries are absent from the pinned native reference manifest, so enabling the required HWM14 reference suite intentionally fails configuration. New independently reviewed HWM14 reference data are needed; disabling that reference gate is not a fix.

The year-dependent production pole model still supports 1920-2025. The IGRF-14 oracle is not connected to production geometry. Full IGRF field alignment can violate the existing orthogonal-grid assumptions; a coordinated metric/basis/driver design and domain-specific checks are required before 2026-2030 production use.

Restart trajectory agreement is measured, not presumed bitwise. Neutrals/driver caches, transverse species drifts, mode-specific electrodynamic history and scheduling must be enumerated for every admitted mode. Frozen-background correction and float64 core state reduce two risks but do not prove full-state equivalence.

Application approval, full conservation/source/flux and multi-equation convergence evidence, independent physics/observation validation, non-GNU and multi-node coverage, target hardware performance, complete advisory/license approval and hosted CI are still explicit gates. Corpus, assimilation, surrogate/UQ, physical control and coupling remain dependent on a defined application and the later integration phase.
