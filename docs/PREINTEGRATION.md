# GEMINI branch handoff - 17 September 2026

## 22 September 2026 audit remediation

This addendum supersedes the historical delivery/status statements below for
the new source-audit work. The starting revision was
`bc9d1a17ff89074586d7f6064b21ceab638dbbb3`. A file inventory is not an exhaustive
code review or scientific certification.

The implementation plan now includes a local installation stage: provision
native prerequisites explicitly, create an isolated Python/build/install tree,
test before declaring installation successful, and check the recorded
environment before model execution. See
[local installation](Readme_prereqs.md#verified-local-environment).
Debian/Ubuntu and macOS use the shell entry point; Windows uses Ubuntu WSL through
the PowerShell entry point. Native Windows and additional Linux distributions
are not automatically provisioned by these scripts.

Implemented components:

- A1–A3: valid Cartesian coordinate lifetimes, completed MPI corner receives,
  corrected magnetic boundary integration and boundary-current cleanup.
- A4–A7: guarded neutral extrapolation and current-background winds,
  nondestructive launching, explicit diffusion solver errors, and a 10 ms
  minimum cadence for legacy timestamped files.
- A8–A10: owned input-data/grid teardown, restored neutral C API entry points
  and explicit interoperable integer kinds, valid odd/small-grid MPI selection.
- A11–A13: required qualification dependencies, hosted qualification definitions,
  and strict conservative-remapping validity masks.
- Exact-restart schema 2 binds sidecars to rank, layout, root/generation, inputs
  and executable, and rejects nonfinite state. This intentionally rejects old
  schema-1 exact-state checkpoints; it is not cryptographic payload authentication
  or a new power-loss durability guarantee.
- Targeted CTest regressions, source/dependency inventories, installation
  records, and existing scientific/optional-mode guards are retained.

Validation and remaining acceptance:

- The local isolated Python packages installed successfully; dependency checks
  and nine installer orchestration regressions passed. The standard full native
  installer stopped when `mumps-solver.org` could not resolve for the pinned
  MUMPS 5.9.1 archive. No successful installation record was emitted.
- Eleven isolated GNU numerical/runtime checks passed, including one-, two- and
  four-rank MPI and diffusion convergence; restoring the original defects caused
  their targeted regressions to fail. These are not full native reference runs.
- Thirty-four qualification/dependency/mask tests passed. Workflow syntax and
  shell checks are definition validation, not hosted execution evidence.
- Additional **diagnostic, non-qualifying** builds exercised 11 contract tests
  (40,248 native/CMake partition cases, 45 cadence cases, 16 launcher cases,
  neutral ABI behavior and public-symbol linking), and 72 synthetic restart
  checks on two/four ranks, including 62 expected corruption rejections.
  The contract build used system MUMPS 5.6.2; restart/lifecycle builds used a
  version-5.9.1 mirror not verified against the pinned archive. These results
  must be repeated with approved source provenance.
- Lifecycle diagnostics passed 23 Debug checks and 18 ASan/UBSan/leak checks,
  including 800 input-object states and 40 meshes, with optional physical models
  disabled. Review identified an additional temperature-placeholder regression;
  the correction admits zero-filled cells only until extrapolation/filling and
  validates the resulting temperatures. Four zero-placeholder follow-up cases
  (ascending, descending, closed and below-ground) passed with floating-point traps.
- LeakSanitizer accepted the clean control and detected the intentional
  1237-byte leak. Actual full-application leak qualification remains separate.
  Kernel-memory measurement was blocked by unavailable cgroup delegation.
- Complete pinned builds, all nine original reference cases, final-candidate
  restart/numerical-accounting campaigns, hosted Debug/Release/current-HDF5/
  sanitizer jobs, macOS/WSL/Intel execution, and target-host performance remain
  required before release acceptance.
- The first pushed candidate's hosted runs returned `action_required`, with no
  jobs or failure logs; maintainers must approve workflow execution. This is not
  a passing or failing test result.
  Separate local-installation smoke jobs exercise the installer and subsequent
  runtime preflight on Ubuntu, macOS and Windows WSL; their execution is also
  required before claiming those installation paths verified.
- Automated combined review/security validation timed out. Do not interpret
  that as a clean security result; complete it on a capable CI host.
- Independent physical review, current-epoch geometry, full coupled conservation
  and refinement, optional-mode reference evidence, and target-specific
  dependency/license/advisory approval remain open in the
  [qualification register](qualification/REMAINING_WORK.md).

Reference archives, hashes and comparison tolerances were not relaxed. No
software-only result closes the independent scientific gates.

### Follow-up installation and acceptance verification

The follow-up started from `0c759dacdba497a74cad50fedc811123c8a42e4a`;
the results below supersede assumptions that the installation or hosted
qualification had already completed:

- Exercised the actual Linux `install-local.sh --system-deps` entry point on
  Ubuntu 24.04. Native prerequisites and all five pinned local Python packages
  installed, and `pip check` passed. Configuration detected GNU 13.3, OpenMPI and
  HDF5 1.10.10, but the approved MUMPS 5.9.1 URL still failed DNS resolution.
  The installer exited unsuccessfully, emitted no success record, and both
  `check` and `run` refused this incomplete environment. No replacement archive
  or weakened provenance check was used.
- The PowerShell installer now forwards local-root, jobs, Debug/Release,
  reference-test and offline-source-cache options. Its orchestration regression
  checks argument boundaries with spaces, invalid options, and a failed WSL
  installation. All ten local-environment tests passed on Linux with PowerShell.
  This is wrapper validation, not a successful installation inside Windows WSL.
- The Windows smoke workflow now invokes the actual PowerShell entry point.
  macOS and Windows end-to-end execution remain pending on their target hosts.
- Fixed the remaining unquoted manual-grid/start/end-time values in the native
  launcher and expanded its preservation/argument regression from 16 to 36
  cases. Running that regression remains blocked on the native build. The
  existing native and CMake partition probes passed 34,848 and 5,400 cases,
  respectively; review of the other A1–A10/restart corrections found no further
  high-confidence defects in this pass.
- Preintegration now explicitly uses Bash with pipeline failure propagation,
  so a failed build or test cannot be hidden by successful log capture.
  Hosted evidence also requires Debug restart/numerical-budget execution, not
  just Release. Fifty-four tests passed across fourteen standalone Python
  qualification modules, including all eight hosted-evidence tests, with no
  skips. These are component/contract checks, not native simulation evidence.
- Recorded an inventory of 355 source files and a partial native dependency/
  toolchain inventory. Hash inventories identify the reviewed source but are
  not exhaustive line coverage, complete build output, or scientific approval.
- Hosted installation and preintegration runs for the starting candidate reported
  `action_required`; the installation run had zero jobs and no failure logs.
  Maintainer approval is needed before these can supply execution evidence.
- LeakSanitizer again accepted the clean control and detected the intentional
  1237-byte leak. The delegated-cgroup probe was blocked by permissions.
  Neither control substitutes for full-application leak or performance results.
- The advisory lookup found no reported vulnerabilities for the five pinned
  local Python requirements. Native/transitive and target-distribution advisory
  and license approval remain open.

Release acceptance is therefore still blocked on an approved, reachable MUMPS
archive/cache, complete native reference/restart/accounting and sanitizer runs,
approved hosted platform execution, and the independent scientific gates above.

## Historical handoff

Branch `fix/preintegration-validation-2026-09-16` now merges upstream `9320b8912343ebd9bfb10a68c1e987f2278ff147` and retains the earlier remediation commit `ec5a1382a69e616e698aa36538958d38eb423e2d`. The delivered verification manifest identifies the final candidate commit, tree, test results, artifact hashes and any environment-limited checks.

The new [qualification profile](qualification/README.md) is the authoritative continuation guide. It covers strict driver coverage, timestamp preflight, precision-preserving core restarts, the frozen-background restart correction, curved-gradient and diffusion checks, capacitance-mode validation, an IGRF-14 oracle, GNU sanitizers, a C17 consumer, HDF5 2.2.0 and actual environment inventories. The nine pinned reference cases, archive hashes and comparison tolerances are retained. All seven saved species remain mandatory in full-output comparisons.

The [remaining-work register](qualification/REMAINING_WORK.md) retains all 18 gates with completed work and specific closure requirements. A [provisional model card](qualification/MODEL_CARD.md) makes the missing application decisions concrete. [Primary sources](qualification/RESEARCH.md) explain the current dependency and scientific evidence choices.

This is a branch-review handoff. Full scientific qualification and parallel-model integration are not approved by local software tests. Hosted CI, independent physical benchmarks, current-epoch production geometry, full-system source/flux budgets and complete-state restart acceptance, platform/target qualification and distribution approval remain explicit. HWM14 and optional-mode smoke results are not promoted to independently reviewed reference truth.

To review the supplied bundle, clone its named branch and inspect the two-parent merge and changes against `origin/main` at the pinned upstream commit. The package includes a patch against that upstream base, source snapshot, dependency/reference archives, notices, evidence and commands. No remote push has been performed.
