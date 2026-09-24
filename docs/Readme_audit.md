# GEMINI pre-integration development profile

This branch merges upstream `9320b8912343ebd9bfb10a68c1e987f2278ff147`, fetched on 16 September 2026. It is a review candidate. The required native reference gate passes locally after correction of the configured density floor. See [the remediation handoff](PREINTEGRATION.md) for the remaining work.

## Behavior and compatibility

- Recompile C/C++ callers: boolean flags at the C boundary use `int`, including the `get_config_vars_C` output. Null configuration/work/grid deallocation is idempotent. `c_params.fortran_nml=0` is rejected because that structure cannot describe a complete GEMINI configuration.
- The tested working precision is real64; other precision settings are rejected. Normal analysis fields remain float32. The new `/restart_core` record preserves float64 plasma fields and full 3D potential; mode-specific auxiliary state is not certified complete.
- The configured `mindens_userval` now controls the physical ion-density floor. A hidden universal 1,000 m^-3 ESF floor introduced upstream on 13 August 2026 is removed. Cases requiring that stabilization must explicitly set `mindens_userval=1e3`; this is a numerical behavior change that needs case-owner review. Null-cell and division floors remain distinct.
- All seven species participate in input/output finite checks and temperature checks. Spatial ghost cells are excluded. Finite nonnegative low-density partitions are accepted; a per-worker maximum-density threshold incorrectly rejected underground restart partitions.
- Required calendar dates, fields, finite values, modes and positive driver cadences are checked. Optional namelist defaults are reset on each read. Deprecated INI input is outside this qualification.
- Calendar search now handles month/year/leap boundaries without an unbounded stepping loop. Driver priming uses elapsed seconds across midnight, and updates catch up over multiple data cadences before interpolation.
- Time steps stop at output/final times and enabled file-driver knots. The existing initialization microstep is capped at one microsecond and never raises a smaller stability-limited step. A global post-electrodynamics CFL check aborts before fluid advancement if the updated state is unstable. Automatic retry/rollback is not implemented.
- Longitude normalization accepts equivalent negative/positive angles; inverse trigonometric arguments are clamped against roundoff; 3D-to-2D coordinate calls reset their dimension state.
- The year-dependent magnetic pole table still covers 1920-2025 only. The guard rejects newer epochs. Legacy fixed-pole mode is still available, but this branch does not implement or qualify IGRF-14 geometry.
- `potsolve=2` is rejected as unimplemented. Modes 0, 1 and 3 have different scientific contracts. Restart discovery rejects legacy mode-3 files with slab-only potential; the new core record carries full 3D potential.
- Restart discovery requires density, parallel velocity, temperature, potential and date/time datasets, and checks `/time/ymd` and `/time/UThour` against the filename. This detects common incomplete/corrupt checkpoints, but does not establish complete-state or bitwise restart equivalence.
- One-point `interp1` axes are invariant dimensions. Driver objects now reject uncovered target sites by default. Explicit opt-in allows warned zero fill with in-memory masks; application-specific handling and pole singularities still need qualification.
- The GLOW data path separator is corrected. MSIS2 still requires its parameter file in the working directory; use the launcher or run from the build/install directory with the supplied data.
- CMake now verifies hashes for pinned dependency archives and offers `scripts/offline_libraries.cmake` to prepare source dependencies, including nested MUMPS. System libraries/toolchains and test-data caches remain separate prerequisites.

## Tested local profile

Linux x86-64; GNU Fortran 13.2; GCC/G++ 13.3; OpenMPI 4.1.6; HDF5 1.10.10 compatibility and 2.2.0 current-library profiles; CMake 3.31.10; Ninja 1.13; real64; GLOW and MSIS2 enabled; HWM14 disabled in the pinned nine-case reference profile, with separate HWM14 unit and smoke checks. Debug adds Fortran bounds/runtime checks. This records the tested environment, not a claim that these are the latest distribution package versions.

The sandbox required OpenMPI `pml=ob1`, `btl=self,vader`, and `btl_vader_single_copy_mechanism=none`. These are local transport settings, not required production defaults. Multi-node communication is untested.

## Build and verification

Install compiler, MPI, HDF5, BLAS, LAPACK and ScaLAPACK prerequisites. Python native probes use the standard library; optional output analysis needs numpy/h5py/matplotlib.

```sh
cmake -Doffline_dir=/absolute/dependency-cache -P scripts/offline_libraries.cmake
cmake -S . -B build -G Ninja -C /absolute/dependency-cache/sources.cmake \
  -DCMAKE_BUILD_TYPE=Release -Dgemini3d_realbits=64 -DMPIEXEC_MAX_NUMPROCS=2
cmake --build build --parallel 2
ctest --test-dir build --output-on-failure --output-junit native-results.xml
```

The source preparation step requires network unless all seven verified archives are supplied in `-Darchive_dir=/absolute/archives`. Configuration then uses prepared source overrides with `FETCHCONTENT_FULLY_DISCONNECTED=ON`. Run CTest once with network to populate reference data, or provide its verified cache and set `GEMINI3D_OFFLINE=1`. A complete release run must keep `gemini3d_test_simulations=ON`; unit-only runs do not satisfy that gate. Do not relax tolerances or replace references with current output merely to pass.

## Saved state

Species: O+, NO+, N2+, O2+, N+, H+, e-. HDF5 field order: `(species,x3,x2,x1)`. `Phiall` is a boundary slab; perpendicular velocity outputs are weighted averages. Coordinates, basis, units, masks and epoch require explicit interpretation before downstream use. This branch contains no CVTWIN or coupling changes.


See [the qualification continuation](qualification/README.md) for all new corrections, reproduction commands, the IGRF-14 oracle, optional-mode limits, instrumentation and remaining gates.
