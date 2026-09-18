# Pinned dependencies and distribution scope

`cmake/libraries.json` records source URLs, commit archives and SHA-256 values. `dependency-inventory.json` records observed license notices. `scripts/offline_libraries.cmake` verifies seven archives, including MUMPS 5.9.1 beneath the pinned MUMPS superbuild, and emits a source-cache initializer. A fresh Debug build was completed using these prepared sources with FetchContent downloads disabled.

This is source dependency closure for the tested Linux profile. It does not package the compiler, MPI, system HDF5/BLAS/LAPACK/ScaLAPACK, Python environment, operating-system patches, or reference archives. Those versions must also be recorded for a release. It is not a complete software-bill-of-materials or a security-advisory clearance.

| Component | Observed license/notice | Local build |
| --- | --- | --- |
| GEMINI | Apache-2.0 | Enabled |
| ffilesystem | MIT | Enabled |
| h5fortran | BSD-3-Clause | Enabled |
| MUMPS superbuild | MIT | Enabled |
| MUMPS upstream | CeCILL-C; bundled exceptions and PORD notices | Enabled |
| MSIS wrapper | Apache-2.0; retain underlying model notices | Enabled |
| HWM14 | Apache-2.0 | Separate unit/native-smoke profile; no pinned reference cases |
| GLOW | GLOW Open Source Academic Research License Agreement | Enabled |

GLOW's bundled agreement limits the granted purposes and includes restrictions on sale/fee transfer and derivative-work conditions. Its terms differ from GEMINI's Apache-2.0 license. Preserve the agreement and review the intended distribution/use before a commercial product release. This is an observed dependency restriction, not an assertion that the whole combined package has one license. The exact agreement and MUMPS/PORD notices are retained in `dependency-notices/` for review.

GitHub Actions are pinned to resolved commit hashes. The former `Vampire/setup-wsl@v7` reference did not resolve; it now uses the commit behind `v7.0.0`. Runner families are explicit. The new `preintegration.yml` runs Debug and Release with native reference tests enabled, and preserves diagnostics on failure. Hosted CI has not been executed from this local branch. Distribution-provided package updates still need normal release maintenance.


## Qualification continuation, 17 September 2026

The current-library profile builds HDF5 2.2.0 from the publisher's verified source archive with Fortran, HL and deflate enabled. `scripts/build_qualification_hdf5.cmake` reproduces that profile and rejects corrupt archives or reused build trees. The prior HDF5 1.10.10 profile is retained only as tested compatibility evidence. HDF5's notice is included in `dependency-notices/HDF5-2.2.0-LICENSE.txt`.

`scripts/qualification/inventory.py` records actual CMake library paths/hashes, available package archive versions, Python distributions and model/input data. The issued package contains local inventories for both profiles. These are local-environment inventories, not a substitute for the target deployment inventory or dependency-license approval. The targeted advisory review is in [qualification research](qualification/RESEARCH.md).

Hosted validation now includes Debug/Release native tests, current HDF5 and GNU sanitizer jobs with a separate standalone leak check. Local process-inspection restrictions prevent a successful LeakSanitizer run here. None of the hosted job definitions is represented as a remotely executed result.
