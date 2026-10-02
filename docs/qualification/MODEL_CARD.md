# Authorized GEMINI research model card

Application scope and budget selection were delegated by the user's 17 September 2026 message. The concrete decisions and their provenance are in `ACCEPTANCE_BUDGETS.json` and `APPROVED_PATH.md`. This authorization is not independent scientific review.

| Contract | Selected scope |
| --- | --- |
| Application | Reproducible idealized regional auroral response and numerical-method development before integration |
| Primary workloads | Pinned mini2dns_fang and mini3d_fang, 300 s, 20 February 2013; shifted-midnight calendar case |
| Domain | Archived idealized Cartesian grids (nominal setup 67.11 N, 212.95 E), nominal 80–1000 km; actual hashed grid coordinates define the domain. No matched radar site/event is asserted. |
| Species | O+, NO+, N2+, O2+, N+, H+, electrons; native order preserved |
| Quantities / units | Species density m^-3, temperature K, velocity m/s, potential V; numerical particle/mass/transport budgets |
| Drivers | Complete electric and precipitation HDF5 streams at declared 60 s cadence; original empirical neutral/solar configuration |
| Solver modes | real64, potsolve=1, flagcap=0, standard energy, FANG, frozen empirical neutrals, fixed grid, full output |
| Platform | Linux x86-64, GNU Fortran 13.2/GCC 13.3, OpenMPI 4.1.6, one node, 1/2 ranks; native HDF5 2.2.0 plus 1.10.10 compatibility |
| Restart | Opt-in complete-state profile, same inputs/executable/layout; duration can extend; atomic final-name publication |
| Coverage / failure | Strict spatial/time coverage; no imputation, missing-site mode or automatic retry |
| Regression | All nine pinned native cases and original comparison tolerances retained |
| Performance | Local 300 s 3D batch workload, 3 repetitions, 600 s wall and 2048 MiB resident-memory budget per run; no production percentile claim |
| Use | Offline research and branch review; physical control excluded; no operational forecast accuracy claimed |

The 2026–2030 production geometry upgrade remains R03. Updating a pole table alone would not validate the field direction, metric/basis and forcing transforms. Keep the production guard and independent IGRF-14 oracle. Retain orthogonal meshes for Stage A; before Stage B choose either a quantified orthogonal approximation with field-error limits over the specific modern domain or a coordinated generalized geometry implementation. Do not connect the full field oracle to only one part of the geometry.

Independent instrument/observation evaluation, optional-mode physics validation, multi-node/non-GNU support, corpus/assimilation/surrogate/UQ and CVTWIN integration remain later program gates. Initial physical-error targets in the budgets are requirements for a future matched-regime held-out study, not measured performance or confidence intervals.
