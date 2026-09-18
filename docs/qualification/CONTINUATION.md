# Remaining-gate implementation

The latest request authorizes implementation and evidence collection across the fourteen open applicable gates. The original eighteen-item register and its prerequisite relationships remain in force. Independent review, held-out observations and hosted execution cannot be substituted with a software test. `STAGE_B_PROTOCOL.json` records this authorization and the new component contracts; `ACCEPTANCE_BUDGETS.json` is unchanged.

## Numerical correction and accounting

`ETD_uncoupled` previously used `P/L*(1-exp(-L*dt))` above a small threshold and dropped the loss below it, including all negative coefficients. Near-zero cancellation also corrupted the production term. The replacement evaluates the analytic phi function by a convergent series near zero. An independent 80-digit Decimal oracle exercises 107 positive, negative, zero, threshold, equilibrium and scaled cases. The normalized error budget is 5e-14. This is a numerical solver correction, not a recalibration of reaction coefficients.

The optional native source ledger records the aggregate production and analytically integrated loss of each frozen-coefficient ODE: six species-number, six parallel-momentum and seven internal-energy equations. Ionization production and electron external heating are separately identified within the aggregate production. Loss is integrated from the coefficients and prior state, never defined as the observed update residual. Physical cells use the existing metric-volume weights. Standard-energy temperature floors record their internal-energy additions. The checker requires every quantity/species on every step and matching MPI/continuity clocks. Transport and continuity checkers now reject missing per-step coverage that a complete earlier step could conceal.

This advances R04 but does not finish it. Individual reactions and collision partners, pressure/electric/gravity work, conductive boundary fluxes, compression, kinetic/internal-energy conversion and complete coupled space/time refinement still require a consistent ledger. The current floor records are the explicit standard-energy diffusion/source clamps, not every possible cleanup in optional modes. No whole-plasma energy-conservation claim is made.

## Geometry decision

The full-grid screen evaluates the pinned Cartesian field against IGRF-14 at 2026, 2027, 2028, 2029 and 2030, using all physical cells. The existing implementation imposes a vertical -50 uT field. The proposed limits are 5% magnitude error and 1 degree direction error, fixed before this screen. A failure retains the production epoch guard.

The implementation path is an explicitly quantified orthogonal geometry approximation or a generalization of the entire geometry/transport system. Editing only the pole epoch or field magnitude cannot establish field alignment. A new domain with an acceptable approximation must be demonstrated, not selected by hiding failed cells. Current idealized historical cases remain valid as numerical benchmarks without being presented as realistic magnetic geometry.

## Data, observations and learned models

`campaign.py` verifies actual file hashes, exact solver identity, all four event splits, common event/forcing/initial-state/content lineage and training-only normalization. Scientific corpus admission checks R02, R04 and R07 evidence. It deliberately remains blocked while the physics and numerical prerequisites are open. The matched-rollout evaluator requires complete shapes, finite predictions, provenance and measured runtime, and includes persistence. It does not train a surrogate or turn caller-supplied timings into a deployment certificate.

`observation_operator.py` implements normalized radar line-of-sight quadrature for explicit ECEF velocity and look vectors, rejects invalid support, and provides a linear Gaussian analysis with Joseph-form covariance and normalized innovation diagnostics. These are building blocks. A calibrated real instrument, beam/range/time response, noise correlations, representativeness, event provenance and independent validation must still be supplied. Public AMISR data access is a source route, not an automatically matched observation of the idealized reference forcing.

Event-level conformal calibration requires one maximum normalized rollout score per independent event. Insufficient finite order statistics, too-wide intervals, population mismatch and excursions beyond the declared training-feature envelope produce abstention. Calibration and test events must be disjoint. Exchangeability is an assumption requiring evidence; adjacent time windows are not independent events. Coverage output is empirical, not a population confidence guarantee. No model has been qualified or control enabled by these component tests.

## Interchange and coupling

`state_exchange.py` exports the nine supported full-output science fields from the selected Cartesian, unit-metric profile with exact stored precision, SI units, UTC, native species/axes, physical mask, native basis, coordinates and input/output hashes. It does not silently substitute the more precise restart state. The native basis is explicitly expressed in geomagnetic Cartesian coordinates; it is not mislabeled as geographic ECEF. The native filename rounds time to 10 milliseconds; the interchange checks that convention, including calendar rollover, and separately preserves the finer stored UTC hour and its exact floating-point representation. The output is an interchange record, not a runnable checkpoint. Round-trip restoration reconstructs the supported science fields, not runtime sidecars or user-defined output channels.

Direct conversion to CVTWIN's six-channel nondimensional periodic incompressible MHD state raises an error. GEMINI's imposed magnetic background is not an evolved perturbation, and density-averaged perpendicular velocities are not every species' velocity. Required dimensional scales, closure, boundary conditions and exchanged sources are unresolved. The two codebases may share validation, provenance and data interfaces before a physical coupling exists.

`conservative_exchange.py` implements overlap-volume remapping of Cartesian cell averages on matching domains and an exact passive two-reservoir exchange step. Manufactured tests require conserved integrals, constants, contraction and the exact stiff exchange solution. Curved metrics, masked cut cells and mismatched domain bounds are rejected. This is component verification, not a stable coupled GEMINI/CVTWIN simulation.

## Runtime, distribution and external closure

The current-HDF5 workflow now requests three complete verified launches in fresh cgroup v2 children with kernel memory.peak accounting, the unchanged 600 s / 2048 MiB budgets, and no sampling fallback. A host without writable delegation is reported as blocked. Charged cgroup memory includes cache/kernel allocations and is distinguished from sampled aggregate process RSS. Workflow definitions are not hosted results.

The exact-restart profile now explicitly rejects the alternate-energy and density/potential-only executables before simulation. This closes a profile enforcement gap; it does not qualify those optional physics modes. R17 requires an explicit independent benchmark for each newly admitted mode, compiler or multi-node layout. Existing GLOW/MSIS2 goldens are regression evidence, while the separate HWM14 build remains a component/smoke profile without matched full-simulation references.

The numeric HDF5 profile rejects virtual datasets, external raw storage, indirect paths, non-allowlisted filters and oversized declared allocations before loading payloads. Python plugin search paths are cleared for these qualification tools, and verified native launches disable plugin loading. The profile is intended for a dedicated qualification process and changes that process's HDF5 plugin policy; it is not a thread-safe policy switch for an existing shared service. Native parser isolation and a target-specific advisory/license review remain separate.

R01 needs an actual maintainer review of the ESF floor and the new ETD correction. R07 needs an independent physics assessment and held-out matched measurements/model data. R08 needs the user-designated writable GitHub destination and successful required jobs, including leak checking on a capable host. R09 needs the distribution owner's actual inventory/notice/advisory review. Templates and hashes cannot create those approvals.
