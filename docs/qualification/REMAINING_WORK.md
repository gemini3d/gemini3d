# Remaining qualification work

This file records program obligations and component progress. The issued package GATES.json and verification report determine actual candidate closure. The latest user request authorizes work across the fourteen open applicable gates; independent review and physical evidence are not waived. See CONTINUATION.md and STAGE_B_PROTOCOL.json.

## R01 - open_with_component_progress

Completed components: Nine original native reference cases preserved; stable ETD analytic probes and review notes added.

Required next: Maintainer reviews both the case-specific ESF floor and stable ETD correction; approval must identify the actual branch.

Owner: Numerical maintainer. Acceptance: Keep all nine native reference cases green with unchanged tolerances and archives.

## R02 - closed

Completed components: Application scope and budgets selected under the latest user authorization; see MODEL_CARD.md and ACCEPTANCE_BUDGETS.json.

Required next: No additional scope/budget approval requested. Independent physics review remains R07.

Owner: Project owner and plasma physicist. Acceptance: A specific ionospheric use case and quantitative acceptance budgets are signed off.

## R03 - open_with_component_progress

Completed components: Full-grid IGRF-14 screen added for 2026-2030, with proposed 5% magnitude / 1 degree direction limits; production epoch guard retained.

Required next: Replace the inaccurate vertical-field approximation with a reviewed, coordinated field/basis/metric/driver formulation and validate the complete domain.

Owner: Geomagnetic and grid specialist. Acceptance: Independent official-calculator comparison over the declared 2026-2030 domain.

## R04 - open_with_component_progress

Completed components: Stable ETD correction; independent aggregate number/momentum/internal-energy production and loss ledgers; explicit standard-energy temperature-floor additions; strict per-step coverage checks.

Required next: Finish individual chemistry/collision channels, work and conductive boundary terms, all cleanup/kinetic/internal-energy accounting, and full coupled space/time refinement.

Owner: Numerical lead. Acceptance: Predeclared mass/species/charge/momentum/energy and convergence budgets pass.

## R05 - maintain_selected_profile_qualification

Completed components: Strict efield/precip required-field schemas, axis ordering, categorical flags, positive precipitation, full input hashes, native coverage and declared MPI layouts are covered by the research profile.

Required next: Maintain the unchanged scoped contract; exact-candidate execution evidence is supplied by the issued package register.

Owner: Solver engineer. Acceptance: Missing coverage cannot become an unlabelled physical zero; tests cover target MPI layouts.

## R06 - maintain_selected_profile_qualification

Completed components: Atomic root publication plus rank-local full-step state; saved fluid/electric/interface arrays, dt/deadlines; reconstructed frozen backgrounds and immutable drivers. Native trajectory/fault matrix targets frozen bounds.

Required next: Maintain the unchanged scoped contract; exact-candidate execution evidence is supplied by the issued package register.

Owner: I/O engineer and numerical lead. Acceptance: Restart trajectories meet declared budgets; field-resolved modes have full-potential checkpoints.

## R07 - open_with_component_progress

Completed components: Observation and lineage tooling, unchanged physical targets, public AMISR access route and current model-validation publications reviewed.

Required next: Independent specialist must review closures/rates and compare actual held-out matched events or independent models with calibration and uncertainty metadata.

Owner: Independent plasma-physics reviewer. Acceptance: Independent review and held-out benchmarks meet the approved application budgets.

## R08 - open_with_component_progress

Completed components: Current-HDF5, native, sanitizer and standalone leak workflows retained; kernel-memory workload jobs added.

Required next: Execute required workflows at the user-designated writable GitHub destination, including LeakSanitizer on a capable host.

Owner: Build engineer. Acceptance: Green required hosted jobs without skipped native scenarios; sanitizer and embedding reports.

## R09 - open_with_component_progress

Completed components: Current inventories/notices retained; September 2026 HDF5 audit reviewed; numeric HDF5 policy rejects VDS, external storage, indirect paths, unsupported filters and excessive allocation.

Required next: Distribution owner reviews the actual target inventory, relevant advisories and component/data notices; parser containment is a separate deployment requirement.

Owner: Release owner. Acceptance: Approved inventory/notices/advisory record for the actual distribution.

## R10 - open_with_component_progress

Completed components: Event/forcing/initial-state/content lineage, artifact hashes, four nonempty splits and training-only normalization validators implemented; prerequisite admission fails closed.

Required next: Close R04/R07 and generate the actual physical-regime corpus with complete verified metadata before scientific training.

Owner: Data/ML engineer. Acceptance: Disjoint event/forcing splits and complete manifests.

## R11 - open_with_component_progress

Completed components: Masked normalized ECEF radar-LOS quadrature and Joseph-form linear analysis with innovation diagnostics implemented and analytically tested.

Required next: Bind the operator to a calibrated actual instrument, its beam/range/time response and independent event data; validate innovations and physical budgets.

Owner: Instrument/assimilation lead. Acceptance: Independent observation and innovation diagnostics pass.

## R12 - open_with_component_progress

Completed components: Complete-rollout and measured-cost benchmark evaluator with persistence and provenance requirements implemented.

Required next: Train the specified surrogate on an admitted corpus and run fixed-budget matched numerical/simple baseline and held-out rollout comparisons.

Owner: ML lead. Acceptance: Matched error/latency/compute and held-out rollout criteria pass.

## R13 - open_with_component_progress

Completed components: Independent-event conformal quantiles, held-out coverage and width/population/feature-excursion abstention implemented; insufficient finite calibration is rejected.

Required next: Collect independent target-population calibration/test events, fix coverage and useful-width targets before evaluation, and demonstrate the targets under distribution-shift challenges.

Owner: UQ lead. Acceptance: Coverage/width criteria pass on the target population.

## R14 - open_with_component_progress

Completed components: Fresh delegated-cgroup kernel memory.peak runner and three-run hosted workload job implemented, with no sampled-memory fallback.

Required next: Run on a capable delegated host under the unchanged 600 s / 2048 MiB bounds; this local cgroup mount is read-only.

Owner: HPC/runtime engineer. Acceptance: Tail-latency and memory budgets met under declared workload.

## R15 - scope_excluded

Completed components: Physical control explicitly excluded from the research-only path selected under the latest user authorization.

Required next: No physical-control implementation or qualification in this scope; reopen if a controller is proposed.

Owner: Systems owner. Acceptance: Application-specific physical/closed-loop evidence and independent safing review.

## R16 - open_with_component_progress

Completed components: Overlap-volume Cartesian conservative remap and exact passive two-reservoir exchange component tests implemented.

Required next: Specify a physically justified two-way GEMINI/CVTWIN exchange law including reservoirs, closure and boundaries; then validate a stable coupled model benchmark.

Owner: Coupling specialist. Acceptance: Manufactured exchange, conservation and stable coupled benchmark.

## R17 - open_with_component_progress

Completed components: Exact-restart admission now rejects alternate-energy and density/potential-only executables before work. Existing mode regression evidence retained.

Required next: Supply independent optional-mode references including HWM14 full simulations, then qualify each newly admitted geometry/BC, non-GNU compiler and multi-node layout.

Owner: Maintainer and specialist reviewers. Acceptance: Explicit supported-feature matrix and corresponding targeted validation.

## R18 - open_with_component_progress

Completed components: Versioned Cartesian full-output exchange schema preserves nine science fields, native precision, SI units, UTC, axes/species, mask and explicit native basis. Direct unsupported MHD6 conversion is rejected.

Required next: Validate the agreed physical adapter after upstream scientific gates; no evolved magnetic perturbation or species-resolved perpendicular velocity may be fabricated.

Owner: Platform engineer. Acceptance: Schema/unit/basis/mask contracts and supported-output round trips pass.
