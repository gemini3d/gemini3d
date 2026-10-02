# Energy and evidence continuation - 17 September 2026

This continuation starts from a412224317eb57010cd4b2ed25889a1729636828. All earlier acceptance budgets and the fourteen open gates remain in force. Approval to perform the work does not substitute for independent review, actual matched observations or hosted execution.

## Native changes

Both 3D diffusion wrappers previously assigned only the interior of their result arrays. The caller copied the entire result into temperature, including undefined halos. The wrappers now preserve the incoming halo values explicitly. This changes no interior PDE formula or boundary type.

Optional parabolic diagnostics expose five integrated temperature increments: linear reaction, thermal drift, conduction, explicit source and prescribed endpoint reservoir change. For backward Euler, the first four are `dt * L(Tfinal)`. For the existing half-step TR/BDF formula, eliminating the intermediate state gives `dt/3 * (L(Tinitial) + L(Tmid) + L(Tfinal))`. The conduction term is computed from the stage face differences. Endpoint values are imposed algebraically; their measured changes are labeled as boundary reservoir changes and are never described as a conductive-flux measurement.

Compression records pressure work and artificial-viscous work separately, using the same midpoint state and divergence as the applied explicit update. The optional native energy-operator ledger integrates those terms over physical-cell volumes and records pre/post values and a separately recomputed residual. `check_energy_operators.py` requires each of seven species and both stages at every continuity time step and every rank. Missing later steps, duplicates, invalid scales and corrupt residuals fail.

The research matrix now checks this ledger for the five admitted MPI layouts and the C++ frontend. This is additional evidence for R04. Individual reaction/collision rates, coupled kinetic/internal-energy exchange, full momentum and energy closure, and coupled three-level refinement remain separate requirements. An exact split-operator identity is not independent plasma validation.

## Observation and uncertainty evidence

The radar operator rejects non-Boolean/non-binary validity masks, including NaN values that NumPy otherwise converts to true. Empty observations and empty rollout arrays are rejected. Calibration coverage requires the actual unique calibration-event inventory and a nonnegative finite radius.

`qualify_coverage` adds an exact one-sided Clopper-Pearson binomial lower confidence bound and a simultaneous width limit. Targets must be declared before the held-out evaluation. With all events covered, ten independent test events give only a 0.7411 lower bound at 95% confidence; 59 give 0.9505. This is a sample-size calculation, not a measured result for GEMINI. Event independence, target-population matching and instrument calibration cannot be verified from score arrays alone. The function uses SciPy 1.17.0 in the tested qualification environment.

## Hosted and host evidence

`host_capabilities.py` requires a leak-free C control to pass and an intentional 1,237-byte leak to be detected with the expected report. A LeakSanitizer runtime crash fails the capability probe. The optional cgroup probe records a blocked outcome when delegation is unavailable; it changes no mount or privilege settings.

The hosted workflow preserves all numerical ledgers, runs the leak capability controls, and adds an unsuppressed leak-enabled native application run. Any MPI/dependency leak remains visible and must be triaged; no broad suppression or disabled leak check can close R08.

`hosted_evidence.py` reads GitHub's API through an already authenticated `gh` client. It verifies the repository, exact candidate commit, workflow path, run attempt, mandatory jobs and mandatory steps. It does not push, dispatch a workflow, post comments or supply external approval. The report is only a hosted-execution record; required artifacts and scientific reviews remain separate.

## Remaining hard boundaries

The 2026-2030 geometry failure is not fixed by replacing a pole table or changing field magnitude alone. The general degree-13 IGRF field requires a coordinated basis/metric/derivative/driver treatment. No production guard has been removed, domain shrunk or acceptance budget relaxed.

No observed event, calibration, scientific corpus, trained surrogate or physical GEMINI/CVTWIN coupling is fabricated. Their prerequisites remain visible in the root gate register. The numerical maintainer and independent scientific/distribution reviewers must authenticate acceptance records themselves. A GitHub destination is required before hosted branch results can exist.
