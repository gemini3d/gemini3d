# GEMINI branch handoff - 17 September 2026

Branch `fix/preintegration-validation-2026-09-16` now merges upstream `9320b8912343ebd9bfb10a68c1e987f2278ff147` and retains the earlier remediation commit `ec5a1382a69e616e698aa36538958d38eb423e2d`. The delivered verification manifest identifies the final candidate commit, tree, test results, artifact hashes and any environment-limited checks.

The new [qualification profile](qualification/README.md) is the authoritative continuation guide. It covers strict driver coverage, timestamp preflight, precision-preserving core restarts, the frozen-background restart correction, curved-gradient and diffusion checks, capacitance-mode validation, an IGRF-14 oracle, GNU sanitizers, a C17 consumer, HDF5 2.2.0 and actual environment inventories. The nine pinned reference cases, archive hashes and comparison tolerances are retained. All seven saved species remain mandatory in full-output comparisons.

The [remaining-work register](qualification/REMAINING_WORK.md) retains all 18 gates with completed work and specific closure requirements. A [provisional model card](qualification/MODEL_CARD.md) makes the missing application decisions concrete. [Primary sources](qualification/RESEARCH.md) explain the current dependency and scientific evidence choices.

This is a branch-review handoff. Full scientific qualification and parallel-model integration are not approved by local software tests. Hosted CI, independent physical benchmarks, current-epoch production geometry, full-system source/flux budgets and complete-state restart acceptance, platform/target qualification and distribution approval remain explicit. HWM14 and optional-mode smoke results are not promoted to independently reviewed reference truth.

To review the supplied bundle, clone its named branch and inspect the two-parent merge and changes against `origin/main` at the pinned upstream commit. The package includes a patch against that upstream base, source snapshot, dependency/reference archives, notices, evidence and commands. No remote push has been performed.
