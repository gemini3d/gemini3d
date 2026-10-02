# Research and dependency review

Reviewed 16-17 September 2026, concentrating on the outstanding GEMINI qualification items. This is a targeted primary-source update, not a claim to have exhaustively reviewed every publication.

| Primary source | Relevance and disposition |
| --- | --- |
| [GEMINI upstream](https://github.com/gemini3d/gemini3d) | Fetched and merged main commit 9320b8912343ebd9bfb10a68c1e987f2278ff147. Repository documentation describes orthogonal curvilinear geometry; retain that assumption when designing magnetic upgrades. |
| [NOAA/NCEI IGRF](https://www.ncei.noaa.gov/products/international-geomagnetic-reference-field) | Use official IGRF-14 coefficients and synthesis code as the independent field oracle. A degree-13 field is richer than the legacy centered-dipole pole orientation; replacing the pole table alone cannot qualify field-aligned metrics and driver transforms. |
| [Official IGRF-14 coefficients](https://www.ngdc.noaa.gov/IAGA/vmod/coeffs/igrf14coeffs.txt) and [Fortran routine](https://www.ngdc.noaa.gov/IAGA/vmod/igrf14.f) | Hashes and extraction provenance are committed under test/qualification/data. The independent implementation is tested over 840 geocentric cases. |
| [HDF5 2.2.0 release](https://support.hdfgroup.org/releases/hdf5/2.2.0/downloads/) | Publisher dates this release 29 July 2026. Source archive SHA256 is pinned; build uses Fortran, HL and deflate, with tools/examples disabled. This updates the tested dependency profile without claiming a complete security clearance. |
| [Ubuntu advisory](https://ubuntu.com/security/CVE-2026-19023) and [Red Hat CVE-2026-19023](https://access.redhat.com/security/cve/cve-2026-19023) | Advisory concerns h5dump variable-length string binary dumping in HDF5 before 2.1.1. Do not equate presence of an older solver library with demonstrated reachability of that utility path. Retain baseline-version triage; deploy only an approved target inventory. |
| [GNU instrumentation documentation](https://gcc.gnu.org/onlinedocs/gcc/Instrumentation-Options.html) | Supplies the address/undefined-behavior instrumentation mechanism. Local leak checking is limited by process-inspection restrictions and remains open. |
| [Huyghebaert et al., 2025](https://angeo.copernicus.org/articles/43/99/2025/angeo-43-99-2025.html) | GEMINI produces a 50 m Kelvin-Helmholtz test volume for a synthetic EISCAT_3D forward/noise/inversion experiment with regularization. This is a concrete pattern for a selected instrument operator. Synthetic reconstruction does not establish independent physical validation of this branch or the proposed application. |
| [Scalable Mesh Coupling for Atmospheric Wave Simulation, 2026 preprint](https://arxiv.org/html/2603.02971v1) | Relevant literature for a future mesh-coupling design. No coupling implementation, conservative exchange claim or validation result is imported into this GEMINI-only change. |

The newer radar work supports implementing and testing an explicit observation operator after an instrument is selected. It does not justify using raw simulator state as if it were a direct measurement. The mesh-coupling preprint is retained for the later integration review; it does not remove the need for units, basis, masks, exchange conservation and stability checks.

## Remaining-gate research update, 17 September 2026

The public GEMINI main reference was rechecked and remains `9320b8912343ebd9bfb10a68c1e987f2278ff147`. The following targeted additions inform the new components. No claim of an exhaustive literature search or external qualification is made.

| Primary source | Compare/contrast and implemented consequence |
| --- | --- |
| [Hu and Dong, CAM-NET v4, 21 July 2026](https://arxiv.org/abs/2506.19340v4) | A geometry-aware WACCM-X emulator with spherical Fourier structure and held-out autoregressive tests. Its approximately 100 km, three-hour global representation is different from GEMINI's regional auroral grid and seconds-scale stepping. The authors distinguish emulation from observational accuracy. Adopt complete held-out rollouts and training-only normalization; do not import its reported skill as GEMINI evidence. |
| [HDF Group, 2 September 2026 audit](https://www.hdfgroup.org/2026/09/02/hdf5-beyond-the-file-findings-from-our-safety-security-and-privacy-audit/) | Risks depend on enabled features, surrounding privileges and deployment. The new qualification policy checks external storage, VDS, filter authorization and allocation bounds before payload access. It does not certify the parser or constitute a complete advisory scan. |
| [AMISR data access](https://amisr.com/amisr/links/data-access/) and [user manual](https://amisr.github.io/amisr_user_manual/src/intro.html) | Public ISR data and instrument interpretation guidance are available. Archive availability does not establish matched driving conditions, calibration/uncertainty, or independence for the idealized GEMINI cases. No mismatched observations were used to claim R07/R11 closure. |
| [Barber and Pananjady, arXiv:2510.02471](https://arxiv.org/abs/2510.02471) | Time dependence and predictor memory affect conformal coverage. Use independent event-level scores, separate calibration/test lineage and explicit exchangeability assumptions. A finite quantile and empirical coverage are insufficient to certify the target population. |
| [Linux cgroup v2 documentation](https://www.kernel.org/doc/html/latest/admin-guide/cgroup-v2.html) | Kernel memory.peak supports a full descendant-tree high-water measurement in a fresh delegated group. Read-only host delegation is a real capability blocker; process samples and the shared root cgroup's historical peak are not substituted. |
| [Scalable Mesh Coupling, 2026](https://arxiv.org/html/2603.02971v1) | Mesh search, coordinate transforms and exchanged atmospheric/plasma terms are useful architectural precedents. The new overlap-volume/passive-exchange tests establish only their manufactured component contracts. No physical law is inferred between periodic incompressible MHD and multi-species GEMINI. |

The updated [continuation implementation](CONTINUATION.md) maps these choices to code, limitations and gate prerequisites.

## Energy/evidence continuation: primary-source recheck

The GitHub API was checked on 17 September 2026 and still returned upstream commit `9320b8912343ebd9bfb10a68c1e987f2278ff147`. No new upstream commit was merged in this continuation.

- [HDF Group 2.2.0 release announcement](https://forum.hdfgroup.org/t/release-of-hdf5-2-2-0-newsletter-210/13861) identifies fixes for CVE-2025-9274, CVE-2026-17572 and CVE-2026-17574, and retention of the earlier CVE-2026-17573 fix. The compatibility HDF5 1.10.10 test is not a security-clearance result. These are specific publisher claims, not a complete advisory inventory or a claim about backports in a distribution package.
- [NOAA/IAGA IGRF-14](https://www.ncei.noaa.gov/products/international-geomagnetic-reference-field) remains the official field reference. It is a degree-13 main-field model; it does not supply a GEMINI-compatible orthogonal grid or storm-time external magnetic perturbations.
- [AMISR remote-access manual](https://amisr.github.io/amisr_user_manual/src/madrigal_database/madrigal_remote_access.html) documents PFISR code 61, experiment/file queries and required attribution identity for downloads. The example data are not automatically matched to the 2013 idealized forcing. No registration identity, calibration or event match was invented.
- [NIMO validation and demonstration, 2026](https://angeo.copernicus.org/articles/44/303/2026/) provides an independent-model validation precedent, not a measured error for this GEMINI branch.
- [Scalable Mesh Coupling for Atmospheric Wave Simulation](https://arxiv.org/html/2603.02971v1) couples neutral-state perturbations and plasma feedback in compatible MAGIC/GEMINI formulations. That precedent does not provide the missing exchange law or compressible neutral state for CVTWIN's existing incompressible MHD representation.
- [SciPy exact binomial interval documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats._result_classes.BinomTestResult.proportion_ci.html) documents the Clopper-Pearson method. The new UQ checker uses a one-sided exact bound; empirical coverage alone cannot close R13.
- [GitHub workflow job API](https://docs.github.com/en/rest/actions/workflow-jobs) supplies attempt-specific job and step metadata. The read-only collector requires success of named steps on the exact candidate, avoiding acceptance of stale or conditionally skipped required work.

These updates support the implemented diagnostics and evidence checks. They do not close the outstanding scientific gates.
