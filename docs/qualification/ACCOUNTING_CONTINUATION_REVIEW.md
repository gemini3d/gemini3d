# Qualification review continuation, 10 October 2026

The existing acceptance budgets, native comparator, and reference archives remain
unchanged. This work adds evidence checks and optional diagnostics. No scientific
gate is closed by these component changes.

Deleting the same middle time step from every native ledger previously left all
four checkers green. The checkers now require contiguous float64 clocks, a zero
start, and the final `tdur` from the retained input configuration. Cross-ledger
coverage remains required. Clock comparisons allow only 64 float64 ULPs for
clock accumulation; no mass, charge, transport, source, or energy budget changes.

The temperature conversion after compression enforces a 100 K lower floor but
previously recorded no addition for this stage. Optional diagnostics now record
stage 3, using actual species density for physical internal-energy capacity.
The source checker requires all seven species at each of the three standard
energy floor stages. A native activation oracle checks a known 50-to-100 K
addition over 183 m3, excludes a null cell and ghosts, and compares the complete
temperature/internal-energy arrays with diagnostics enabled and disabled.

This does not complete R04. Individual chemistry/collision channels, coupled
kinetic/internal-energy exchanges, density/velocity cleanup and conversion
effects, independently measured conductive boundary power, and full coupled
space/time refinement remain distinct obligations. Endpoint reservoir changes
are not reclassified as conductive flux. Numerical regression does not supply
R01 maintainer acceptance or R07 independent matched-event physics validation.
