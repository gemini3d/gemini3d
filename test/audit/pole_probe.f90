! Audit addition 2026-09-16. Apache-2.0.
program audit_pole
use geomagnetic, only: set_magnetic_pole
implicit none
call set_magnetic_pole(2026)
end program
