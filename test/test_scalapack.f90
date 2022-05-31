program test_scalapack
! minimal Scalapack demo
use, intrinsic :: iso_fortran_env, only: wp => real64

implicit none (type, external)

integer :: ictxt, myid, nprocs, mycol, myrow, npcol, nprow
real(wp) :: eps
real(wp), external :: pdlamch

external :: blacs_pinfo, blacs_get, blacs_gridinit, blacs_gridinfo, blacs_gridexit, blacs_exit

! arbitrary test parameters
npcol = 2
nprow = 2

call blacs_pinfo(myid, nprocs)
call blacs_get(-1, 0, ictxt)
call blacs_gridinit(ictxt, "C", nprocs, 1)

call blacs_gridinfo(ictxt, nprow, npcol, myrow, mycol)

eps = pdlamch(ictxt, 'E')

if(myrow == mycol) print '(A, F10.6)', "OK: Scalapack Fortran  eps=", eps

call blacs_gridexit(ictxt)
call blacs_exit(0)

end program
