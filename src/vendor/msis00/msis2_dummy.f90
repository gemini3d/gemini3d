module msis_calc

use, intrinsic :: iso_fortran_env, only : real64, real32, stderr=>error_unit
implicit none (type, external)

interface msiscalc
  procedure msiscalc_r32, msiscalc_r64
end interface msiscalc

contains

subroutine msiscalc_r32(day,utsec,z,lat,lon,sfluxavg,sflux,ap,tn,dn,tex)

real(real32) :: day,utsec,z,lat,lon,sfluxavg,sflux,ap(7)
real(real32) :: tn, dn(10)
real(real32) :: tex

write(stderr,*) 'ERROR: to use MSIS 2.0 requires "cmake -B build -Dmsis20=yes"'
error stop 20

end subroutine msiscalc_r32


subroutine msiscalc_r64(day,utsec,z,lat,lon,sfluxavg,sflux,ap,tn,dn,tex)

real(real64) :: day,utsec,z,lat,lon,sfluxavg,sflux,ap(7)
real(real64) :: tn, dn(10)
real(real64) :: tex

write(stderr,*) 'ERROR: to use MSIS 2.0 requires "cmake -B build -Dmsis20=yes"'
error stop 20

end subroutine msiscalc_r64

end module msis_calc


module msis_init

use, intrinsic :: iso_fortran_env, only : real32, stderr=>error_unit
implicit none (type, external)

contains

subroutine msisinit(parmpath,parmfile,iun,switch_gfn,switch_legacy, &
                      lzalt_type,lspec_select,lmass_include,lN2_msis00)

integer, parameter :: nspec=11, maxnbf=512

character(*), optional :: parmpath, parmfile
integer, optional :: iun
logical, optional :: switch_gfn(0:maxnbf-1)
real(real32), optional :: switch_legacy(25)
logical, optional :: lzalt_type,lspec_select(nspec-1), lmass_include(nspec-1), lN2_msis00

write(stderr,*) 'ERROR: to use MSIS 2.0 requires "cmake -B build -Dmsis20=yes"'
error stop 20

end subroutine msisinit

end module msis_init
