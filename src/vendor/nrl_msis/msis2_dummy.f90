module msis_calc

use, intrinsic :: iso_fortran_env, only : real32, stderr=>error_unit
implicit none (type, external)

contains

subroutine msiscalc(day,utsec,z,lat,lon,sfluxavg,sflux,ap,tn,dn,tex)

class(*) :: day,utsec,z,lat,lon,sfluxavg,sflux,ap(7), tn, dn(10), tex

write(stderr,*) 'ERROR: to use MSIS 2.x requires "cmake -Dmsis2=yes"'
error stop 20

end subroutine msiscalc


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

write(stderr,*) 'ERROR: to use MSIS 2.x requires "cmake -B build -Dmsis2=yes"'
error stop 20

end subroutine msisinit

end module msis_init
