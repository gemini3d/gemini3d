module reader
!! simple file reading procedures

use, intrinsic :: iso_fortran_env, only: stderr=>error_unit
use phys_consts, only: wp, debug

implicit none
private
public :: get_simsize3, get_simsize2, get_grid2, get_Efield, get_precip, get_neutral2, get_neutral3

interface ! reader_{raw,hdf5,nc4}.f90
module subroutine get_simsize2(path, llon, llat)
character(*), intent(in) :: path
integer, intent(out) :: llon, llat
end subroutine get_simsize2

module subroutine get_simsize3(path, lx1, lx2all, lx3all)
character(*), intent(in) :: path
integer, intent(out) :: lx1, lx2all
integer, intent(out), optional :: lx3all
end subroutine get_simsize3

module subroutine get_grid2(path, mlonp, mlatp)
character(*), intent(in) :: path
real(wp), dimension(:), intent(out) :: mlonp, mlatp
end subroutine get_grid2

module subroutine get_Efield(path, flagdirich,E0xp,E0yp,Vminx1p,Vmaxx1p,Vminx2pslice,Vmaxx2pslice,Vminx3pslice,Vmaxx3pslice)
character(*), intent(in) :: path
integer, intent(out) :: flagdirich
real(wp), dimension(:,:), intent(out) :: E0xp,E0yp,Vminx1p,Vmaxx1p
real(wp), dimension(:), intent(out) :: Vminx2pslice,Vmaxx2pslice,Vminx3pslice,Vmaxx3pslice
end subroutine get_Efield

module subroutine get_precip(path, Qp, E0p)
character(*), intent(in) :: path
real(wp), dimension(:,:), intent(out) :: Qp, E0p
end subroutine get_precip

module subroutine get_neutral2(path, dnO,dnN2,dnO2,dvnrho,dvnz,dTn)
character(*), intent(in) :: path
real(wp), dimension(:,:,:), intent(out) :: dnO,dnN2,dnO2,dvnrho,dvnz,dTn
end subroutine get_neutral2

module subroutine get_neutral3(path, dnOall,dnN2all,dnO2all,dvnxall,dvnrhoall,dvnzall,dTnall)
character(*), intent(in) :: path
real(wp), dimension(:,:,:), intent(out) :: dnOall,dnN2all,dnO2all,dvnxall,dvnrhoall,dvnzall,dTnall
end subroutine get_neutral3
end interface

contains


end module reader