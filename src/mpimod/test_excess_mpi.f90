module mpimod
!! dummy parent module to avoid needing MPI to test range of parameters
use, intrinsic:: iso_fortran_env, only: stderr=>error_unit

implicit none (type, external)

integer :: lid, lid2,lid3
integer :: myid=0, myid2, myid3

interface ! mpi_grid.f90
module subroutine mpigrid(lx2all,lx3all)
integer, intent(in) :: lx2all,lx3all
end subroutine mpigrid

module subroutine mpi_manualgrid(lx2all,lx3all,lid2in,lid3in)
integer, intent(in) :: lx2all,lx3all, lid2in,lid3in
end subroutine mpi_manualgrid

module integer function grid2id(i2,i3)
integer, intent(in) :: i2,i3
end function grid2id

module function ID2grid(ID)
integer, dimension(2) :: ID2grid
integer, intent(in) :: ID
end function id2grid
end interface


end module mpimod


program test_excess_mpi

use mpimod
implicit none (type, external)


integer, parameter :: N(*) = [1,2,3,4,5,6,8,18,28,64] !< number of fake CPU

!> 2D
call checker(N, lx2all=1, lx3all=40, &
  rx2=[1,1,1,1,1,1,1,1,1,1], &
  rx3=[1,2,3,4,5,6,8,18,28,40])
call checker(N, lx2all=40, lx3all=1, &
  rx2=[1,1,1,1,1,1,1,1,1,1], &
  rx3=[1,2,3,4,5,6,8,18,28,40])

!> 3D
!> NOTE: these are not the most efficient CPU use possible, we
!> can make the auto-grid algorithm more CPU-effective.
call checker(N, lx2all=8, lx3all=12, &
  rx2=[1,1,1,2,2,2,2,4,4,4], &
  rx3=[1,2,3,2,2,3,3,3,3,3])
call checker(N, lx2all=12, lx3all=8, &
  rx2=[1,1,1,2,2,2,2,2,2,2], &
  rx3=[1,2,2,2,2,2,4,4,4,4])

call checker(N, lx2all=6, lx3all=10, &
  rx2=[1,1,1,1,1,1,1,2,2,2], &
  rx3=[1,2,2,2,5,5,5,5,5,5])
call checker(N, lx2all=10, lx3all=6, &
  rx2=[1,1,1,1,1,2,2,2,2,2], &
  rx3=[1,2,3,3,3,3,3,3,3,3])

call checker(N, lx2all=4, lx3all=4, &
  rx2=[1,1,1,2,2,2,2,2,2,2], &
  rx3=[1,2,2,2,2,2,2,2,2,2])


print *, 'OK: auto process grid'


contains


subroutine checker(N, lx2all, lx3all, rx2, rx3)

integer, intent(in) :: N(:), rx2(:), rx3(:), lx2all, lx3all
integer :: i

do i = 1,size(N)
  lid = N(i)
  call mpigrid(lx2all,lx3all)
  if (lid2 /= rx2(i) .or. lid3 /= rx3(i)) then
    write(stderr,'(A,5I4)') 'failed: lx2all,lx3all,lid,N:',lx2all,lx3all,lid,N(i)
    write(stderr,*) 'expected lid2,lid3', rx2(i), rx3(i), 'but got:',lid2,lid3
    error stop
  end if
end do


end subroutine checker

end program
