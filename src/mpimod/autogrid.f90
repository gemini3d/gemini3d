module autogrid

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
implicit none (type, external)

private
public :: grid_auto, test_process_number, max_mpi

contains


subroutine grid_auto(lx2all,lx3all, lid, lid2, lid3)
!! Automatically determine the PROCESS GRID
!! sets value of lid2,lid3 globally
!!
!! lid: total number of MPI processes, set by mpi_comm_size (from CLI mpiexec -n).
!!      Consider lid read-only or weird behavior can result (deadlock, crash)

integer, intent(in) :: lx2all,lx3all,lid
integer, intent(out) :: lid2, lid3

integer :: inds(2), i,j,N

if (lx3all==1) then
  !! 2D simulation in x2, SWAP x2 to x3
  lid3 = gcd(lid, lx2all)
  lid2 = 1
elseif (lx2all==1) then
  !! 2D simulation in x3
  lid3 = gcd(lid, lx3all)
  lid2 = 1
else
  !! 3D simulation
  lid2 = 1
  lid3 = 1
  N = 1
  do i = gcd(lid, lx2all),1,-1
    do j = gcd(lid, lx3all),1,-1
      if (i*j > lid) cycle
      if (i*j > N) then
        N = i*j
        lid2 = i
        lid3 = j
      endif
    enddo
  enddo
end if

call check_partition(lx2all, lx3all, lid, lid2, lid3)

end subroutine grid_auto


pure subroutine check_partition(lx2all, lx3all, lid, lid2, lid3)
!! checks grid partitioning for MPI
integer, intent(in) :: lx2all,lx3all,lid, lid2, lid3
character(6) :: s1,s2

if (lx2all > 1 .and. modulo(lx2all, lid2) /= 0) then
  write(s1, '(I6)') lid2
  write(s2, '(I6)') lx2all
  error stop 'autogrid:grid_auto MPI x2 image count ' // s1 // ' not a factor of lx2 ' // s2
endif
if (lx3all > 1 .and. modulo(lx3all, lid3) /= 0) then
  write(s1, '(I6)') lid3
  write(s2, '(I6)') lx3all
  error stop 'autogrid:grid_auto MPI x3 image count ' // s1 // ' not a factor of lx3 ' // s2
endif

if (lx3all > 1 .and. lid3 > lx3all) error stop "lid3 cannot be greater than lx3"
if (lx2all > 1 .and. lid2 > lx2all) error stop "lid2 cannot be greater than lx2"

if (modulo(lid, lid2) /= 0) then
  write(s1, '(I6)') lid2
  write(s2, '(I6)') lid
  error stop "autogrid:grid_auto MPI x2 image count " // s1 // " not a factor of lid " // s2
endif

if (modulo(lid, lid3) /= 0)  then
  write(s1, '(I6)') lid3
  write(s2, '(I6)') lid
  error stop "autogrid:grid_auto MPI x3 image count " // s1 // " not a factor of lid " // s2
endif

if (lid2*lid3 /= lid) then
  write(s1, '(I6)') lid
  write(s2, '(I6)') lid2*lid3
  error stop "autogrid:grid_auto MPI image count " // s1 // " not a factor of x2*x3 " // s2
endif

end subroutine check_partition


integer function gcd(a, b)

integer, intent(in) :: a, b
integer :: x,y,z

if (a < 1 .or. b < 1) error stop "autogrid:gcd positive integers only"

x = a
y = b
z = modulo(x, y)
do while (z /= 0)
  x = y
  y = z
  z = modulo(x, y)
end do
gcd = y

end function gcd


integer function max_mpi(lx2, lx3, max_cpu)
!! goal is to find the highest x2 + x3 to maximum CPU core count
integer, intent(in) :: lx2, lx3, max_cpu

if (lx3 == 1) then
  max_mpi = max_gcd(lx2, max_cpu)
elseif (lx2 == 1) then
  max_mpi = max_gcd(lx3, max_cpu)
else
  max_mpi = max_gcd2(lx2, lx3, max_cpu)
end if

end function max_mpi


integer function max_gcd(L, M)
!! find the Greatest Common Factor to evenly partition the simulation grid
!! Output range is [M, 1]
integer, intent(in) :: L,M
integer :: i

if (M < 1) error stop "autogrid:max_gcd CPU count must be at least one"

max_gcd = 1
do i = M, 1, -1
  max_gcd = max(gcd(L, i), max_gcd)
  if (i < max_gcd) exit
end do

end function max_gcd


integer function max_gcd2(lx2, lx3, M)
!! find the Greatest Common Factor to evenly partition the simulation grid
!! Output range is [M, 1]
!!
!!    1. find factors of each dimension
!!    2. choose partition that yields highest CPU count usage

integer, intent(in) :: lx2, lx3, M
integer :: f2, f3, i,j

if (M < 1) error stop "autogrid:max_gcd2 CPU count must be at least one"

max_gcd2 = 1
do i = M,1,-1
  do j = M,1,-1
    f2 = max_gcd(lx2, i)
    f3 = max_gcd(lx3, j)
    if ((M >= f2 * f3) .and. (f2 * f3 > max_gcd2)) then
      ! print *, f2, f3
      max_gcd2 = f2 * f3
    endif
  end do
end do

end function max_gcd2

logical function test_process_number(N, lx2all, lx3all, rx2, rx3) result (ok)
!! this is only for testing
integer, intent(in) :: N, rx2, rx3, lx2all, lx3all

integer :: lid, lid2, lid3

lid = N
call grid_auto(lx2all, lx3all, lid, lid2, lid3)

ok = (lid2 == rx2 .and. lid3 == rx3)

if (.not. ok) then
  write(stderr,'(A,I0,A1,I0,A1,I0,A1,I0)') 'failed: lx2all,lx3all,lid,Nimg: ', lx2all,' ', lx3all,' ', lid,' ',N
  write(stderr,'(A,I0,A1,I0,A,I0,A1,I0)') 'expected lid2,lid3 ', rx2,' ', rx3, ' but got: ', lid2,' ', lid3
end if

end function test_process_number

end module autogrid
