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

integer :: i,j,N

if(lid < 1) error stop "MPI image count must be positive"
if(lid == 1) then
  lid2 = 1
  lid3 = 1
  return
endif

if (lx3all==1) then
!  !! 2D simulation in x2, SWAP x2 to x3
!  lid3 = gcd(lid, lx2all)
!  lid2 = 1
  lid3 = 1
  lid2 =  gcd(lid, lx2all)
elseif (lx2all==1) then
  !! 2D simulation in x3
  lid3 = gcd(lid, lx3all)
  lid2 = 1
else
  ! print *, 'lx2all,lx3all,lid', lx2all,lx3all,lid
  !! 3D simulation
  !! in 3D sims, must have at least 2 images on each axis to avoid MPI/MUMPS memory / vader errors
  lid2 = 0
  lid3 = huge(0)
  N = 0

  do i = gcd(lid, lx2all), 1, -1
    do j = gcd(lid, lx3all), 1, -1
      if (i*j /= lid) cycle
      if (i*j < N) cycle
      if (modulo(lx2all, i) /= 0) cycle
      if (modulo(lx3all, j) /= 0) cycle
      ! print *, "trying ", i, j,  N
      if (lx2all / i == 1 .or. lx3all / j == 1) cycle
      if (abs(i-j) > abs(lid2-lid3)) cycle

      N = i*j
      lid2 = i
      lid3 = j
      ! print *, "found ", i,j, " trying to find better"
    enddo
  enddo
  if (N==0) then
    write(stderr,'(A,3I4)') "grid_auto: no non-singular factorization found for MPI image partition. lx2all,lx3all,lid: ",&
      lx2all,lx3all,lid
    error stop
  endif
end if

call check_partition(lx2all, lx3all, lid, lid2, lid3)

end subroutine grid_auto


subroutine check_partition(lx2all, lx3all, lid, lid2, lid3)
!! checks grid partitioning for MPI
integer, intent(in) :: lx2all,lx3all,lid, lid2, lid3

character(6) :: s1,s2
integer :: lx2, lx3

if (lid < 1 .or. lid2 < 1 .or. lid3 < 1) error stop "MPI image count must be positive"

if (lid2*lid3 /= lid) then
  write(s1, '(I6)') lid
  write(s2, '(I6)') lid2*lid3
  error stop "autogrid:check_partition MPI image count " // s1 // " not a factor of x2*x3 " // s2
endif

if (lx2all > 1 .and. lx3all > 1) then
  lx2 = lx2all / lid2
  lx3 = lx3all / lid3
  if(lx2 == 1 .or. lx3 == 1) then
    write(stderr,'(A,/,A,8I4)') "ERROR: 3-D grid cannot have singular MPI axis:","lx2,lx3,lx2all,lx3all,lid2,lid3,lid: ", &
      lx2, lx3,lx2all,lx3all,lid2,lid3,lid
    error stop
  endif
endif

if (lx2all > 1) then
  if (lid2 > lx2all) error stop "lid2 cannot be greater than lx2"

  if (modulo(lx2all, lid2) /= 0) then
    write(s1, '(I6)') lid2
    write(s2, '(I6)') lx2all
    error stop 'autogrid:check_partition MPI x2 image count ' // s1 // ' not a factor of lx2 ' // s2
  endif
endif

if (lx3all > 1) then
  if (lid3 > lx3all) error stop "lid3 cannot be greater than lx3"

  if (modulo(lx3all, lid3) /= 0) then
    write(s1, '(I6)') lid3
    write(s2, '(I6)') lx3all
    error stop 'autogrid:check_partition MPI x3 image count ' // s1 // ' not a factor of lx3 ' // s2
  endif
endif

if (modulo(lid, lid2) /= 0) then
  write(s1, '(I6)') lid2
  write(s2, '(I6)') lid
  error stop "autogrid:check_partition MPI x2 image count " // s1 // " not a factor of lid " // s2
endif

if (modulo(lid, lid3) /= 0)  then
  write(s1, '(I6)') lid3
  write(s2, '(I6)') lid
  error stop "autogrid:check_partition MPI x3 image count " // s1 // " not a factor of lid " // s2
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

!> divide by 2 to ensure MPI partition has at least 2 per axis--for 2D and 3D
if (lx3 == 1) then
  max_mpi = max_gcd(lx2/2, max_cpu)
elseif (lx2 == 1) then
  max_mpi = max_gcd(lx3/2, max_cpu)
else
  max_mpi = max_gcd2(lx2/2, lx3/2, max_cpu)
end if

end function max_mpi


integer function max_gcd(L, M)
!! find the Greatest Common Factor to evenly partition the simulation grid
!! Output range is [M, 1]

integer, intent(in) :: L,M
integer :: i

if (M < 1) error stop "autogrid:max_gcd CPU count must be at least one"

max_gcd = 1
do i = M, 2, -1
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
integer :: f2, f3, i,j, t2, t3

if (M < 1) error stop "autogrid:max_gcd2 CPU count must be at least one"

max_gcd2 = 1
t2 = 1
t3 = huge(0)
x2 : do i = M,2,-1
  x3 : do j = M,2,-1
    f2 = max_gcd(lx2, i)
    f3 = max_gcd(lx3, j)
    if (M < f2 * f3) cycle x3 !< too big for CPU count
    if (f2 * f3 < max_gcd2) exit x3 !< too small for max CPU count
    if (lx2 / i == 1) cycle x2
    if (lx3 / j == 1) cycle x3
    if (abs(f2-f3) > abs(t2-t3)) cycle x3

    t2 = f2
    t3 = f3
    ! print *, "lid2, lid3: ", f2, f3
    max_gcd2 = f2 * f3
  end do x3
  if (i*M < max_gcd2) exit x2
end do x2

! print *, "final: lid2, lid3: ", t2, t3

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
