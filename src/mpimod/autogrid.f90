module autogrid

implicit none (type, external)

private
public :: grid_auto, ID2grid

contains


pure subroutine grid_auto(lx2all,lx3all, lid, myid, lid2, lid3, myid2, myid3)
!! Automatically determine the PROCESS GRID
!! sets value of lid2,lid3 globally
!!
!! lid: total number of MPI processes, set by mpi_comm_size (from CLI mpiexec -n).
!!      Consider lid read-only or weird behavior can result (deadlock, crash)

integer, intent(in) :: lx2all,lx3all,lid, myid
integer, intent(out) :: lid2, lid3, myid2, myid3

character(:), allocatable :: s1,s2
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

  if (modulo(lx2all, lid2) /= 0) then
    write(s1, '(I0)') lid2
    write(s2, '(I0)') lx2all
    error stop 'autogrid:grid_auto MPI x2 image count ' // s1 // ' not a factor of lx2 ' // s2
  endif
  if (modulo(lx3all, lid3) /= 0) then
    write(s1, '(I0)') lid3
    write(s2, '(I0)') lx3all
    error stop 'autogrid:grid_auto MPI x3 image count ' // s1 // ' not a factor of lx3 ' // s2
  endif
end if

!> checks
if (lx3all > 1 .and. lid3 > lx3all) error stop "lid3 cannot be greater than lx3"
if (lx2all > 1 .and. lid2 > lx2all) error stop "lid2 cannot be greater than lx2"

if (modulo(lid, lid2) /= 0) then
  write(s1, '(I0)') lid2
  write(s2, '(I0)') lid
  error stop "autogrid:grid_auto MPI x2 image count " // s1 // " not a factor of lid " // s2
endif

if (modulo(lid, lid3) /= 0)  then
  write(s1, '(I0)') lid3
  write(s2, '(I0)') lid
  error stop "autogrid:grid_auto MPI x3 image count " // s1 // " not a factor of lid " // s2
endif

if (lid2*lid3 /= lid) then
  write(s1, '(I0)') lid
  write(s2, '(I0)') lid2*lid3
  error stop "autogrid:grid_auto MPI image count " // s1 // " not a factor of x2*x3 " // s2
endif

!> THIS PROCESS' LOCATION ON THE GRID
inds = ID2grid(myid, lid2)
myid2 = inds(1)
myid3 = inds(2)

end subroutine grid_auto


pure integer function gcd(a, b)

integer, intent(in) :: a, b
integer :: x,y,z

if (a < 1 .or. b < 1) error stop "positive integers only"

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


pure function ID2grid(ID, lid2)
!! COMPUTES GRID LOCATION FROM A PROCESS ID
integer, dimension(2) :: ID2grid
integer, intent(in) :: ID, lid2

ID2grid(2) = ID / lid2
!! x3 index into process grid
ID2grid(1) = ID - ID2grid(2) * lid2
!! x2 index into process grid

end function ID2grid


end module autogrid
