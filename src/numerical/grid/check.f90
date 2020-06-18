submodule (grid) grid_checker

implicit none (type, external)

contains

module procedure grid_check
!! Help avoid confusing errors or bad simulations
!! Must be called after both grid size and mpi gridding have been established...

real(wp) :: tol_inc

!> check correct number of MPI images along x2 and x3
if (lx2all > 1) then
  if (modulo(lx2all, lid2) /= 0) then
    write(stderr,'(/,A,I6,A,I6,/)') 'ERROR:grid_size_root: Number of MPI images along x2', lid2, &
                                    ' not a factor of lx2all: ', lx2all
    error stop
  endif
endif

if (lx3all > 1) then
  if (modulo(lx3all, lid3) /= 0) then
    write(stderr,'(/,A,I6,A,I6,/)') 'ERROR:grid_size_root: Number of MPI images along x3', lid3, &
                                    ' not a factor of lx3all: ', lx3all
    error stop
  endif
endif

!> check for monotonic increasing

!tol_inc = 0.1
tol_inc = 1e-6     !< MZ - patched to function with dipole grids
call is_monotonic_increasing(x%x1, tol_inc, 'x1')
call is_monotonic_increasing(x%x1i, tol_inc, 'x1i')
call is_monotonic_increasing(x%x2, tol_inc, 'x2')
call is_monotonic_increasing(x%x2i, tol_inc, 'x2i')
call is_monotonic_increasing(x%x3, tol_inc, 'x3')
call is_monotonic_increasing(x%x3i, tol_inc, 'x3i')

end procedure grid_check


subroutine is_monotonic_increasing(A, tol, name)

real(wp), intent(in) :: A(:), tol
character(*), intent(in) :: name

integer :: i

!> probably better than creating a temporary array
do i = 2,size(A)
  if(A(i) - A(i-1) < tol) then
    write(stderr,*) 'ERROR:is_monotonic_increasing: ',name,' is not monotonic increasing'
    error stop
  endif
end do


end subroutine is_monotonic_increasing

end submodule grid_checker
