submodule (grid:grid_mpi) grid_checker_mpi

implicit none (type, external)

contains

module procedure grid_check
  !! Help avoid confusing errors or bad simulations
  !! Must be called after both grid size and mpi gridding have been established...

  real(wp) :: tol_inc, tol_inc_big, tol_big

  !> check correct number of MPI images along x2 and x3
  if (lx2all > 1) then
    if (modulo(lx2all, mpi_cfg%lid2) /= 0) then
      write(stderr,'(/,A,I6,A,I6,/)') 'ERROR:grid_size_root: Number of MPI images along x2', mpi_cfg%lid2, &
                                      ' not a factor of lx2all: ', lx2all
      error stop
    endif
  endif

  if (lx3all > 1) then
    if (modulo(lx3all, mpi_cfg%lid3) /= 0) then
      write(stderr,'(/,A,I6,A,I6,/)') 'ERROR:grid_size_root: Number of MPI images along x3', mpi_cfg%lid3, &
                                      ' not a factor of lx3all: ', lx3all
      error stop
    endif
  endif

  !> check for monotonic increasing

  !tol_inc = 0.1
  tol_inc = 1e-6     !< MZ - patched to function with dipole grids
  tol_inc_big = 1d6
  tol_big = 1e9
  call is_monotonic_increasing(x%x1, tol_inc, tol_inc_big, tol_big, 'x1')
  call is_monotonic_increasing(x%x1i, tol_inc, tol_inc_big, tol_big, 'x1i')
  call is_monotonic_increasing(x%x2, tol_inc, tol_inc_big, tol_big, 'x2')
  call is_monotonic_increasing(x%x2i, tol_inc, tol_inc_big, tol_big, 'x2i')
  call is_monotonic_increasing(x%x3, tol_inc, tol_inc_big, tol_big, 'x3')
  call is_monotonic_increasing(x%x3i, tol_inc, tol_inc_big, tol_big, 'x3i')

  !> geo lat/lon

  if (any(x%glat > 90) .or. any(x%glat < -90)) error stop 'geo latitude outside expected range'

  if (any(x%glon < 0) .or. any(x%glon > 360)) error stop 'geo longitude outside expected range'
end procedure grid_check


subroutine is_monotonic_increasing(A, tol_inc, tol_inc_big, tol_big, name)
  real(wp), intent(in) :: A(:), tol_inc, tol_inc_big, tol_big
  character(*), intent(in) :: name

  integer :: i

  !> probably better than creating a temporary array
  do i = 2,size(A)
    if(A(i) - A(i-1) < tol_inc) then
      error stop 'grid_check: ' // name // ' is not monotonic increasing'
    endif
    if(A(i) - A(i-1) > tol_inc_big) then
      error stop 'grid_check: ' // name // ' has too large derivative'
    endif
  end do

  !> not too big value
  if (any(abs(A) > tol_big)) then
    write(stderr,*) 'ERROR:grid_check: ',name,' has too large values: ', maxval(A)
    error stop
  end if
end subroutine is_monotonic_increasing

end submodule grid_checker_mpi
