submodule (compare_h5) compare_grid_h5

implicit none (type, external)

contains

module procedure check_grid

type(hdf5_file) :: hnew, href

character(:), allocatable :: new_file, ref_file

!! read grid functions are being refactored, so for now we just iterate over a priori list of variables

character(5), parameter :: var1(*) = [character(5) :: &
"/x1", "/x1i", "/dx1b", "/dx1h", &
"/x2", "/x2i", "/dx2b", "/dx2h", &
"/x3", "/x3i", "/dx3b", "/dx3h"]

character(2), parameter :: var2(*) = [character(2) :: '/I']

character(8), parameter :: var3(*) = [character(8) :: &
"/h1", "/h2", "/h3", &
'/h1x1i', '/h2x1i', '/h3x1i', &
'/h1x2i', '/h2x2i', '/h3x2i', &
'/h1x3i', '/h2x3i', '/h3x3i', &
'/gx1', '/gx2', '/gx3', &
'/alt', '/glat', '/glon', &
'/Bmag', '/nullpts', &
'/r', '/theta', '/phi']

character(7), parameter :: var4(*) = [character(7) :: &
'/e1', '/e2', '/e3', &
'/er', '/etheta', '/ephi']

real, allocatable :: A1(:), B1(:), A2(:,:), B2(:,:), A3(:,:,:), B3(:,:,:), A4(:,:,:,:),B4(:,:,:,:)
integer(int64), allocatable :: dims(:)

integer :: i, bad

new_file = new_path // "/inputs/simgrid.h5"
ref_file = ref_path // "/inputs/simgrid.h5"

if(.not. is_file(new_file)) error stop "compare:check_grid: new grid file not found: " // new_file
if(.not. is_file(ref_file)) error stop "compare:check_grid: ref grid file not found: " // ref_file

bad = 0

call hnew%open(new_file, action='r')
call href%open(ref_file, action='r')

do i = 1,size(var1)
  call href%shape(var1(i), dims)
  allocate(A1(dims(1)), B1(dims(1)))

  if(P%debug) print '(A,1X,I0,1X,I0)', var1(i), dims, shape(A1)

  call hnew%read(var1(i), A1)
  call href%read(var1(i), B1)

  if(.not.all(isclose(B1, A1, real(rtol), real(atol)))) then

    bad = bad + 1

    write(stderr,'(A,/,A,ES12.3,A,2ES12.3,A,2ES12.3)') "MISMATCH: " // file_name(new_file) // " " // var1(i), &
    ' max diff:', maxval(abs(B1 - A1)), ' max & min ref:', maxval(B1), minval(B1), ' max & min new:', maxval(A1), minval(A1)
  endif

  deallocate(A1, B1)
end do


do i = 1,size(var2)
  call href%shape(var2(i), dims)
  allocate(A2(dims(1), dims(2)), B2(dims(1), dims(2)))

  if(P%debug) print '(A,1X,2I6,1X,2I6)', var2(i), dims, shape(A2)

  call hnew%read(var2(i), A2)
  call href%read(var2(i), B2)

  if(.not.all(isclose(B2, A2, real(rtol), real(atol)))) then

    bad = bad + 1

    write(stderr,'(A,/,A,ES12.3,A,2ES12.3,A,2ES12.3)') "MISMATCH: " // file_name(new_file) // " " // var2(i), &
    ' max diff:', maxval(abs(B2 - A2)), ' max & min ref:', maxval(B2), minval(B2), ' max & min new:', maxval(A2), minval(A2)
  endif

  deallocate(A2, B2)
end do


do i = 1,size(var3)
  call href%shape(var3(i), dims)
  allocate(A3(dims(1), dims(2), dims(3)), B3(dims(1), dims(2), dims(3)))

  if(P%debug) print '(A,1X,3I6,1X,3I6)', var3(i), dims, shape(A3)

  call hnew%read(var3(i), A3)
  call href%read(var3(i), B3)

  if(.not.all(isclose(B3, A3, real(rtol), real(atol)))) then

    bad = bad + 1

    write(stderr,'(A,/,A,ES12.3,A,2ES12.3,A,2ES12.3)') "MISMATCH: " // file_name(new_file) // " " // var3(i), &
    ' max diff:', maxval(abs(B3 - A3)), ' max & min ref:', maxval(B3), minval(B3), ' max & min new:', maxval(A3), minval(A3)
  endif

  deallocate(A3, B3)
end do


do i = 1,size(var4)
  call href%shape(var4(i), dims)
  allocate(A4(dims(1), dims(2), dims(3), dims(4)), B4(dims(1), dims(2), dims(3), dims(4)))

  if(P%debug) print '(A,1X,4I6,1X,4I6)', var4(i), dims, shape(A4)

  call hnew%read(var4(i), A4)
  call href%read(var4(i), B4)

  if(.not.all(isclose(B4, A4, real(rtol), real(atol)))) then

    bad = bad + 1

    write(stderr,'(A,/,A,ES12.3,A,2ES12.3,A,2ES12.3)') "MISMATCH: " // file_name(new_file) // " " // var4(i), &
    ' max diff:', maxval(abs(B4 - A4)), ' max & min ref:', maxval(B4), minval(B4), ' max & min new:', maxval(A4), minval(A4)
  endif

  deallocate(A4, B4)
end do


call hnew%close()
call href%close()

check_grid = bad == 0


end procedure check_grid

end submodule compare_grid_h5
