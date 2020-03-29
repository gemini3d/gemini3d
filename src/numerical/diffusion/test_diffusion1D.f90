program test_diffusion1d
!! Need program statement for FORD
!! Solve a time-dependent heat equation in 1D.  See GEMINI-docs repo for
!! a description of the specific problem solved here

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use phys_consts, only : wp,pi
use PDEparabolic, only : backEuler1D,TRBDF21D
use h5fortran, only : hdf5_file

implicit none (type, external)

type(hdf5_file) :: h5f
integer, parameter :: npts=256,lt=20*5

character(256) :: argv
character(:), allocatable :: outfile
character(4) :: ic

real(wp), parameter :: abs_tol = 0.05_wp

real(wp), dimension(npts) :: v1,dx1i
real(wp), dimension(-1:npts+2) :: x1
real(wp), dimension(-1:npts+2, 0:lt) :: TsEuler,TsBDF2,Tstrue
real(wp), dimension(npts) :: lambda,A,B,C,D,E
real(wp), dimension(npts+1) :: x1i
real(wp), dimension(0:npts+2) :: dx1
integer :: lx1,i,ix1, flagdirichTop, flagdirichBottom
real(wp) :: t(lt), dt
real(wp) :: Tsminx1,Tsmaxx1

real(wp), dimension(npts-2) :: errorEuler,errorBDF2
logical :: failed = .false.

call get_command_argument(1, argv, status=i)
if(i==0) then
  read(argv,'(I2)') flagdirichBottom
else
  flagdirichBottom = 1
endif
call get_command_argument(2, argv, status=i)
if(i==0) then
  read(argv,'(I2)') flagdirichTop
else
  flagdirichTop = 1
endif

call get_command_argument(3, argv, status=i)
if(i==0) then
  outfile = trim(argv)
else
  outfile = 'test_diffusion1d.h5'
endif

!! create a grid for the calculation
x1=[ (real(ix1-1) / real(npts-1), ix1=-1,npts+2) ]
lx1=npts   !exclude ghost cells in count
dx1=x1(0:lx1+2)-x1(-1:lx1+1)
x1i(1:lx1+1)=0.5*(x1(0:lx1)+x1(1:lx1+1))
dx1i=x1i(2:lx1+1)-x1i(1:lx1)


!! write the time, space length adn spatial grid to a file
print *,'writing ',outfile
call h5f%initialize(outfile, status='replace', action='write')
call h5f%write('/lx1', lx1)
call h5f%write('/x1', x1)
call h5f%write('/flagdirichTop', flagdirichTop)
call h5f%write('/flagdirichBottom', flagdirichBottom)

!! initial conditions
TsEuler(-1:lx1+2, 0) = sin(2*pi*x1(-1:lx1+2)) + sin(8*pi*x1(-1:lx1+2))
TsBDF2(:, 0) = TsEuler(:, 0)
lambda(:) = 1     !< thermal conductivity


!! typical diffusion time, make our time step a fraction of this
dt = 0.05 * 1/8.0_wp**2 / pi**2 / maxval(lambda)
t(1) = 0
do i = 2,lt
  t(i) = t(i-1) + dt
end do
call h5f%write('/t', t)

print '(A3, 3A12)', 'i','dt  ', 'Error Euler', 'Error TBDF'
!! time interations

do i = 1,lt
  !! compute analytical solution to compare
  Tstrue(1:lx1, i) = exp(-4*pi**2*lambda* (t(i) + dt)) * sin(2*pi*x1(1:lx1)) + &
                     exp(-64*pi**2*lambda* (t(i) + dt)) * sin(8*pi*x1(1:lx1))
end do

do i=1,lt
  write(ic, '(I4.4)') i
  !boundary values
  Tsminx1 = 0
  Tsmaxx1 = 0

  !solve using two different numerical schemes
  A(:)= 0
  B(:)= 0
  C(:)= 1
  D(:)=lambda(:)
  E(:)= 0
  TsEuler(1:lx1, i) = backEuler1D(TsEuler(1:lx1, i-1),A,B,C,D,E,Tsminx1,Tsmaxx1,dt,dx1,dx1i, flagdirichBottom, flagdirichTop)
  TsBDF2(1:lx1, i) = TRBDF21D(TsBDF2(1:lx1, i-1),A,B,C,D,E,Tsminx1,Tsmaxx1,dt,dx1,dx1i, flagdirichBottom, flagdirichTop)

  if(flagdirichBottom /= 1 .or. flagdirichTop /= 1) cycle
  !! these checks only for Dirichlet BCS for now.

  !! first, last x1 is miniscule value, meaningless to compare
  !! check the validity of the numerical solutions at this time step
  !! absolute error comparison due to smaller values
  errorEuler = (TsEuler(2:lx1-1, i) - Tstrue(2:lx1-1, i))
  errorBDF2  = (TsBDF2(2:lx1-1, i) -  Tstrue(2:lx1-1, i))

  if (mod(i,5) == 0) then
    print '(I3,3EN12.3)',i,dt,maxval(abs(errorEuler)),maxval(abs(errorBDF2))
  end if
  if (maxval(abs(errorEuler)) > abs_tol .or. maxval(abs(errorBDF2)) > abs_tol) then
    failed = .true.
    write(stderr, '(I3,3EN12.3,A)') i,dt,maxval(abs(errorEuler)),maxval(abs(errorBDF2)),' <-- error'
  end if
end do

!> output
call h5f%write('/Euler/Ts', TsEuler(1:lx1, 1:))
call h5f%write('/BDF2/Ts', TsBDF2(1:lx1, 1:))
call h5f%write('/True/Ts', Tstrue(1:lx1, 1:))

call h5f%finalize()

if (failed) error stop 'Excessive error, check time step perhaps'

if(flagdirichBottom == 1 .and. flagdirichTop == 1) print *, "OK: diffusion1D:dirichlet"

end program
