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
integer, parameter :: Npts=256, Nt=20*5

character(4) :: ic

real(wp), parameter :: abs_tol = 0.05_wp

real(wp), dimension(npts) :: v1,dx1i
real(wp), dimension(-1:npts+2) :: x1
real(wp), dimension(-1:npts+2, 0:Nt) :: TsEuler,TsBDF2,TsTrue
real(wp), dimension(npts) :: lambda,A,B,C,D,E
real(wp), dimension(npts+1) :: x1i
real(wp), dimension(0:npts+2) :: dx1
integer :: lx1,i,ix1, flagdirichTop, flagdirichBottom
real(wp) :: t(Nt), dt
real(wp) :: Tsminx1,Tsmaxx1

real(wp), dimension(npts-2) :: errorEuler,errorBDF2
logical :: failed = .false.
character(256) :: argv
character(:), allocatable :: outfile

flagdirichBottom = 1  !< Dirichlet always for bottom, allow top to be Dirichlet or Neumann

flagdirichTop = 1
call get_command_argument(1, argv, status=i)
if(i==0) read(argv,'(I2)') flagdirichTop

outfile = 'test_diffusion1d.h5'
call get_command_argument(2, argv, status=i)
if(i==0) outfile = trim(argv)

!! create a grid for the calculation
x1=[ (real(ix1-1) / real(npts-1), ix1=-1,npts+2) ]
lx1=npts   !exclude ghost cells in count
dx1=x1(0:lx1+2)-x1(-1:lx1+1)
x1i(1:lx1+1)=0.5*(x1(0:lx1)+x1(1:lx1+1))
dx1i=x1i(2:lx1+1)-x1i(1:lx1)

!! initial conditions
TsEuler(-1:lx1+2, 0) = sin(2*pi*x1(-1:lx1+2)) + sin(8*pi*x1(-1:lx1+2))
TsBDF2(:, 0) = TsEuler(:, 0)
lambda = 1     !< thermal conductivity

!! typical diffusion time, make our time step a fraction of this
dt = 0.05 / 64 / pi**2 / maxval(lambda)
t(1) = 0
do i = 2, Nt
  t(i) = t(i-1) + dt
end do

if (flagdirichTop == 1) then
  print *, "Dirichlet top boundary condition"
else
  print *, "Neumann top boundary condition"
endif
print '(A3, 3A12)', 'i','dt  ', 'Error Euler', 'Error TBDF'
!! time interations

call ideal(t, dt, x1, lx1, flagdirichTop, TsTrue)

do i = 1, Nt
  write(ic, '(I4.4)') i
  !boundary values
  Tsminx1 = 0
  Tsmaxx1 = 0

  !solve using different numerical schemes
  A(:) = 0
  B(:) = 0
  C(:) = 1
  D(:) = lambda
  E(:) = 0

  TsEuler(1:lx1, i) = backEuler1D(TsEuler(1:lx1, i-1), A,B,C,D,E,Tsminx1,Tsmaxx1,dt,dx1,dx1i, flagdirichBottom, flagdirichTop)
  TsBDF2(1:lx1, i) = TRBDF21D(TsBDF2(1:lx1, i-1), A,B,C,D,E,Tsminx1,Tsmaxx1,dt,dx1,dx1i, flagdirichBottom, flagdirichTop)

  !! first, last x1 is miniscule value, meaningless to compare
  !! check the validity of the numerical solutions at this time step
  !! absolute error comparison due to smaller values
  errorEuler = (TsEuler(2:lx1-1, i) - TsTrue(2:lx1-1, i))
  errorBDF2  = (TsBDF2(2:lx1-1, i) -  TsTrue(2:lx1-1, i))

  if (maxval(abs(errorEuler)) > abs_tol .or. maxval(abs(errorBDF2)) > abs_tol) then
    failed = .true.
    write(stderr, '(I3,3EN12.3,A)') i,dt,maxval(abs(errorEuler)),maxval(abs(errorBDF2)),' <-- error'
    cycle
  end if

  if (mod(i,5) == 0) then
    print '(I3,3EN12.3)',i,dt,maxval(abs(errorEuler)),maxval(abs(errorBDF2))
  end if

end do

call write_file(outfile)

if (failed) error stop 'Excessive error, check time step perhaps'

if (flagdirichTop == 1) then
  print *, "OK: diffusion1D:dirichlet"
else
  print *, "OK: diffusion1D:neumann"
endif


contains


subroutine write_file(outfile)
!! write the time, space length and spatial grid to a file

character(*), intent(in) :: outfile

print *,'writing ',outfile
call h5f%initialize(outfile, status='replace', action='write')
call h5f%write('/lx1', lx1)
call h5f%write('/x1', x1)
call h5f%write('/t', t)
call h5f%write('/flagdirichTop', flagdirichTop)
call h5f%write('/flagdirichBottom', flagdirichBottom)
call h5f%write('/Euler/Ts', TsEuler(1:lx1, 1:))
call h5f%write('/BDF2/Ts', TsBDF2(1:lx1, 1:))
call h5f%write('/True/Ts', TsTrue(1:lx1, 1:))
call h5f%finalize()
end subroutine write_file


pure subroutine ideal(t,dt, z,Lz, flagdirichTop, TsTrue)

integer, intent(in) :: lz, flagdirichTop
real(wp), intent(in) :: dt, t(:), z(-1:lz+2)
real(wp), dimension(-1:lz+2, 0:size(t)), intent(out) :: TsTrue

integer :: i
!! compute analytical solution to compare
!! equations from https://github.com/gemini3d/GEMINI-docs/blob/master/test_descriptions/GEMINItests.pdf
do i = 1,size(t)

  if (flagdirichTop==1) then
    TsTrue(1:Lz, i) = exp(-4 *pi**2 * lambda * (t(i) + dt)) * sin(2*pi * z(1:lz)) + &
                            exp(-64*pi**2 * lambda * (t(i) + dt)) * sin(8*pi * z(1:lz))
    !! Equation (17)
  else
    TsTrue(1:Lz, i) = exp(-25*pi**2 / 4 * lambda * (t(i) + dt)) * sin(5*pi / 2 * z(1:lz)) + &
                            exp(-289*pi**2 / 4 * lambda * (t(i) + dt)) * sin(17*pi / 2 * z(1:lz))
    !! Equation (25)
  endif
end do

end subroutine ideal

end program
