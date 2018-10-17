program test_diffusion1D

!------------------------------------------------------------
!-------SOLVE HEAT EQUATION IN 1D
!------------------------------------------------------------
use phys_consts, only: wp
use diffusion
implicit none

integer, parameter :: npts=25,lt=100
character(*), parameter :: outfn='test_diffusion1d.dat'

real(wp), dimension(npts) :: v1,dx1i
real(wp), dimension(-1:npts+2) :: x1,Ts
real(wp), dimension(npts) :: lambda,A,B,C,D,E
real(wp), dimension(npts+1) :: x1i
real(wp), dimension(0:npts+2) :: dx1
integer :: lx1,it,ix1,u
real(wp) :: t=0,dt
real(wp) :: Tsminx1,Tsmaxx1


!include ghost cell locations
x1=[ (dble(ix1)/dble(npts), ix1=-2,npts+1) ]
!x1=[ (ix1, ix1=-2,npts+1) ]
lx1=npts   !exclude ghost cells in count

  print *,'writing ',outfn
open(newunit=u,file=outfn,status='replace')
write(u,*) lt
write(u,*) lx1
call writearray(u,x1)


!IC's
Ts(:)=0.0
Ts(11:15)=1.0
lambda(:)=1.0

!backward diffs
dx1=x1(0:lx1+2)-x1(-1:lx1+1)


!interface data
x1i(1:lx1+1)=0.5*(x1(0:lx1)+x1(1:lx1+1))
dx1i=x1i(2:lx1+1)-x1i(1:lx1)


!typical diffusion time
!dt=0.1*maxval(dx1)**2/maxval(lambda)
dt=0.04/lt


do it=1,lt
  !time step
  t=t+dt

  !boundary values
  Tsminx1=0.0
  Tsmaxx1=0.0

  !diffuse
  A(:)=0.0
  B(:)=0.0
  C(:)=1.0
  D(:)=lambda(:)
  E(:)=0.0
  Ts(1:lx1)=backEuler1D(Ts(1:lx1),A,B,C,D,E,Tsminx1,Tsmaxx1,dt,dx1,dx1i)

  !output
  write(u,*) t
  call writearray(u,Ts(1:lx1))
end do

 close(u)


contains

  subroutine writearray(u,array)
    integer, intent(in) :: u
    real(wp), dimension(:), intent(in) :: array
    
    integer :: k

    do k=1,size(array)
      write(u,*) array(k)
    end do
  end subroutine writearray


  subroutine write2Darray(u,array)
    integer, intent(in) :: u
    real(wp), dimension(:,:), intent(in) :: array
    
    integer :: k1,k2

    do k1=1,size(array,1)
      write(u,'(f8.0)') (array(k1,k2), k2=1,size(array,2))
    end do
  end subroutine write2Darray


end program test_diffusion1D
