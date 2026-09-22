program audit_diffusion_singular
use phys_consts, only: wp
use PDEparabolic, only: backEuler1D,TRBDF21D
implicit none
integer, parameter :: n=5
real(wp) :: u(n),a(n),zero(n),one(n),dx(0:n+2),dxi(n),v(n)
character(16) :: mode
call get_command_argument(1,mode)
u=1; zero=0; one=1; dx=1; dxi=1
select case (trim(mode))
case ('euler')
  a=1
  v=backEuler1D(u,a,zero,zero,one,zero,1._wp,1._wp,1._wp,[0,0],dx,dxi)
case ('tr')
  a=4
  v=TRBDF21D(u,a,zero,zero,one,zero,1._wp,1._wp,1._wp,[0,0],dx,dxi)
case ('bdf2')
  a=3
  v=TRBDF21D(u,a,zero,zero,one,zero,1._wp,1._wp,1._wp,[0,0],dx,dxi)
case default
  error stop 'Unknown singular diffusion case'
end select
print *, 'Unexpected success of singular solve: ',v
end program
