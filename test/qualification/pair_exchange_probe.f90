program pair_exchange_probe
use phys_consts, only: wp,lsp,kB,ms,gammas
use collisions, only: coulomb_colln
use pair_exchange, only: shared_pair_coefficients,pair_exchange_rates,advance_pair_exchange
implicit none(type,external)
real(wp) :: ns(-1:3,-1:3,-1:3,lsp),Ts(-1:3,-1:3,-1:3,lsp),v(-1:3,-1:3,-1:3,lsp)
real(wp) :: nu(1,1,1),nurev(1,1,1),phi(1,1,1),psi(1,1,1),phirev(1,1,1),psirev(1,1,1)
real(wp) :: rhoa,rhob,Ca,Cb,fraction,drag,K,dt,rate,Ua,Ub,va(3),vb(3),va0(3),vb0(3)
real(wp) :: Fa(3),Fb(3),Qa,Qb,drag0,K0,energy_residual,momentum_residual,thermal_residual,fa0(3),fb0(3),q0a,q0b
integer :: a,b,drift,level,u,it,nt
character(4096) :: filename
call get_command_argument(1,filename)
if(len_trim(filename)==0) error stop 'pair probe requires CSV path'
open(newunit=u,file=trim(filename),status='replace')
write(u,'(a)') 'a,b,level,dt,rhoa,rhob,Ca,Cb,fraction_a,drag,K,Ta0,Tb0,'// &
 'va0x,va0y,va0z,vb0x,vb0y,vb0z,vax,vay,vaz,vbx,vby,vbz,Ta,Tb,'// &
 'instantaneous_momentum,instantaneous_energy,instantaneous_thermal'
do a=1,lsp-1
 do b=a+1,lsp
  ns=0;Ts=900;v=0
  ns(:,:,:,a)=2e11_wp;ns(:,:,:,b)=8e10_wp
  Ts(:,:,:,a)=600;Ts(:,:,:,b)=1800
  rhoa=2e11_wp*ms(a);rhob=8e10_wp*ms(b)
  Ca=2e11_wp*kB/(gammas(a)-1);Cb=8e10_wp*kB/(gammas(b)-1)
  fraction=ms(b)/(ms(a)+ms(b))
  do drift=0,1
   va0=real(drift,wp)*[1200._wp,300._wp,-600._wp]
   vb0=real(drift,wp)*[-300._wp,-100._wp,200._wp]
   v(:,:,:,a)=va0(1);v(:,:,:,b)=vb0(1)
   call coulomb_colln(a,b,ns,Ts,v,nu,phi,psi)
   call coulomb_colln(b,a,ns,Ts,v,nurev,phirev,psirev)
   if(maxval(abs(phi-phirev))+maxval(abs(psi-psirev))>1e-13_wp) error stop 'pair response not symmetric'
   call shared_pair_coefficients(ms(a),ms(b),2e11_wp,8e10_wp,nu(1,1,1),nurev(1,1,1), &
     phi(1,1,1),psi(1,1,1),drag0,K0)
   do level=0,8
    drag=drag0;K=K0
    if(level==6) K=2*drag*(1/rhoa+1/rhob)*Ca*(Cb/(Ca+Cb))
    if(level==7) K=0
    if(level==8) drag=0
    rate=max(drag*(1/rhoa+1/rhob),K*(1/Ca+1/Cb))
    select case(level)
     case(0);dt=0
     case(1);dt=1e-12_wp/rate
     case(2);dt=.01_wp/rate
     case(3);dt=.2_wp/rate
     case(4);dt=2._wp/rate
     case(5);dt=50._wp/rate
     case default;dt=.3_wp/rate
    end select
    call pair_exchange_rates(ms(a),ms(b),drag,K,600._wp,1800._wp,va0,vb0,Fa,Fb,Qa,Qb)
    call pair_exchange_rates(ms(a),ms(b),drag,K,600._wp,1800._wp,[0._wp,0._wp,0._wp], &
      [0._wp,0._wp,0._wp],fa0,fb0,q0a,q0b)
    momentum_residual=maxval(abs(Fa+Fb))/max(maxval(abs(Fa)),tiny(1._wp))
    energy_residual=abs(Qa+Qb+sum(Fa*va0+Fb*vb0))/ &
      max(abs(Qa)+abs(Qb)+sum(abs(Fa*va0)+abs(Fb*vb0)),tiny(1._wp))
    thermal_residual=abs(q0a+q0b)/max(abs(q0a)+abs(q0b),tiny(1._wp))
    va=va0;vb=vb0;Ua=Ca*600;Ub=Cb*1800
    call advance_pair_exchange(rhoa,rhob,Ca,Cb,fraction,drag,K,dt,va,vb,Ua,Ub)
    if(dt==0) then
      if(any(va/=va0).or.any(vb/=vb0).or.Ua/=Ca*600.or.Ub/=Cb*1800) error stop 'zero-step pair changed state'
    endif
    write(u,'(3(i0,","),*(es25.17,:,","))') a,b,10*drift+level,dt,rhoa,rhob,Ca,Cb,fraction,drag,K, &
      600._wp,1800._wp,va0,vb0,va,vb,Ua/Ca,Ub/Cb,momentum_residual,energy_residual,thermal_residual
   enddo
  enddo
 enddo
enddo
! The original closed thermal-exchange fixture, now using one shared transfer.
Ca=2e11_wp*kB/(5._wp/3-1);Cb=8e10_wp*kB/(5._wp/3-1);K=.2_wp*Ca
do level=0,2
 nt=5*2**level;dt=1._wp/nt
 Ua=Ca*300;Ub=Cb*1800;va=0;vb=0
 do it=1,nt
  call advance_pair_exchange(1._wp,1._wp,Ca,Cb,.5_wp,0._wp,K,dt,va,vb,Ua,Ub)
 enddo
 write(u,'(3(i0,","),*(es25.17,:,","))') 0,0,level,1._wp,1._wp,1._wp,Ca,Cb,.5_wp,0._wp,K, &
   300._wp,1800._wp,[0._wp,0._wp,0._wp],[0._wp,0._wp,0._wp],va,vb,Ua/Ca,Ub/Cb,0._wp,0._wp,0._wp
enddo
close(u)
end program pair_exchange_probe
