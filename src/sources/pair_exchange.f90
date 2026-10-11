module pair_exchange
! Constant-rate, closed two-reservoir model. SI: rho [kg/m^3], Cv [J/m^3/K],
! drag [kg/m^3/s], conductance [J/m^3/s/K]. Rates are explicitly supplied:
! this operator enforces laws, but does not validate their physical magnitudes.
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use phys_consts, only: wp,kB
implicit none(type,external)
private
public :: shared_pair_coefficients, pair_exchange_rates, advance_pair_exchange
contains
pure subroutine shared_pair_coefficients(ma,mb,na,nb,nuab,nuba,phi,psi,drag,conductance)
  real(wp), intent(in) :: ma,mb,na,nb,nuab,nuba,phi,psi
  real(wp), intent(out) :: drag,conductance
  real(wp) :: reciprocal_rate
  if(any(.not.ieee_is_finite([ma,mb,na,nb,nuab,nuba,phi,psi]))) error stop 'nonfinite pair coefficient'
  if(min(ma,mb)<=0.or.min(na,nb,nuab,nuba,phi,psi)<0) error stop 'invalid pair coefficient'
  ! Explicit model proposal: average the two rounded directional momentum
  ! rates once, then use this identical rate in both species equations.
  reciprocal_rate=(ma*na*nuab+mb*nb*nuba)/2
  drag=reciprocal_rate*phi
  ! Three translational degrees of freedom determine the exchanged POWER.
  ! Species-specific heat capacities determine temperature response separately.
  conductance=3*kB*reciprocal_rate*psi/(ma+mb)
end subroutine shared_pair_coefficients

pure subroutine pair_exchange_rates(ma,mb,drag,conductance,Ta,Tb,va,vb,Fa,Fb,Qa,Qb)
  real(wp), intent(in) :: ma,mb,drag,conductance,Ta,Tb,va(3),vb(3)
  real(wp), intent(out) :: Fa(3),Fb(3),Qa,Qb
  real(wp) :: heating,transfer,fraction_a
  if(min(ma,mb)<=0.or.min(drag,conductance)<0) error stop 'invalid pair rate'
  fraction_a=mb/(ma+mb)
  Fa=drag*(vb-va);Fb=-Fa
  heating=drag*sum((va-vb)**2)
  transfer=conductance*(Tb-Ta)
  Qa=transfer+fraction_a*heating
  Qb=-transfer+(1-fraction_a)*heating
end subroutine pair_exchange_rates

pure subroutine advance_pair_exchange(rhoa,rhob,Ca,Cb,fraction_a,drag,conductance,dt,va,vb,Ua,Ub)
  real(wp), intent(in) :: rhoa,rhob,Ca,Cb,fraction_a,drag,conductance,dt
  real(wp), intent(inout) :: va(3),vb(3),Ua,Ub
  real(wp) :: reduced_rho, reduced_C, rate_v,rate_T,delta(3),dP(3),kinetic_heat
  real(wp) :: difference0,friction0,integral,thermal_transfer
  real(wp) :: heat_a,heat_b,friction_difference,thermal_fraction
  if(any(.not.ieee_is_finite([rhoa,rhob,Ca,Cb,fraction_a,drag,conductance,dt,va,vb,Ua,Ub]))) &
    error stop 'nonfinite pair state'
  if(min(rhoa,rhob,Ca,Cb)<=0.or.min(drag,conductance,dt,Ua,Ub)<0.or. &
     fraction_a<0.or.fraction_a>1) error stop 'invalid pair state'
  if(dt==0.or.(drag==0.and.conductance==0)) return
  reduced_rho=rhoa*(rhob/(rhoa+rhob))
  reduced_C=Ca*(Cb/(Ca+Cb))
  rate_v=drag/reduced_rho
  rate_T=conductance/reduced_C
  delta=va-vb
  difference0=Ua/Ca-Ub/Cb
  dP=-reduced_rho*delta*decay_fraction(rate_v*dt)
  va=va+dP/rhoa;vb=vb-dP/rhob
  ! Account for the actual finite-step loss of relative bulk kinetic energy,
  ! rather than integrating frozen velocities as if the loss stayed constant.
  kinetic_heat=(reduced_rho*sum(delta**2)/2)*decay_fraction(2*rate_v*dt)
  friction0=drag*sum(delta**2)
  ! Stable convolution of e^(-rate_T*t) with e^(-2*rate_v*t).
  integral=dt*exp(-min(rate_T,2*rate_v)*dt)*phi_negative(abs(rate_T-2*rate_v)*dt)
  friction_difference=friction0*(fraction_a/Ca-(1-fraction_a)/Cb)*integral
  thermal_fraction=decay_fraction(rate_T*dt)
  ! A single opposite transfer retains small changes without reconstructing
  ! each large reservoir from its equilibrium value.
  thermal_transfer=reduced_C*(-difference0*thermal_fraction+friction_difference)
  heat_a=Ca/(Ca+Cb)*kinetic_heat+thermal_transfer
  heat_b=kinetic_heat-heat_a
  Ua=Ua+heat_a;Ub=Ub+heat_b
end subroutine advance_pair_exchange

pure elemental function phi_negative(z) result(phi)
  real(wp), intent(in) :: z
  real(wp) :: phi
  if(z<1e-3_wp) then
    phi=1-z*(1._wp/2-z*(1._wp/6-z*(1._wp/24-z*(1._wp/120-z*(1._wp/720-z/5040)))))
  else
    phi=(1-exp(-z))/z
  endif
end function phi_negative

pure elemental function decay_fraction(z) result(fraction)
  real(wp), intent(in) :: z
  real(wp) :: fraction
  fraction=z*phi_negative(z)
end function decay_fraction
end module pair_exchange
