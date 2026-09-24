! SPDX-License-Identifier: Apache-2.0
module transport_audit
! Actual numerical face fluxes and physical-cell integrals for split advection.
! Rank-local interface faces cancel when matched across the MPI decomposition.
use phys_consts, only: wp
use meshobj, only: curvmesh
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
implicit none (type,external)
private
public :: audit_transport_init,audit_step,audit_begin,audit_slice,audit_finish,audit_transport_enabled
public :: audit_mass_start,audit_mass_source,audit_mass_cleanup,audit_mass_finish
public :: audit_etd_source,audit_temperature_floor,audit_energy_operator
logical :: mass_active=.false.
real(wp), allocatable :: volume(:,:,:),mass_old(:,:,:,:),cleanup_old(:,:,:,:)
real(wp) :: mass_source(6)=0,mass_flux(6)=0,mass_cleanup(6)=0
integer :: mass_unit,source_unit,temperature_unit,energy_unit
logical :: audit_transport_enabled=.false.
integer :: unit,quantity,species,axis
real(wp) :: time=0,step=0,terms(7)=0
contains
subroutine audit_transport_init(outdir,rank)
  character(*), intent(in) :: outdir
  integer, intent(in) :: rank
  character(16) :: flag,suffix
  character(:), allocatable :: path
  logical :: exists
  integer :: status
  call get_environment_variable('GEMINI_NUMERICAL_AUDIT',flag,status=status)
  audit_transport_enabled=status==0.and.trim(flag)=='1'
  if(.not.audit_transport_enabled) return
  write(suffix,'(I8.8)') rank
  path=outdir//'/transport-r'//trim(suffix)//'.csv'
  inquire(file=path,exist=exists)
  open(newunit=unit,file=path,status='unknown',position='append',action='write')
  if(.not.exists) write(unit,'(A)') &
    't,dt,quantity,species,axis,before,after,delta,outward,min_face,max_face,scale,residual,normalized_residual'
  path=outdir//'/continuity-r'//trim(suffix)//'.csv'
  inquire(file=path,exist=exists)
  open(newunit=mass_unit,file=path,status='unknown',position='append',action='write')
  if(.not.exists) write(mass_unit,'(A)') &
    't,dt,species,before,after,source,outward,cleanup,residual,scale,normalized_residual,charge_fraction'
  path=outdir//'/sources-r'//trim(suffix)//'.csv'
  inquire(file=path,exist=exists)
  open(newunit=source_unit,file=path,status='unknown',position='append',action='write')
  if(.not.exists) write(source_unit,'(A)') &
    't,dt,quantity,species,before,after,delta,production,external_production,integrated_loss,residual,scale'
  path=outdir//'/temperature-floor-r'//trim(suffix)//'.csv'
  inquire(file=path,exist=exists)
  open(newunit=temperature_unit,file=path,status='unknown',position='append',action='write')
  if(.not.exists) write(temperature_unit,'(A)') 't,dt,stage,species,energy_added_J'
  path=outdir//'/energy-operators-r'//trim(suffix)//'.csv'
  inquire(file=path,exist=exists)
  open(newunit=energy_unit,file=path,status='unknown',position='append',action='write')
  if(.not.exists) write(energy_unit,'(A)') &
    't,dt,stage,species,before_J,after_J,delta_J,term1_J,term2_J,term3_J,term4_J,term5_J,scale_J,residual_J'
end subroutine
subroutine audit_energy_operator(stage,s,before,after,increments,capacity)
  ! Stage 1: pressure compression, artificial-viscous work, zero, zero, zero.
  ! Stage 2: linear reaction, thermal drift, conductive divergence, explicit
  ! heating, imposed endpoint reservoir change. No interior residual source.
  integer, intent(in) :: stage,s
  real(wp), intent(in) :: before(:,:,:),after(:,:,:),increments(:,:,:,:),capacity(:,:,:)
  real(wp) :: b,a,d,contribution(5),scale,r
  integer :: n
  if(.not.mass_active) return
  if(size(increments,4)/=5) error stop 'Invalid energy operator channel count'
  b=sum(before*capacity*volume);a=sum(after*capacity*volume)
  d=sum((after-before)*capacity*volume)
  do n=1,5
    contribution(n)=sum(increments(:,:,:,n)*capacity*volume)
  enddo
  scale=max(sum(abs(before)*capacity*volume),sum(abs(after)*capacity*volume), &
            sum(abs(contribution)),tiny(1._wp))
  r=d-sum(contribution)
  if(.not.all(ieee_is_finite([b,a,d,contribution,scale,r]))) error stop 'Nonfinite energy operator ledger'
  write(energy_unit,'(2(ES24.16E3,","),2(I0,","),9(ES24.16E3,","),ES24.16E3)') &
    time,step,stage,s,b,a,d,contribution,scale,r
  flush(energy_unit)
end subroutine
subroutine audit_etd_source(q,s,before,after,production,loss,dt,external_production)
  ! Independently integrate P-L*y for frozen coefficients. No term is defined
  ! from the observed state difference. Channels are aggregate rates, not an
  ! independent validation of individual reactions or collision coefficients.
  integer, intent(in) :: q,s
  real(wp), intent(in) :: before(:,:,:),after(:,:,:),production(:,:,:),loss(:,:,:),dt
  real(wp), intent(in), optional :: external_production(:,:,:)
  real(wp) :: z,phi1,phi2,term1,term2,integrated_loss,b,a,d,p,ext,scale,r,vals(8)
  integer :: i,j,k,n
  if(.not.mass_active) return
  integrated_loss=0
  do k=1,size(volume,3)
    do j=1,size(volume,2)
      do i=1,size(volume,1)
        if(volume(i,j,k)==0) cycle
        z=loss(i,j,k)*dt
        if(abs(z)<0.1_wp) then
          phi1=1;phi2=0.5_wp;term1=1;term2=0.5_wp
          do n=1,12
            term1=-term1*z/real(n+1,wp);phi1=phi1+term1
            term2=-term2*z/real(n+2,wp);phi2=phi2+term2
          enddo
        else
          phi1=(1-exp(-z))/z;phi2=(1-phi1)/z
        endif
        integrated_loss=integrated_loss+volume(i,j,k)* &
          (before(i,j,k)*z*phi1+production(i,j,k)*dt*z*phi2)
      enddo
    enddo
  enddo
  b=sum(before*volume);a=sum(after*volume);d=sum((after-before)*volume)
  p=sum(production*volume)*dt;ext=0
  if(present(external_production)) ext=sum(external_production*volume)*dt
  r=d-p+integrated_loss
  scale=max(sum(abs(before)*volume),sum(abs(after)*volume),abs(p),abs(integrated_loss),tiny(1._wp))
  vals=[b,a,d,p,ext,integrated_loss,r,scale]
  if(.not.all(ieee_is_finite(vals))) error stop 'Nonfinite source ledger'
  write(source_unit,'(2(ES24.16E3,","),2(I0,","),7(ES24.16E3,","),ES24.16E3)') time,dt,q,s,vals
  flush(source_unit)
end subroutine
subroutine audit_temperature_floor(stage,s,before,after,heat_capacity)
  ! Physical-cell internal energy added by an explicitly applied temperature floor.
  integer, intent(in) :: stage,s
  real(wp), intent(in) :: before(:,:,:),after(:,:,:),heat_capacity(:,:,:)
  real(wp) :: addition
  if(.not.mass_active) return
  addition=sum((after-before)*heat_capacity*volume)
  if(.not.ieee_is_finite(addition)) error stop 'Nonfinite temperature floor ledger'
  write(temperature_unit,'(2(ES24.16E3,","),2(I0,","),ES24.16E3)') time,step,stage,s,addition
  flush(temperature_unit)
end subroutine
subroutine audit_step(t,dt)
  real(wp), intent(in) :: t,dt
  time=t;step=dt
end subroutine
subroutine audit_begin(q,s,a)
  integer, intent(in) :: q,s,a
  quantity=q;species=s;axis=a;terms=0
end subroutine
subroutine audit_slice(before,after,weights,flux,physical,transverse)
  real(wp), intent(in) :: before(:),after(:),weights(:),flux(:),transverse
  logical, intent(in) :: physical(:)
  real(wp) :: b,a,d,f,scale,lo,hi
  integer :: n
  n=size(before)
  if(size(flux)/=n+1) error stop 'Audit face count differs from cell count'
  b=sum(before*weights,mask=physical)*transverse
  a=sum(after*weights,mask=physical)*transverse
  d=sum((after-before)*weights,mask=physical)*transverse
  f=sum(flux(2:n+1)-flux(1:n),mask=physical)*transverse
  lo=0;hi=0
  if(physical(1)) lo=-flux(1)*transverse
  if(physical(n)) hi=flux(n+1)*transverse
  scale=max(sum(abs(before)*weights,mask=physical),sum(abs(after)*weights,mask=physical))*transverse
  terms=terms+[b,a,d,f,lo,hi,scale]
end subroutine
subroutine audit_finish()
  real(wp) :: residual,relative
  if(.not.audit_transport_enabled) return
  if(.not.all(ieee_is_finite(terms))) error stop 'Nonfinite transport budget'
  if(mass_active.and.quantity==1.and.species<=6) mass_flux(species)=mass_flux(species)+terms(4)
  residual=terms(3)+terms(4)
  relative=abs(residual)/max(terms(7),abs(terms(4)),tiny(1._wp))
  write(unit,'(2(ES24.16E3,","),3(I0,","),8(ES24.16E3,","),ES24.16E3)') &
    time,step,quantity,species,axis,terms,residual,relative
  flush(unit)
end subroutine
subroutine audit_mass_start(ns,x)
  real(wp), intent(in) :: ns(-1:,-1:,-1:,:)
  class(curvmesh), intent(in) :: x
  integer :: i,j,k
  if(.not.audit_transport_enabled) return
  if(.not.allocated(volume)) then
    allocate(volume(x%lx1,x%lx2,x%lx3),mass_old(x%lx1,x%lx2,x%lx3,6),cleanup_old(x%lx1,x%lx2,x%lx3,6))
    do k=1,x%lx3
      do j=1,x%lx2
        do i=1,x%lx1
          volume(i,j,k)=0
          if(.not.x%nullpts(i,j,k)) &
            volume(i,j,k)=x%h1(i,j,k)*x%h2(i,j,k)*x%h3(i,j,k)*x%dx1i(i)*x%dx2i(j)*x%dx3i(k)
        enddo
      enddo
    enddo
  endif
  mass_old=ns(1:x%lx1,1:x%lx2,1:x%lx3,1:6)
  mass_source=0;mass_flux=0;mass_cleanup=0;mass_active=.true.
end subroutine
subroutine audit_mass_source(species,before,after)
  integer, intent(in) :: species
  real(wp), intent(in) :: before(:,:,:),after(:,:,:)
  if(.not.mass_active) return
  if(species<1.or.species>6) error stop 'Invalid ion continuity diagnostic species'
  mass_source(species)=mass_source(species)+sum((after-before)*volume)
end subroutine
subroutine audit_mass_cleanup(ns,before)
  real(wp), intent(in) :: ns(-1:,-1:,-1:,:)
  logical, intent(in) :: before
  integer :: i,n1,n2,n3
  if(.not.mass_active) return
  n1=size(volume,1);n2=size(volume,2);n3=size(volume,3)
  if(before) then
    cleanup_old=ns(1:n1,1:n2,1:n3,1:6)
  else
    do i=1,6
      mass_cleanup(i)=mass_cleanup(i)+sum((ns(1:n1,1:n2,1:n3,i)-cleanup_old(:,:,:,i))*volume)
    enddo
  endif
end subroutine
subroutine audit_mass_finish(ns)
  real(wp), intent(in) :: ns(-1:,-1:,-1:,:)
  integer :: i,n1,n2,n3
  real(wp) :: b,a,delta,residual,scale,relative,charge_fraction,values(9)
  if(.not.mass_active) return
  n1=size(volume,1);n2=size(volume,2);n3=size(volume,3)
  charge_fraction=sum(abs(ns(1:n1,1:n2,1:n3,7)-sum(ns(1:n1,1:n2,1:n3,1:6),4))*volume) / &
    max(sum(abs(ns(1:n1,1:n2,1:n3,7))*volume),tiny(1._wp))
  do i=1,6
    b=sum(mass_old(:,:,:,i)*volume);a=sum(ns(1:n1,1:n2,1:n3,i)*volume)
    delta=sum((ns(1:n1,1:n2,1:n3,i)-mass_old(:,:,:,i))*volume)
    residual=delta-mass_source(i)+mass_flux(i)-mass_cleanup(i)
    scale=max(sum(abs(mass_old(:,:,:,i))*volume),sum(abs(ns(1:n1,1:n2,1:n3,i))*volume), &
              abs(mass_source(i)),abs(mass_flux(i)),abs(mass_cleanup(i)),tiny(1._wp))
    relative=abs(residual)/scale
    values=[b,a,mass_source(i),mass_flux(i),mass_cleanup(i),residual,scale,relative,charge_fraction]
    if(.not.all(ieee_is_finite(values))) error stop 'Nonfinite continuity ledger'
    write(mass_unit,'(2(ES24.16E3,","),I0,",",8(ES24.16E3,","),ES24.16E3)') time,step,i,values
  enddo
  flush(mass_unit)
  mass_active=.false.
end subroutine
end module
