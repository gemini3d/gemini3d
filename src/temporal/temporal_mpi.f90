module temporal_mpi

!DO NOT FIX THESE WARNINGS - THEY ARE FOR UNUSED VARIABLES THAT MAY BE LEVERAGED
!IN LATER RELEASES
!/home/zettergm/zettergmdata/GEMINI/temporal/temporal.f90:65:0: warning: unused
!parameter ‘ns’ [-Wunused-parameter]
!   pure subroutine
!dt_calc(tcfl,ns,Ts,vs1,vs2,vs3,B1,B2,B3,dx1i,dx2i,dx3i,potsolve,cour1,cour2,cour3,dt) ^
!/home/zettergm/zettergmdata/GEMINI/temporal/temporal.f90:65:0: warning: unused parameter ‘b1’ [-Wunused-parameter]
!/home/zettergm/zettergmdata/GEMINI/temporal/temporal.f90:65:0: warning: unused parameter ‘b2’ [-Wunused-parameter]
!/home/zettergm/zettergmdata/GEMINI/temporal/temporal.f90:65:0: warning: unused parameter ‘b3’ [-Wunused-parameter]
!/home/zettergm/zettergmdata/GEMINI/temporal/temporal.f90:65:0: warning: unused parameter ‘potsolve’ [-Wunused-parameter]

use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use, intrinsic :: iso_fortran_env, only: int64
use temporal, only: cflcalc
use phys_consts, only: kB,mu0,ms,lsp,pi, wp, debug
use mpimod, only: mpi_realprec, tag=>gemini_mpi, mpi_cfg
use meshobj, only:  curvmesh
use gemini3d_config, only: gemini_cfg

use mpi_f08, only: MPI_COMM_WORLD,MPI_STATUS_IGNORE, mpi_send, mpi_recv, MPI_Allreduce, MPI_MAX

implicit none (type, external)

private
public :: dt_comm, enforce_post_update_cfl

contains
  pure real(wp) function to_next_knot(t,cadence) result(dt)
    real(wp), intent(in) :: t,cadence
    if (.not.all(ieee_is_finite([t,cadence]))) error stop "driver knot: nonfinite time"
    if (cadence<=0) error stop "driver knot: cadence must be positive"
    dt=real(floor(t/cadence,kind=int64)+1_int64,wp)*cadence-t
    if (dt<=epsilon(t)*max(1._wp,abs(t))*8) dt=cadence
  end function to_next_knot

  !> Fail before fluid advancement if newly solved drifts invalidate the selected step.
  !> Retrying requires restoring all stateful driver/potential updates and is not implicit.
  subroutine enforce_post_update_cfl(Ts,vs1,vs2,vs3,x,dt)
    real(wp), intent(in) :: Ts(-1:,-1:,-1:,:),vs1(-1:,-1:,-1:,:), &
                           vs2(-1:,-1:,-1:,:),vs3(-1:,-1:,-1:,:)
    class(curvmesh), intent(in) :: x
    real(wp), intent(in) :: dt
    real(wp) :: local_cfl,global_cfl
    call cflcalc(Ts,vs1,vs2,vs3,x%dl1i,x%dl2i,x%dl3i,dt,local_cfl)
    call MPI_Allreduce(local_cfl,global_cfl,1,mpi_realprec,MPI_MAX,MPI_COMM_WORLD)
    if (global_cfl>1._wp+64*epsilon(global_cfl)) then
      if (mpi_cfg%myid==0) print *, 'Post-update CFL exceeds one: ',global_cfl, ' dt=',dt
      error stop "Updated forcing invalidates the time step; reduce tcfl or refine driver cadence."
    endif
  end subroutine enforce_post_update_cfl

  subroutine dt_comm(t,tout,tglowout,cfg,ns,Ts,vs1,vs2,vs3,B1,B2,B3,x,dt)
    real(wp), intent(in) :: t,tout,tglowout
    type(gemini_cfg), intent(in) :: cfg
    real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: ns,Ts,vs1,vs2,vs3
    real(wp),  dimension(-1:,-1:,-1:), intent(in) :: B1,B2,B3
    class(curvmesh), intent(in) :: x
    real(wp), intent(out) :: dt
    real(wp), dimension(lsp) :: cour1,cour2,cour3
    integer :: iid,isp
    real(wp) :: dttmp

    call dt_calc(cfg%tcfl,ns,Ts,vs1,vs2,vs3,B1,B2,B3,x%dl1i,x%dl2i,x%dl3i,cfg%potsolve,cour1,cour2,cour3,dt)

    if (mpi_cfg%myid/=0) then
      call mpi_send(dt,1,mpi_realprec,0,tag%dt,MPI_COMM_WORLD)
      !! send what I think dt should be
      call mpi_recv(dt,1,mpi_realprec,0,tag%dt,MPI_COMM_WORLD,MPI_STATUS_IGNORE)
      !! receive roots decision
    else
      !> FIGURE OUT GLOBAL DT REQUIRED FOR STABILITY
      do iid=1,mpi_cfg%lid-1
        call mpi_recv(dttmp,1,mpi_realprec,iid,tag%dt,MPI_COMM_WORLD,MPI_STATUS_IGNORE)

        if (dttmp < dt) dt=dttmp
      end do

      !CHECK WHETHER WE'D OVERSTEP OUR TARGET OUTPUT TIME
      !GLOW OUTPUT HAS PRIORITY SINCE IT WILL OUTPUT MORE OFTEN
      if ((cfg%flagglow/=0).and.(tglowout>t).and.(t+dt>tglowout)) then
        dt=tglowout-t
        print *, 'GLOW is throttling dt...'
      end if

      if (tout>t .and. t+dt>tout) then
        dt=tout-t
        if (debug) print *, 'Slowing down for an output...'
      end if

      ! The existing frontend initializes output deadlines to t.  Its first
      ! electrodynamic/fluid update populates derived output quantities before
      ! writing that frame.  Preserve this bootstrap with a CAP, never a floor:
      ! a smaller stability-limited step must not be increased to 1 microsecond.
      if (tout==t .or. (cfg%flagglow/=0 .and. tglowout==t)) dt=min(dt,1e-6_wp)

      ! Do not cross the final time or a supplied driver knot.
      dt=min(dt,cfg%tdur-t)
      if (cfg%flagprecfile/=0) dt=min(dt,to_next_knot(t,cfg%dtprec))
      if (cfg%flagE0file/=0) dt=min(dt,to_next_knot(t,cfg%dtE0))
      if (cfg%flagsolfluxfile/=0) dt=min(dt,to_next_knot(t,cfg%dtsolflux))
      if (cfg%flagneutralBGfile/=0) dt=min(dt,to_next_knot(t,cfg%dtneuBGfile))
      if (cfg%flagdneu/=0) dt=min(dt,to_next_knot(t,cfg%dtneu))
      ! Never increase a stability-limited step to an arbitrary minimum.
      if (.not.ieee_is_finite(dt)) error stop "dt_comm: nonfinite time step"
      if (dt<=0 .or. t+dt<=t) error stop "dt_comm: time step cannot advance the simulation"

      !! SEND GLOBAL DT TO ALL WORKERS
      do iid=1,mpi_cfg%lid-1
        call mpi_send(dt,1,mpi_realprec,iid,tag%dt,MPI_COMM_WORLD)
      end do

      if (debug) then
        print *, 'dt figured to be:  ',dt
        print *, 'x1,x2,x3 courant numbers (root process only!):  '
        do isp=1,lsp
          print '(a,f4.2,a,f4.2,a,f4.2)', '    ',cour1(isp),', ',cour2(isp),', ',cour3(isp)
          !! these are roots courant numbers
        end do
        print *, 'Min and max density:  ',minval(pack(ns(:,:,:,7),.true.)),maxval(pack(ns(:,:,:,7),.true.))
      endif
    end if
  end subroutine dt_comm


  pure subroutine dt_calc(tcfl,ns,Ts,vs1,vs2,vs3,B1,B2,B3,dx1i,dx2i,dx3i,potsolve,cour1,cour2,cour3,dt)
    !------------------------------------------------------------
    !-------COMPUTE TIME STEP SUCH THAT STABILITY CONDITION IS
    !-------SATISFIED.  NOTE THAT THE DIFFERENTIALS ARE ASSUMED
    !-------TO HAVE UNITS OF DISTANCE, SO THEY MUST IMPLICITLY
    !-------INCLUDE THE METRIC FACTORS.
    !------------------------------------------------------------
    real(wp), intent(in) :: tcfl
    real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: ns,Ts,vs1,vs2,vs3
    real(wp),  dimension(-1:,-1:,-1:), intent(in) :: B1,B2,B3
    real(wp), dimension(:,:,:), intent(in) :: dx1i
    real(wp), dimension(:,:,:), intent(in) :: dx2i
    real(wp), dimension(:,:,:), intent(in) :: dx3i
    integer, intent(in) :: potsolve
    real(wp), dimension(lsp), intent(out) :: cour1,cour2,cour3
    real(wp), intent(out) :: dt
    real(wp), dimension(lsp) :: gridrate1,gridrate2,gridrate3
    real(wp) :: vsnd
    real(wp) :: rhom,Bmag,vA
    integer :: lx1,lx2,lx3,ix1,ix2,ix3,isp

    lx1=size(Ts,1)-4
    lx2=size(Ts,2)-4
    lx3=size(Ts,3)-4

    gridrate1=0._wp
    gridrate2=0._wp
    gridrate3=0._wp

    !EVALUATE TIME STEP AGAINST LOCAL SOUND SPEED AND ADVECTION
    do isp=1,lsp
      do ix3=1,lx3
        do ix2=1,lx2
          do ix1=1,lx1
            if (isp<lsp) then
              vsnd=sqrt(kB*Ts(ix1,ix2,ix3,isp)/ms(isp)+5._wp/3._wp*kB*Ts(ix1,ix2,ix3,lsp)/ms(isp))
            else
              vsnd=0._wp
            end if

            gridrate1(isp)=max((vsnd+abs(vs1(ix1,ix2,ix3,isp)))/dx1i(ix1,ix2,ix3),gridrate1(isp))
            gridrate2(isp)=max(abs(vs2(ix1,ix2,ix3,isp))/dx2i(ix1,ix2,ix3),gridrate2(isp))
            gridrate3(isp)=max(abs(vs3(ix1,ix2,ix3,isp))/dx3i(ix1,ix2,ix3),gridrate3(isp))
          end do
        end do
      end do
    end do

    !    !CHECK GRIDRATE MAX AGAINST LOCAL ALFVEN SPEED (IF THIS SIMULATION IS INDUCTIVE)
    !    !NOTE THAT THIS SHOULD REALLY INCLUDE MAGNETOSONIC (FAST) MODES IN GRIDRATE DETERMINATION
    !    !AS OF 2/10/2016 THIS CODE IS NOT USED AT ALL, BUT IS KEPT FOR FUTURE DEVELOPMENT
    !    if (potsolve == 2) then
    !      do ix3=1,lx3
    !        do ix2=1,lx2
    !          do ix1=1,lx1
    !            rhom=0._wp
    !            do isp=1,lsp
    !              rhom=rhom+ms(isp)*ns(ix1,ix2,ix3,isp)
    !            end do
    !            Bmag=sqrt(B1(ix1,ix2,ix3)**2+B2(ix1,ix2,ix3)**2+B3(ix1,ix2,ix3)**2)
    !
    !            vA=Bmag/sqrt(rhom*mu0)
    !
    !            do isp=1,lsp
    !              gridrate1(isp)=max(vA/dx1i(ix1),gridrate1(isp))
    !              gridrate2(isp)=max(vA/dx2i(ix2),gridrate2(isp))
    !              gridrate3(isp)=max(vA/dx3i(ix3),gridrate3(isp))
    !            end do
    !          end do
    !        end do
    !      end do
    !    end if

    !ENFORCE A MINIMUM VALUE FOR THE GRIDRATE (WHICH IS CFL/DT, GRID POINTS PER SECOND)
    gridrate1=max(gridrate1, 1e-10_wp)
    gridrate2=max(gridrate2, 1e-10_wp)
    gridrate3=max(gridrate3, 1e-10_wp)

    dt=tcfl*min(minval(1._wp/gridrate1),minval(1._wp/gridrate2),minval(1._wp/gridrate3))
    cour1=dt*gridrate1
    cour2=dt*gridrate2
    cour3=dt*gridrate3
  end subroutine dt_calc
end module temporal_mpi
