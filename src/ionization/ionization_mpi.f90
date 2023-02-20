module ionization_mpi

use phys_consts, only: wp,debug
use grid, only: lx1,lx2,lx3
use mpimod, only: mpi_realprec, mpi_cfg, tag=>gemini_mpi
use neutral, only: neutral_info

use mpi_f08, only : MPI_COMM_WORLD,MPI_STATUS_IGNORE, mpi_send,mpi_recv

implicit none (type, external)

private
public :: get_gavg_Tinf

contains
  !> Query workers to get a single value for Tninf and gavg (or makes something up)
  subroutine get_gavg_Tinf(atmos,gavg,Tninf)
    type(neutral_info), intent(in) :: atmos
    real(wp), intent(out) :: gavg,Tninf
    real(wp) :: Tninftmp
    integer :: iid

    ! use an average value for the gravitational field; FIXME: perhaps should be done via averaging over all workers???
    !allocate(g(1:lx1,1:lx2,1:lx3))
    !g=sqrt(g1**2+g2**2+g3**2)
    !    gavg=sum(g)/(lx1*lx2*lx3)    !single average value for computing column dens.  Interestingly this is a worker average...  Do we need root grav vars. grid mod to prevent tearing?  Should be okay as long as the grid is only sliced along the x3-dimension; problematic for x2divisions...
    gavg=8._wp

    Tninf=maxval(atmos%Tnmsis)   !set exospheric temperature based on the max value of the background MSIS atmosphere; note this is a worker max

    !both g and Tinf need to be computed as average over the entire grid...
    if (mpi_cfg%myid==0) then     !root
      do iid=1,mpi_cfg%lid-1
          call mpi_recv(Tninftmp,1,mpi_realprec,iid,tag%Tninf,MPI_COMM_WORLD,MPI_STATUS_IGNORE)
          if (Tninf < Tninftmp) Tninf=Tninftmp
      end do

      do iid=1,mpi_cfg%lid-1
        call mpi_send(Tninf,1,mpi_realprec,iid,tag%Tninf,MPI_COMM_WORLD)
      end do

      if (debug) print *, 'Exospheric temperature used for photoionization:  ',Tninf
    else                  !workders
      call mpi_send(Tninf,1,mpi_realprec,0,tag%Tninf,MPI_COMM_WORLD)
      !! send what I think Tninf should be

      call mpi_recv(Tninf,1,mpi_realprec,0,tag%Tninf,MPI_COMM_WORLD,MPI_STATUS_IGNORE)
      !! receive roots decision
    end if
    !deallocate(g)
  end subroutine
end module ionization_mpi
