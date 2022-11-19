!> This modules contains subroutines for input and output meant to be used from applications which are
!    mpi-agnostic
module io_nompi

use phys_consts, only : lsp,wp,comp_lvl,mindens
use interpolation, only : interp3
use grid, only : gridflag,lx1,lx2,lx3,get_grid3_coords_hdf5
use timeutils, only : date_filename
use h5fortran, only: hdf5_file
use reader, only : get_simsize3

private
public :: interp_file2subgrid, plasma_output_nompi

contains
  !> Interpolate initial conditions onto "local" subgrid; we assume that the input data grid is specified
  !    by the input file, whereas the target grid *could* be different, e.g. due to refinement or some other
  !    custom arrangement.  The entire input file will be read by each worker calling this procedure.
  !  Since this is only performing spatial interpolation it is easiest to just use the interpolation module
  !    directly rather than create a type extension for inputdata (which inherently wants to also do time interpolation)
  !    and then overriding the interp to space-only.
  subroutine interp_file2subgrid(indatsize,indatfile,indatgrid,x1,x2,x3,ns,vs1,Ts,Phi)
    character(*), intent(in) :: indatsize,indatfile,indatgrid
    real(wp), dimension(-1:) :: x1
    real(wp), dimension(-1:) :: x2
    real(wp), dimension(-1:) :: x3
    real(wp), dimension(-1:,-1:,-1:,:), intent(inout) :: ns,vs1,Ts
    real(wp), dimension(-1:,-1:,-1:), intent(inout) :: Phi
    real(wp), dimension(:,:,:,:), allocatable :: nsall,vs1all,Tsall
    real(wp), dimension(:,:,:), allocatable :: Phiall
    real(wp), dimension(:), allocatable :: parmflat
    integer :: lx1,lx2,lx3
    integer :: lx1in,lx2in,lx3in
    real(wp), dimension(:), allocatable :: x1in,x2in,x3in
    real(wp) :: glatctr=0._wp, glonctr=0._wp
    integer :: isp
    real(wp), dimension(:,:,:), allocatable :: x1imat,x2imat,x3imat    ! variables for interpolation sites
    real(wp), dimension(:), allocatable :: x1i,x2i,x3i    ! variables for interpolation sites
    integer :: ix1,ix2,ix3

    ! set initially to zero to cover ghost cell data so that it is initialized with random junk
    ns=0._wp
    vs1=0._wp
    Ts=0._wp
    Phi=0._wp

    ! convenience
    lx1=size(x1)-4; lx2=size(x2)-4; lx3=size(x3)-4;

    ! read in the ICs size and allocate data
    print*, 'indatsize;  ',indatsize
    call get_simsize3(indatsize,lx1in,lx2in,lx3in)
    allocate(x1in(-1:lx1in+2),x2in(-1:lx2in+2),x3in(-1:lx3in+2))
    allocate(nsall(-1:lx1in+2,-1:lx2in+2,-1:lx3in+2,1:lsp), &
              vs1all(-1:lx1in+2,-1:lx2in+2,-1:lx3in+2,1:lsp), &
              Tsall(-1:lx1in+2,-1:lx2in+2,-1:lx3in+2,1:lsp), &
              Phiall(-1:lx1in+2,-1:lx2in+2,-1:lx3in+2))
    allocate(parmflat(lx1*lx2*lx3))

    ! allocate space for the target coordinates
    allocate(x1imat(lx1,lx2,lx3))
    allocate(x2imat,x3imat,mold=x1imat)
    allocate(x1i(lx1*lx2*lx3))
    allocate(x2i,x3i,mold=x1i)
    do ix3=1,lx3
      do ix2=1,lx2
         do ix1=1,lx1
           x1imat(ix1,ix2,ix3)=x1(ix1)
           x2imat(ix1,ix2,ix3)=x2(ix2)
           x3imat(ix1,ix2,ix3)=x3(ix3)
         end do
      end do
    end do
    x1i=pack(x1imat,.true.)
    x2i=pack(x2imat,.true.)
    x3i=pack(x3imat,.true.)
    deallocate(x1imat,x2imat,x3imat)    ! nuke these as soon as we are done with them

    ! get the input grid coordinates
    print*, 'indatgrid:  ',indatgrid
    call get_grid3_coords_hdf5(indatgrid,x1in,x2in,x3in,glonctr,glatctr)

    !print*, x1in(1:lx1in)
    !print*, '====================================================================================='
    !print*, x1(1:lx1)

    ! we must make sure that the target coordinates do not range outside the input file coordinates
    print*, 'check grid extents'
    if(x1(1)<x1in(1) .or. x1(lx1)>x1in(lx1in)) then
      error stop 'interp_file2grid: x1 target coordinates beyond input grid coords'
    end if
    if(x2(1)<x2in(1) .or. x2(lx2)>x2in(lx2in)) then
      error stop 'interp_file2grid: x2 target coordinates beyond input grid coords'
    end if
    if(x3(1)<x3in(1) .or. x3(lx3)>x3in(lx3in)) then
      error stop 'interp_file2grid: x3 target coordinates beyond input grid coords'
    end if

    ! read in the input initial conditions, only hdf5 files are support for this functionality
    print*, 'read file'
    call getICs_hdf5_nompi(indatsize,indatfile,nsall,vs1all,Tsall,Phiall)

    print*, 'interp_file2subgrid:  error checking complete...'

    ! interpolation input data to mesh sites; do not interpolate to ghost cells
    do isp=1,lsp
      parmflat=interp3(x1in(1:lx1in),x2in(1:lx2in),x3in(1:lx3in),nsall(1:lx1in,1:lx2in,1:lx3in,isp), &
                         x1i(1:lx1*lx2*lx3),x2i(1:lx1*lx2*lx3),x3i(1:lx1*lx2*lx3))
      ns(1:lx1,1:lx2,1:lx3,isp)=reshape(parmflat,[lx1,lx2,lx3])
      parmflat=interp3(x1in(1:lx1in),x2in(1:lx2in),x3in(1:lx3in),vs1all(1:lx1in,1:lx2in,1:lx3in,isp), &
                         x1i(1:lx1*lx2*lx3),x2i(1:lx1*lx2*lx3),x3i(1:lx1*lx2*lx3))
      vs1(1:lx1,1:lx2,1:lx3,isp)=reshape(parmflat,[lx1,lx2,lx3])
      parmflat=interp3(x1in(1:lx1in),x2in(1:lx2in),x3in(1:lx3in),Tsall(1:lx1in,1:lx2in,1:lx3in,isp), &
                         x1i(1:lx1*lx2*lx3),x2i(1:lx1*lx2*lx3),x3i(1:lx1*lx2*lx3))
      Ts(1:lx1,1:lx2,1:lx3,isp)=reshape(parmflat,[lx1,lx2,lx3])
    end do
    parmflat=interp3(x1in(1:lx1in),x2in(1:lx2in),x3in(1:lx3in),Phiall(1:lx1in,1:lx2in,1:lx3in), &
                       x1i(1:lx1*lx2*lx3),x2i(1:lx1*lx2*lx3),x3i(1:lx1*lx2*lx3))
    Phi(1:lx1,1:lx2,1:lx3)=reshape(parmflat,[lx1,lx2,lx3])
    print*, 'interp_file2subgrid:  interpolations complete...'

    ! at this point to be totally safe we should set the ghost cells, use a zero-order hold as a total guess
    call forceinputZOH(ns)
    call forceinputZOH(vs1)
    call forceinputZOH(Ts)

    ! get rid of local vars
    deallocate(x1in,x2in,x3in,nsall,vs1all,Tsall,Phiall,parmflat)
    deallocate(x1i,x2i,x3i)
  end subroutine interp_file2subgrid


  subroutine forceinputZOH(param)
    real(wp), dimension(-1:,-1:,-1:,:), intent(inout) :: param
    integer :: lx1,lx2,lx3
  
    lx1=size(param,1)-4; lx2=size(param,2)-4; lx3=size(param,3)-4;
  
    param(0,:,:,:)=param(1,:,:,:)
    param(-1,:,:,:)=param(1,:,:,:)
    param(lx1+1,:,:,:)=param(lx1,:,:,:)
    param(lx1+2,:,:,:)=param(lx1,:,:,:)
  
    param(:,0,:,:)=param(:,1,:,:)
    param(:,-1,:,:)=param(:,1,:,:)
    param(:,lx2+1,:,:)=param(:,lx2,:,:)
    param(:,lx2+2,:,:)=param(:,lx2,:,:)
  
    param(:,:,0,:)=param(:,:,1,:)
    param(:,:,-1,:)=param(:,:,1,:)
    param(:,:,lx3+1,:)=param(:,:,lx3,:)
    param(:,:,lx3+2,:)=param(:,:,lx3,:)
  end subroutine forceinputZOH


  !> output just a the local subgrid data to a file
  subroutine plasma_output_nompi(outdir,flagoutput,ymd,UTsec,ns,vs1,vs2,vs3,Ts, &
                                                Phi,J1,J2,J3,identifier,x1lims,x2lims,x3lims)
    character(*), intent(in) :: outdir
    integer, intent(in) :: flagoutput
    integer, dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec
    real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(-1:,-1:,-1:), intent(in) :: Phi   ! okay to have ghost cells b/c already resides on root.
    real(wp), dimension(1:,1:,1:), intent(in) :: J1,J2,J3
    integer, intent(in), optional :: identifier
    real(wp), dimension(2), intent(in), optional :: x1lims,x2lims,x3lims

    character(:), allocatable :: filenamefull
    character(10) :: IDstr
    integer :: isp
    type(hdf5_file) :: hout
    real(wp), dimension(1:lx1,1:lx2,1:lx3) :: v2avg,v3avg,v1avg,Tavg,ne,Te

    !> Compute some averages for the output file
    v2avg=sum(ns(1:lx1,1:lx2,1:lx3,1:lsp-1)*vs2(1:lx1,1:lx2,1:lx3,1:lsp-1),4)
    v2avg=v2avg/ns(1:lx1,1:lx2,1:lx3,lsp)    !compute averages for output.
    v3avg=sum(ns(1:lx1,1:lx2,1:lx3,1:lsp-1)*vs3(1:lx1,1:lx2,1:lx3,1:lsp-1),4)
    v3avg=v3avg/ns(1:lx1,1:lx2,1:lx3,lsp)
    if (flagoutput==2 .or. flagoutput==3) then
      ne=ns(1:lx1,1:lx2,1:lx3,lsp)
    end if
    if (flagoutput==2) then
      v1avg=sum(ns(1:lx1,1:lx2,1:lx3,1:lsp-1)*vs1(1:lx1,1:lx2,1:lx3,1:lsp-1),4)
      v1avg=v1avg/ns(1:lx1,1:lx2,1:lx3,lsp)    !compute averages for output.
      Tavg=sum(ns(1:lx1,1:lx2,1:lx3,1:lsp-1)*Ts(1:lx1,1:lx2,1:lx3,1:lsp-1),4)
      Tavg=Tavg/ns(1:lx1,1:lx2,1:lx3,lsp)    !compute averages for output.
      Te=Ts(1:lx1,1:lx2,1:lx3,lsp)
    end if

    !> get filename
    if (.not. present(identifier)) then
      filenamefull = date_filename(outdir,ymd,UTsec) // '.h5'
    else
      write(IDstr,'(I0)') identifier
      filenamefull = date_filename(outdir,ymd,UTsec) // '_' // trim(IDstr) // '.h5'
    end if
    print *, 'HDF5 Output file name:  ', filenamefull

    call hout%open(filenamefull, action='w',comp_lvl=comp_lvl)
    call hout%write("/flagoutput", flagoutput)
    call hout%write('/time/ymd', ymd)
    call hout%write('/time/UThour',   real(UTsec/3600.))
    if (present(x1lims)) call hout%write('/x1lims',real(x1lims))
    if (present(x2lims)) call hout%write('/x2lims',real(x2lims))
    if (present(x3lims)) call hout%write('/x3lims',real(x3lims))
    if (present(identifier)) call hout%write('/patchID',identifier)

    select case (flagoutput)
      case (2)    !output ISR-like average parameters
        call hout%write('neall',    real(ne(1:lx1,1:lx2,1:lx3)))
        call hout%write('v1avgall', real(v1avg(1:lx1,1:lx2,1:lx3)))
        !output of ISR-like parameters (ne,Ti,Te,v1,etc.)
        call hout%write('Tavgall',  real(Tavg(1:lx1,1:lx2,1:lx3)))
        call hout%write('Teall',    real(Te(1:lx1,1:lx2,1:lx3)))
        call hout%write('J1all',    real(J1(1:lx1,1:lx2,1:lx3)))
        call hout%write('J2all',    real(J2(1:lx1,1:lx2,1:lx3)))
        call hout%write('J3all',    real(J3(1:lx1,1:lx2,1:lx3)))
        call hout%write('v2avgall', real(v2avg(1:lx1,1:lx2,1:lx3)))
        call hout%write('v3avgall', real(v3avg(1:lx1,1:lx2,1:lx3)))
      case (3)     !just electron density
        print *, 'INFO:  Input file has selected electron density only output, make sure this is what you really want!'
        call hout%write('neall',    real(ne(1:lx1,1:lx2,1:lx3)))
      case default    !output everything
        print *, 'INFO:  Input file has selected full output or milestones, large files may result!'
        call hout%write('nsall',    real(ns(1:lx1,1:lx2,1:lx3,:)))
        call hout%write('vs1all',   real(vs1(1:lx1,1:lx2,1:lx3,:)))
        !this is full output of all parameters in 3D
        call hout%write('Tsall',    real(Ts(1:lx1,1:lx2,1:lx3,:)))

        call hout%write('J1all',    real(J1(1:lx1,1:lx2,1:lx3)))
        call hout%write('J2all',    real(J2(1:lx1,1:lx2,1:lx3)))
        call hout%write('J3all',    real(J3(1:lx1,1:lx2,1:lx3)))
        call hout%write('v2avgall', real(v2avg(1:lx1,1:lx2,1:lx3)))
        call hout%write('v3avgall', real(v3avg(1:lx1,1:lx2,1:lx3)))
    end select

    if (gridflag==1) then
      print *, 'Writing topside boundary conditions for inverted-type grid...'
      call hout%write('Phiall',       real(Phi(1,1:lx2,1:lx3)))
    else
      print *, 'Writing topside boundary conditions for non-inverted-type grid...'
      call hout%write('Phiall',       real(Phi(lx1,1:lx2,1:lx3)))
    end if

    call hout%close()
  end subroutine plasma_output_nompi


  !> This may only differ by variable names from what is in read_hdf, but I needed a copy in a file
  !    that doesn't depend on mpi libs
  subroutine getICs_hdf5_nompi(indatsize,indatfile,ns,vs1,Ts,Phi)
    character(*), intent(in) :: indatsize, indatfile
    real(wp), dimension(-1:,-1:,-1:,:), intent(inout) :: ns,vs1,Ts
    real(wp), dimension(-1:,-1:,-1:), intent(inout) :: Phi
    type(hdf5_file) :: hf
    integer :: lx1,lx2,lx3,isp
    integer :: ix1
    integer :: lx1in,lx2in,lx3in,u,utrace
    real(wp), dimension(:,:), allocatable :: Phislab
    real(wp), allocatable :: tmp(:,:,:,:), tmpPhi(:), tmpPhi2(:,:)

    !integer ix2,ix3
    !real(wp), dimension(:,:,:,:), allocatable :: tmpread

    !> so that random values (including NaN) don't show up in Ghost cells

    !> SYSTEM SIZES
    lx1=size(ns,1)-4
    lx2=size(ns,2)-4
    lx3=size(ns,3)-4

    allocate(Phislab(1:lx2,1:lx3))  !space to store EFL potential

    !> READ IN FROM FILE, AS OF CURVILINEAR BRANCH THIS IS NOW THE ONLY INPUT OPTION
    call get_simsize3(indatsize, lx1in, lx2in, lx3in)
    print '(2A,3I6)', indatsize,' input dimensions:',lx1in,lx2in,lx3in
    print '(A,3I6)', 'Target (output) grid structure dimensions:',lx1,lx2,lx3

    if (.not. (lx1==lx1in .and. lx2==lx2in .and. lx3==lx3in)) then
      error stop 'ERROR:gemini3d: The input data must be the same size as the grid which you are running the simulation on' // &
           '- use a script to interpolate up/down to the simulation grid'
    end if

    print*, 'opening hdf5 file...'
    call hf%open(indatfile, action='r')

!    print*, 'test setting values...'
!    do ix1=-1,lx1+2
!      do ix2=-1,lx2+2
!        do ix3=-1,lx3+2
!          do isp=1,lsp
!            ns(ix1,ix2,ix3,isp)=-1.0
!            vs1(ix1,ix2,ix3,isp)=-2.0
!            Ts(ix1,ix2,ix3,isp)=-3.0
!          end do
!        end do
!      end do
!    end do
!    Phislab=0.0
!    print*, shape(ns)
!    print*, shape(vs1)
!    print*, shape(Ts)
!    print*, maxval(abs(ns))
!    print*, maxval(abs(vs1))
!    print*, maxval(-1.0*abs(Ts))
!    print*, minval(Ts),maxval(Ts)

!    allocate(tmpread(1:lx1,1:lx2,1:lx3,1:lsp))
!    print*, 'tmp reading in fluid data 1'
!    call hf%read('/nsall', tmpread)
!    print*, 'tmp reading in fluid data 2'
!    call hf%read('/vs1all', tmpread)
!    print*, 'tmp reading in fluid data 3'
!    call hf%read('/Tsall', tmpread)
!    deallocate(tmpread)

    !print*, 'reading in fluid data 1'
    call hf%read('/nsall', ns(1:lx1,1:lx2,1:lx3,1:lsp))
    !print*, 'reading in fluid data 2'
    call hf%read('/vs1all', vs1(1:lx1,1:lx2,1:lx3,1:lsp))
    !print*, 'reading in fluid data 3'
    call hf%read('/Tsall', Ts(1:lx1,1:lx2,1:lx3,1:lsp))

    !print*, 'reading in potential...'
    if (hf%exist('/Phiall')) then
      if (hf%ndim('/Phiall') == 1) then
        if (lx2==1) then
          allocate(tmpPhi(lx3))
        else
          allocate(tmpPhi(lx2))
        end if
        call hf%read('/Phiall', tmpPhi)
        ! FIXME: MH please delete if you are okay with this
        !if (size(Phislab, 1) /= 1) then
        !  write(stderr,*) 'Phislab shape',shape(Phislab)
        !  error stop 'Phislab x2 /= 1'
        !endif
        if (lx2==1) then
          Phislab(1,:) = tmpPhi
        else
          Phislab(:,1)=tmpPhi
        end if
      else
        call hf%read('/Phiall', Phislab)
      endif
    else
      Phislab = 0
    end if

    !print*, 'closing hdf5 file...'
    call hf%close()

    !> Apply EFL approx to compute full grid potential
    !print*, 'apply EFL approximation...'
    do ix1=1,lx1
      Phi(ix1,1:lx2,1:lx3)=Phislab(1:lx2,1:lx3)
    end do

    deallocate(Phislab)    ! explicitly get rid of allocated storage
  end subroutine getICs_hdf5_nompi
end module io_nompi
