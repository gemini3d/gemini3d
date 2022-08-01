submodule (io:plasma_input) plasma_input_hdf5



implicit none (type, external)

contains

  module procedure input_root_mpi_hdf5
    !! READ INPUT FROM FILE AND DISTRIBUTE TO WORKERS.
    !! STATE VARS ARE EXPECTED INCLUDE GHOST CELLS.  NOTE ALSO
    !! THAT RECORD-BASED INPUT IS USED SO NO FILES > 2GB DUE
    !! TO GFORTRAN BUG WHICH DISALLOWS 8 BYTE INTEGER RECORD
    !! LENGTHS.
    real(wp), dimension(-1:size(x1,1)-2,-1:size(x2all,1)-2,-1:size(x3all,1)-2,1:lsp) :: nsall, vs1all, Tsall
    real(wp) :: tstart,tfin

    !> to avoid having garbage in ghost cells
    nsall = 0
    ns = 0
    vs1all= 0
    vs1 = 0
    Tsall = 0
    Ts = 0

    !> read in the full initial conditions files
    call getICs_hdf5(indatsize,indatfile,nsall,vs1all,Tsall,Phiall)

    !> ROOT BROADCASTS IC DATA TO WORKERS
    call cpu_time(tstart)
    call bcast_send(nsall,tag%ns,ns)
    call bcast_send(vs1all,tag%vs1,vs1)
    call bcast_send(Tsall,tag%Ts,Ts)
    !call bcast_send(Phiall,tag%Phi,Phi)
    call bcast_send3D_ghost(Phiall,tag%Phi,Phi)
    call cpu_time(tfin)
    print '(A,ES12.3,A)', 'Sent ICs to workers in', tfin-tstart, ' seconds.'
  end procedure input_root_mpi_hdf5


  !> Read in a full dataset from an input file
  module procedure getICs_hdf5
    type(hdf5_file) :: hf
    integer :: lx1,lx2all,lx3all,isp
    integer :: ix1
    integer :: lx1in,lx2in,lx3in,u,utrace
    real(wp), dimension(:,:), allocatable :: Phislab
    real(wp), allocatable :: tmp(:,:,:,:), tmpPhi(:), tmpPhi2(:,:)

    !> so that random values (including NaN) don't show up in Ghost cells

    !> SYSTEM SIZES
    lx1=size(nsall,1)-4
    lx2all=size(nsall,2)-4
    lx3all=size(nsall,3)-4

    allocate(Phislab(1:lx2all,1:lx3all))  !space to store EFL potential

    !> READ IN FROM FILE, AS OF CURVILINEAR BRANCH THIS IS NOW THE ONLY INPUT OPTION
    call get_simsize3(indatsize, lx1in, lx2in, lx3in)
    print '(2A,3I6)', indatsize,' input dimensions:',lx1in,lx2in,lx3in
    print '(A,3I6)', 'Target (output) grid structure dimensions:',lx1,lx2all,lx3all

    if (.not. (lx1==lx1in .and. lx2all==lx2in .and. lx3all==lx3in)) then
      error stop 'ERROR:gemini3d: The input data must be the same size as the grid which you are running the simulation on' // &
           '- use a script to interpolate up/down to the simulation grid'
    end if

    call hf%open(indatfile, action='r')

    call hf%read('/nsall', nsall(1:lx1,1:lx2all,1:lx3all,1:lsp))
    call hf%read('/vs1all', vs1all(1:lx1,1:lx2all,1:lx3all,1:lsp))
    call hf%read('/Tsall', Tsall(1:lx1,1:lx2all,1:lx3all,1:lsp))
    if (hf%exist('/Phiall')) then
      if (hf%ndim('/Phiall') == 1) then
        if (lx2all==1) then
          allocate(tmpPhi(lx3all))
        else
          allocate(tmpPhi(lx2all))
        end if
        call hf%read('/Phiall', tmpPhi)
        ! FIXME: MH please delete if you are okay with this
        !if (size(Phislab, 1) /= 1) then
        !  write(stderr,*) 'Phislab shape',shape(Phislab)
        !  error stop 'Phislab x2 /= 1'
        !endif
        if (lx2all==1) then
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

    call hf%close()

    !> Apply EFL approx to compute full grid potential
    do ix1=1,lx1
      Phiall(ix1,1:lx2all,1:lx3all)=Phislab(1:lx2all,1:lx3all)
    end do

    deallocate(Phislab)    ! explicitly get rid of allocated storage
  end procedure getICs_hdf5
end submodule plasma_input_hdf5
