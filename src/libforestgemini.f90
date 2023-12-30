!> These are extensions to libgemini that are used for forestgemini
module forestgemini

implicit none (type, external)

public :: checkE1,grid_from_extents_in,interp_file2subgrid_in,forceZOH_all, permute_fluidvars, &
            ipermute_fluidvars, tag4refine_location, tag4refine_vperp, clean_param_after_regrid_in, get_locationsi_in, &
            set_datainow_in, get_datainow_ptr_in, swap_statevars, interp3_in, interp2_in, tag4refine_diff, &
            tag4refine_grad, tag4coarsening_diff, gemini_grid_alloc

contains
  !> sequence of calls to allocate space for grid object and internal variables (analogous to gemini_work_alloc)
  subroutine gemini_grid_alloc(x1lims,x2lims,x3lims,lx1wg,lx2wg,lx3wg,x,xtype,xC)
    real(wp), dimension(2), intent(in) :: x1lims,x2lims,x3lims
    integer, intent(in) :: lx1wg,lx2wg,lx3wg
    class(curvmesh), intent(inout), pointer :: x
    integer, intent(inout) :: xtype
    type(c_ptr), intent(inout) :: xC
    integer :: ix2,ix3
    real(wp), dimension(:), allocatable :: x2,x3
    real(wp), dimension(:), pointer :: x1
    real(wp) :: glonctr,glatctr

    ! retrieve data from grid module
    call get_gridcenter(glonctr,glatctr)

    ! error checking
    if (glatctr<-90._wp .or. glatctr>90._wp) then
      error stop ' grid_from_extents:  prior to calling must use read_size_gridcenter or set_size_gridcenter to assign &
                   &module variables glonctr,glatctr'
    end if

    ! create temp space
    allocate(x2(lx2wg),x3(lx3wg))

    ! make uniformly spaced coordinate arrays; unless a x1 array was already provided by user
    call get_x1coords(x1)
    !x1=[(x1lims(1) + (x1lims(2)-x1lims(1))/(lx1wg-1)*(ix1-1),ix1=1,lx1wg)]
    x2=[(x2lims(1) + (x2lims(2)-x2lims(1))/(lx2wg-1)*(ix2-1),ix2=1,lx2wg)]
    x3=[(x3lims(1) + (x3lims(2)-x3lims(1))/(lx3wg-1)*(ix3-1),ix3=1,lx3wg)]

    ! allocate mesh class and create a C pointer to it
    call meshobj_alloc(x1,x2,x3,x,xtype,xC)

    ! allocate grid internal data arrays/structs
    call grid_internaldata_alloc(x1,x2,x3,x2,x3,glonctr,glatctr,x)

    ! get rid of temp. arrays
    deallocate(x2,x3)
  end subroutine gemini_grid_alloc


  !> Have each worker read the entire input file and then interpolate it onto its own subgrid
  subroutine interp_file2subgrid_in(cfg,x,fluidvars,electrovars)
    type(gemini_cfg), intent(in) :: cfg
    class(curvmesh), intent(in) :: x
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars,electrovars
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:), pointer :: E1,E2,E3,J1,J2,J3,Phi
    logical :: errflag=.false.
    real(wp) :: nlower,nupper,vlower,vupper,Tlower,Tupper

    print*, 'Initiating tiled/interpolate input...'
    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call electrovar_pointers(electrovars,E1,E2,E3,J1,J2,J3,Phi)
    call interp_file2subgrid(cfg%indatsize,cfg%indatfile,cfg%indatgrid,x%x1,x%x2,x%x3,ns,vs1,Ts,Phi)

    ! it's important to note that the file input *does not specify* vs2,vs3 so we need to set some initial
    !   values here, otherwise we may just get existing garbage in memory.
    vs2=0._wp; vs3=0._wp;

    ! the input code also does not assign electrodynamic variables
    E1=0._wp; E2=0._wp; E3=0._wp; J1=0._wp; J2=0._wp; J3=0._wp;

    ! Clean the input data since it has been interpolated and there is a chance we get unacceptable values
    !   This appears to sometimes be needed with dipole grid simulations.
    call clean_param(x,1,ns)
    call clean_param(x,2,vs1)
    call clean_param(x,3,Ts)

    ! this is a good time to do some error checking since each patch only calls this code once per simulation
    nlower=0; nupper=1e14;
    vlower=-1e4; vupper=1e4;
    Tlower=0; Tupper=20000;
    print*, 'In fortran file i/o...'
    errflag=errflag .or. checkarray(ns,nlower,nupper,'>>> Full density corrupted:  ',0)
    errflag=errflag .or. checkarray(vs1,vlower,vupper,'>>> Full velocity corrupted:  ',0)
    errflag=errflag .or. checkarray(Ts,Tlower,Tupper,'>>> Full temperature corrupted:  ',0)
    if (errflag) error stop
    print*, 'Exiting fortran file i/o...'
  end subroutine interp_file2subgrid_in


  !> A somewhat superfluous wrapper for grid generation from known extents, included here to keep with the
  !    pattern of only having applications access procedures from this module and no others.  In the case
  !    of a Cartesian grid we need to somehow set the center point so that there are glat/glon associated
  !    with each location of the mesh.  This will be obtained from the input cfg file in cases where we are
  !    generating a grid from extents.  Note that this is distinctly different from the situation where we
  !    are using read_grid() to input a grid from a file - in that case the parameters glonctr and glatctr
  !    are expected to be kept within that file.
  subroutine grid_from_extents_in(x1lims,x2lims,x3lims,lx1wg,lx2wg,lx3wg,x)
    real(wp), dimension(2), intent(in) :: x1lims,x2lims,x3lims
    integer, intent(in) :: lx1wg,lx2wg,lx3wg
    class(curvmesh), intent(inout) :: x

    call grid_from_extents(x1lims,x2lims,x3lims,lx1wg,lx2wg,lx3wg,x)
  end subroutine grid_from_extents_in


  !> print out min/max values for variables
  subroutine checkE1(fluidvars,fluidauxvars,electrovars,locID)
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidauxvars
    real(wp), dimension(:,:,:,:), pointer, intent(in) :: electrovars
    integer, intent(in) :: locID
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:,:), pointer :: rhovs1,rhoes
    real(wp), dimension(:,:,:), pointer :: rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom
    real(wp), dimension(:,:,:),pointer :: E1,E2,E3,J1,J2,J3,Phi
    real(wp) :: nlower,nupper,vlower,vupper,Tlower,Tupper
    real(wp) :: vplower,vpupper
    real(wp) :: Eparlower,Eparupper,Elower,Eupper,Jlower,Jupper,Philower,Phiupper
    real(wp) :: rhovlower,rhovupper,rhoeslower,rhoesupper,Blower,Bupper
    logical :: errflag=.false.
    integer :: funit

    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)
    call electrovar_pointers(electrovars,E1,E2,E3,J1,J2,J3,Phi)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Check main variables
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! force a check on the interior cells of the domain
    nlower=0; nupper=1e14;
    vlower=-1e4; vupper=1e4;
    vplower=-1e4; vpupper=1e4;
    Tlower=0; Tupper=10000;

    errflag=errflag .or. checkarray(ns(3:lx1+2,3:lx2+2,3:lx3+2,:),nlower,nupper, &
                                     '>>> Interior density data corrupted:  ',locID)
    errflag=errflag .or. checkarray(vs1(3:lx1+2,3:lx2+2,3:lx3+2,:),vlower,vupper, &
                                     '>>> Interior velocity data corrupted:  ',locID)
    errflag=errflag .or. checkarray(Ts(3:lx1+2,3:lx2+2,3:lx3+2,:),Tlower,Tupper, &
                                     '>>> Interior temperature data corrupted:  ',locID)
    errflag=errflag .or. checkarray(vs2(3:lx1+2,3:lx2+2,3:lx3+2,:),vplower,vpupper, &
                                     '>>> Interior v2 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(vs3(3:lx1+2,3:lx2+2,3:lx3+2,:),vplower,vpupper, &
                                     '>>> Interior v3 data corrupted:  ',locID)

    ! now check the bottom ghost cells
    errflag=errflag .or. checkarray(ns(1:2,:,:,:),nlower,nupper, &
                                     '>>> Bottom density data corrupted:  ',locID)
    errflag=errflag .or. checkarray(vs1(1:2,:,:,:),vlower,vupper, &
                                     '>>> Bottom velocity data corrupted:  ',locID)
    errflag=errflag .or. checkarray(Ts(1:2,:,:,:),Tlower,Tupper, &
                                     '>>> Bottom temperature data corrupted:  ',locID)
    errflag=errflag .or. checkarray(vs2(1:2,:,:,:),vplower,vpupper, &
                                     '>>> Bottom v2 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(vs3(1:2,:,:,:),vplower,vpupper, &
                                     '>>> Bottom v3 data corrupted:  ',locID)

    ! check top
    errflag=errflag .or. checkarray(ns(lx1+3:lx1+4,:,:,:),nlower,nupper, &
                                     '>>> Top density data corrupted:  ',locID)
    errflag=errflag .or. checkarray(vs1(lx1+3:lx1+4,:,:,:),vlower,vupper, &
                                     '>>> Top velocity data corrupted:  ',locID)
    errflag=errflag .or. checkarray(Ts(lx1+3:lx1+4,:,:,:),Tlower,Tupper, &
                                     '>>> Top temperature data corrupted:  ',locID)
    errflag=errflag .or. checkarray(vs2(lx1+3:lx1+4,:,:,:),vplower,vpupper, &
                                     '>>> Top v2 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(vs3(lx1+3:lx1+4,:,:,:),vplower,vpupper, &
                                     '>>> Top v3 data corrupted:  ',locID)

    ! check left
    errflag=errflag .or. checkarray(ns(:,1:2,:,:),nlower,nupper, &
                                     '>>> Left density data corrupted:  ',locID)
    errflag=errflag .or. checkarray(vs1(:,1:2,:,:),vlower,vupper, &
                                     '>>> Left velocity data corrupted:  ',locID)
    errflag=errflag .or. checkarray(Ts(:,1:2,:,:),Tlower,Tupper, &
                                     '>>> Left temperature data corrupted:  ',locID)
    errflag=errflag .or. checkarray(vs2(:,1:2,:,:),vplower,vpupper, &
                                     '>>> Left v2 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(vs3(:,1:2,:,:),vplower,vpupper, &
                                     '>>> Left v3 data corrupted:  ',locID)

    ! check right
    errflag=errflag .or. checkarray(ns(:,lx2+3:lx2+4,:,:),nlower,nupper, &
                                     '>>> Right density data corrupted:  ',locID)
    errflag=errflag .or. checkarray(vs1(:,lx2+3:lx2+4,:,:),vlower,vupper, &
                                     '>>> Right velocity data corrupted:  ',locID)
    errflag=errflag .or. checkarray(Ts(:,lx2+3:lx2+4,:,:),Tlower,Tupper, &
                                     '>>> Right temperature data corrupted:  ',locID)
    errflag=errflag .or. checkarray(vs2(:,lx2+3:lx2+4,:,:),vplower,vpupper, &
                                     '>>> Right v2 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(vs3(:,lx2+3:lx2+4,:,:),vplower,vpupper, &
                                     '>>> Right v3 data corrupted:  ',locID)

    ! check bwd
    errflag=errflag .or. checkarray(ns(:,:,1:2,:),nlower,nupper, &
                                     '>>> Bwd density data corrupted:  ',locID)
    errflag=errflag .or. checkarray(vs1(:,:,1:2,:),vlower,vupper, &
                                     '>>> Bwd velocity data corrupted:  ',locID)
    errflag=errflag .or. checkarray(Ts(:,:,1:2,:),Tlower,Tupper, &
                                     '>>> Bwd temperature data corrupted:  ',locID)
    errflag=errflag .or. checkarray(vs2(:,:,1:2,:),vplower,vpupper, &
                                     '>>> Bwd v2 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(vs3(:,:,1:2,:),vplower,vpupper, &
                                     '>>> Bwd v3 data corrupted:  ',locID)

    ! check fwd
    errflag=errflag .or. checkarray(ns(:,:,lx3+3:lx3+4,:),nlower,nupper, &
                                     '>>> Fwd density data corrupted:  ',locID)
    errflag=errflag .or. checkarray(vs1(:,:,lx3+3:lx3+4,:),vlower,vupper, &
                                     '>>> Fwd velocity data corrupted:  ',locID)
    errflag=errflag .or. checkarray(Ts(:,:,lx3+3:lx3+4,:),Tlower,Tupper, &
                                     '>>> Fwd temperature data corrupted:  ',locID)
    errflag=errflag .or. checkarray(vs2(:,:,lx3+3:lx3+4,:),vplower,vpupper, &
                                     '>>> Fwd v2 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(vs3(:,:,lx3+3:lx3+4,:),vplower,vpupper, &
                                     '>>> Fwd v3 data corrupted:  ',locID)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Check electro variables
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Eparlower=-1e-4; Eparupper=1e-4;
    Elower=-0.5; Eupper=0.5;
    Jlower=-1e-3; Jupper=1e-3;
    Philower=-1e5; Phiupper=1e5;

    errflag=errflag .or. checkarray(E1(3:lx1+2,3:lx2+2,3:lx3+2),Eparlower,Eparupper, &
                                     '>>> Interior E1 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(E2(3:lx1+2,3:lx2+2,3:lx3+2),Elower,Eupper, &
                                     '>>> Interior E2 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(E3(3:lx1+2,3:lx2+2,3:lx3+2),Elower,Eupper, &
                                     '>>> Interior E3 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(J1(3:lx1+2,3:lx2+2,3:lx3+2),Jlower,Jupper, &
                                     '>>> Interior J1 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(J2(3:lx1+2,3:lx2+2,3:lx3+2),Jlower,Jupper, &
                                     '>>> Interior J2 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(J3(3:lx1+2,3:lx2+2,3:lx3+2),Jlower,Jupper, &
                                     '>>> Interior J3 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(Phi(3:lx1+2,3:lx2+2,3:lx3+2),Philower,Phiupper, &
                                     '>>> Interior Phi data corrupted:  ',locID)

    errflag=errflag .or. checkarray(E1(1:2,:,:),Eparlower,Eparupper, &
                                     '>>> Bottom E1 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(E2(1:2,:,:),Elower,Eupper, &
                                     '>>> Bottom E2 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(E3(1:2,:,:),Elower,Eupper, &
                                     '>>> Bottom E3 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(J1(1:2,:,:),Jlower,Jupper, &
                                     '>>> Bottom J1 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(J2(1:2,:,:),Jlower,Jupper, &
                                     '>>> Bottom J2 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(J3(1:2,:,:),Jlower,Jupper, &
                                     '>>> Bottom J3 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(Phi(1:2,:,:),Philower,Phiupper, &
                                     '>>> Bottom Phi data corrupted:  ',locID)

    errflag=errflag .or. checkarray(E1(lx1+3:lx1+4,:,:),Eparlower,Eparupper, &
                                     '>>> Top E1 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(E2(lx1+3:lx1+4,:,:),Elower,Eupper, &
                                     '>>> Top E2 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(E3(lx1+3:lx1+4,:,:),Elower,Eupper, &
                                     '>>> Top E3 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(J1(lx1+3:lx1+4,:,:),Jlower,Jupper, &
                                     '>>> Top J1 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(J2(lx1+3:lx1+4,:,:),Jlower,Jupper, &
                                     '>>> Top J2 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(J3(lx1+3:lx1+4,:,:),Jlower,Jupper, &
                                     '>>> Top J3 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(Phi(lx1+3:lx1+4,:,:),Philower,Phiupper, &
                                     '>>> Top Phi data corrupted:  ',locID)

    errflag=errflag .or. checkarray(E1(:,1:2,:),Eparlower,Eparupper, &
                                     '>>> Left E1 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(E2(:,1:2,:),Elower,Eupper, &
                                     '>>> Left E2 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(E3(:,1:2,:),Elower,Eupper, &
                                     '>>> Left E3 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(J1(:,1:2,:),Jlower,Jupper, &
                                     '>>> Left J1 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(J2(:,1:2,:),Jlower,Jupper, &
                                     '>>> Left J2 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(J3(:,1:2,:),Jlower,Jupper, &
                                     '>>> Left J3 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(Phi(:,1:2,:),Philower,Phiupper, &
                                     '>>> Left Phi data corrupted:  ',locID)

    errflag=errflag .or. checkarray(E1(:,lx2+3:lx2+4,:),Eparlower,Eparupper, &
                                     '>>> Right E1 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(E2(:,lx2+3:lx2+4,:),Elower,Eupper, &
                                     '>>> Right E2 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(E3(:,lx2+3:lx2+4,:),Elower,Eupper, &
                                     '>>> Right E3 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(J1(:,lx2+3:lx2+4,:),Jlower,Jupper, &
                                     '>>> Right J1 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(J2(:,lx2+3:lx2+4,:),Jlower,Jupper, &
                                     '>>> Right J2 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(J3(:,lx2+3:lx2+4,:),Jlower,Jupper, &
                                     '>>> Right J3 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(Phi(:,lx2+3:lx2+4,:),Philower,Phiupper, &
                                     '>>> Right Phi data corrupted:  ',locID)

    errflag=errflag .or. checkarray(E1(:,:,1:2),Eparlower,Eparupper, &
                                     '>>> Bwd E1 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(E2(:,:,1:2),Elower,Eupper, &
                                     '>>> Bwd E2 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(E3(:,:,1:2),Elower,Eupper, &
                                     '>>> Bwd E3 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(J1(:,:,1:2),Jlower,Jupper, &
                                     '>>> Bwd J1 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(J2(:,:,1:2),Jlower,Jupper, &
                                     '>>> Bwd J2 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(J3(:,:,1:2),Jlower,Jupper, &
                                     '>>> Bwd J3 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(Phi(:,:,1:2),Philower,Phiupper, &
                                     '>>> Bwd Phi data corrupted:  ',locID)

    errflag=errflag .or. checkarray(E1(:,:,lx3+3:lx3+4),Eparlower,Eparupper, &
                                     '>>> Fwd E1 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(E2(:,:,lx3+3:lx3+4),Elower,Eupper, &
                                     '>>> Fwd E2 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(E3(:,:,lx3+3:lx3+4),Elower,Eupper, &
                                     '>>> Fwd E3 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(J1(:,:,lx3+3:lx3+4),Jlower,Jupper, &
                                     '>>> Fwd J1 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(J2(:,:,lx3+3:lx3+4),Jlower,Jupper, &
                                     '>>> Fwd J2 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(J3(:,:,lx3+3:lx3+4),Jlower,Jupper, &
                                     '>>> Fwd J3 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(Phi(:,:,lx3+3:lx3+4),Jlower,Jupper, &
                                     '>>> Fwd Phi data corrupted:  ',locID)

    errflag=errflag .or. checkarray(E1(:,:,lx3+3:lx3+4),Eparlower,Eparupper, &
                                     '>>> Fwd E1 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(E2(:,:,lx3+3:lx3+4),Elower,Eupper, &
                                     '>>> Fwd E2 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(E3(:,:,lx3+3:lx3+4),Elower,Eupper, &
                                     '>>> Fwd E3 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(J1(:,:,lx3+3:lx3+4),Jlower,Jupper, &
                                     '>>> Fwd J1 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(J2(:,:,lx3+3:lx3+4),Jlower,Jupper, &
                                     '>>> Fwd J2 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(J3(:,:,lx3+3:lx3+4),Jlower,Jupper, &
                                     '>>> Fwd J3 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(Phi(:,:,lx3+3:lx3+4),Philower,Phiupper, &
                                     '>>> Fwd Phi data corrupted:  ',locID)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Check aux variables
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    rhoeslower=0; rhoesupper=1;
    rhovlower=-1; rhovupper=1;
    Blower=-100000e-9; Bupper=100000e-9;

    errflag=errflag .or. checkarray(rhovs1(3:lx1+2,3:lx2+2,3:lx3+2,:),rhovlower,rhovupper, &
                                     '>>> Interior rhovs1 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(rhoes(3:lx1+2,3:lx2+2,3:lx3+2,:),rhoeslower,rhoesupper, &
                                     '>>> Interior rhoes data corrupted:  ',locID)
    errflag=errflag .or. checkarray(B1(3:lx1+2,3:lx2+2,3:lx3+2),Blower,Bupper, &
                                     '>>> Interior B1 data corrupted:  ',locID)

    errflag=errflag .or. checkarray(rhovs1(1:2,:,:,:),rhovlower,rhovupper, &
                                     '>>> Bottom rhovs1 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(rhoes(1:2,:,:,:),rhoeslower,rhoesupper, &
                                     '>>> Bottom rhoes data corrupted:  ',locID)
    errflag=errflag .or. checkarray(B1(1:2,:,:),Blower,Bupper, &
                                     '>>> Bottom B1 data corrupted:  ',locID)

    errflag=errflag .or. checkarray(rhovs1(lx1+3:lx1+4,:,:,:),rhovlower,rhovupper, &
                                     '>>> Top rhovs1 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(rhoes(lx1+3:lx1+4,:,:,:),rhoeslower,rhoesupper, &
                                     '>>> Top rhoes data corrupted:  ',locID)
    errflag=errflag .or. checkarray(B1(lx1+3:lx1+4,:,:),Blower,Bupper, &
                                     '>>> Top B1 data corrupted:  ',locID)

    errflag=errflag .or. checkarray(rhovs1(:,1:2,:,:),rhovlower,rhovupper, &
                                     '>>> Left rhovs1 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(rhoes(:,1:2,:,:),rhoeslower,rhoesupper, &
                                     '>>> Left rhoes data corrupted:  ',locID)
    errflag=errflag .or. checkarray(B1(:,1:2,:),Blower,Bupper, &
                                     '>>> Left B1 data corrupted:  ',locID)

    errflag=errflag .or. checkarray(rhovs1(:,lx2+3:lx2+4,:,:),rhovlower,rhovupper, &
                                     '>>> Right rhovs1 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(rhoes(:,lx2+3:lx2+4,:,:),rhoeslower,rhoesupper, &
                                     '>>> Right rhoes data corrupted:  ',locID)
    errflag=errflag .or. checkarray(B1(:,lx2+3:lx2+4,:),Blower,Bupper, &
                                     '>>> Right B1 data corrupted:  ',locID)

    errflag=errflag .or. checkarray(rhovs1(:,:,1:2,:),rhovlower,rhovupper, &
                                     '>>> Bwd rhovs1 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(rhoes(:,:,1:2,:),rhoeslower,rhoesupper, &
                                     '>>> Bwd rhoes data corrupted:  ',locID)
    errflag=errflag .or. checkarray(B1(:,:,1:2),Blower,Bupper, &
                                     '>>> Bwd B1 data corrupted:  ',locID)

    errflag=errflag .or. checkarray(rhovs1(:,:,lx3+3:lx3+4,:),rhovlower,rhovupper, &
                                     '>>> Fwd rhovs1 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(rhoes(:,:,lx3+3:lx3+4,:),rhoeslower,rhoesupper, &
                                     '>>> Fwd rhoes data corrupted:  ',locID)
    errflag=errflag .or. checkarray(B1(:,:,lx3+3:lx3+4),Blower,Bupper, &
                                     '>>> Fwd B1 data corrupted:  ',locID)

!    if (errflag) then
!      open(newunit=funit,file='error.dat',status='replace',access='stream')
!      write(funit) ns
!      write(funit) vs1
!      write(funit) Ts
!      close(funit)
!
!      error stop
!    end if

    ! FIXME: desperate attempt to fix issues following a refine+time step.  This is obviously not a great
    !   solution but it does work to address many/most of the problems I see on the time step after refine.
    !   I still don't know *why* these are happening but in the meantime this allows us to move forward and
    !   do some testing and basic simulations.  
    where (abs(vs1)>1e4)
      vs1=0._wp
    end where
    where (abs(vs2)>1e4)
      vs2=0._wp
    end where
    where (abs(vs3)>1e4)
      vs3=0._wp
    end where
    where (Ts>1.e4)
      Ts=1.e4
    end where
    print*, minval(vs1),maxval(vs1),minval(vs2),maxval(vs2),minval(vs3),maxval(vs3),minval(Ts),maxval(Ts)
  end subroutine checkE1


  !> just force a zero-order hold for the ghost cells, should probably only be used for debugging purposes
  subroutine forceZOH(param)
    real(wp), dimension(-1:,-1:,-1:,:), intent(inout) :: param
    integer :: lx1,lx2,lx3

    lx1=size(param,1)-4; lx2=size(param,2)-4; lx3=size(param,3)-4;

!    param(0,:,:,:)=param(1,:,:,:)
!    param(-1,:,:,:)=param(1,:,:,:)
!    param(lx1+1,:,:,:)=param(lx1,:,:,:)
!    param(lx1+2,:,:,:)=param(lx1,:,:,:)
!
!    param(:,0,:,:)=param(:,1,:,:)
!    param(:,-1,:,:)=param(:,1,:,:)
!    param(:,lx2+1,:,:)=param(:,lx2,:,:)
!    param(:,lx2+2,:,:)=param(:,lx2,:,:)

    param(:,:,0,:)=param(:,:,1,:)
    param(:,:,-1,:)=param(:,:,1,:)
    param(:,:,lx3+1,:)=param(:,:,lx3,:)
    param(:,:,lx3+2,:)=param(:,:,lx3,:)
  end subroutine forceZOH


  !> force all primary variables to have a zero-order hold extrapolation in ghost cells for testing purposes
  subroutine forceZOH_all(fluidvars)
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts

    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call forceZOH(ns)
    call forceZOH(vs1)
    call forceZOH(vs2)
    call forceZOH(vs3)
    call forceZOH(Ts)
  end subroutine forceZOH_all


  !> permute state variables x1,x2,x3 --> x2,x3,x1.  FIXME: this is now deprecated since we also need to 
  !    swap for density variables when using forestclaw.
  subroutine permute_fluidvars(fluidvars)
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars
    integer lx1,lx2,lx3,leqn,ix1,ix2,ix3,ieqn
    real(wp), dimension(:,:,:,:), allocatable :: fluidvarsT

    ! size and allocate
    lx1=size(fluidvars,1)-4; lx2=size(fluidvars,2)-4; lx3=size(fluidvars,3)-4; leqn=size(fluidvars,4);
    allocate(fluidvarsT(-1:lx2+2,-1:lx3+2,-1:lx1+2,leqn))

    ! place data in permuted array
    do ieqn=1,leqn
      do ix3=-1,lx3+2
        do ix2=-1,lx2+2
          do ix1=-1,lx1+2
            fluidvarsT(ix2,ix3,ix1,ieqn)=fluidvars(ix1+2,ix2+2,ix3+2,ieqn)
          end do
        end do
      end do
    end do

    ! reshape the permuted array to match original, copy into state variable array
    fluidvars=reshape(fluidvarsT,[lx1+4,lx2+4,lx3+4,leqn])

    ! explicitly clear memory
    deallocate(fluidvarsT)
  end subroutine permute_fluidvars


  ! inverse permute state variables x2,x3,x1 --> x1,x2,x3.  FIXME: also deprecated.  
  subroutine ipermute_fluidvars(fluidvars)
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars
    integer lx1,lx2,lx3,leqn,ix1,ix2,ix3,ieqn
    real(wp), dimension(:,:,:,:), allocatable :: fluidvarsT

    ! size and allocate, sizes of input array will be correct but data will be shaped wrongly
    lx1=size(fluidvars,1)-4; lx2=size(fluidvars,2)-4; lx3=size(fluidvars,3)-4; leqn=size(fluidvars,4);
    allocate(fluidvarsT(-1:lx2+2,-1:lx3+2,-1:lx1+2,leqn))
    fluidvarsT=reshape(fluidvars,[lx2+4,lx3+4,lx1+4,leqn])

    ! place data in permuted array
    do ieqn=1,leqn
      do ix3=-1,lx3+2
        do ix2=-1,lx2+2
          do ix1=-1,lx1+2
            fluidvars(ix1+2,ix2+2,ix3+2,ieqn)=fluidvarsT(ix2,ix3,ix1,ieqn)
          end do
        end do
      end do
    end do

    ! explicitly clear memory
    deallocate(fluidvarsT)
  end subroutine ipermute_fluidvars


  !> Tag for refinement based on location
  subroutine tag4refine_location(x,flagrefine)
    class(curvmesh), intent(in) :: x
    logical, intent(inout) :: flagrefine
    integer :: ix1,ix2,ix3
    real(wp) :: mlat,mlon,alt

    flagrefine=.false.
    do ix3=1,x%lx3
      do ix2=1,x%lx2
        do ix1=1,x%lx1
          mlon=x%phi(ix1,ix2,ix3)*180.0/pi
          mlat=90.0-x%theta(ix1,ix2,ix3)*180.0/pi
          alt=x%alt(ix1,ix2,ix3)
!          if ( mlon > 209.5 .and. mlon < 210.5 .and. mlat > 28.0 .and. mlat < 29.0 .and. &
!               alt > 80e3 .and. alt < 300e3) then
!            flagrefine=.true.
!            exit
!          end if
          if ( sqrt( ((mlon-210)/3)**2 + (mlat-28.5)**2) < 0.35 .and. &
               alt > 80e3 .and. alt < 300e3) then
            flagrefine=.true.
            exit
          end if
        end do
      end do
    end do
  end subroutine tag4refine_location


  !> Tag for refinement based on perpendicular velocity
  subroutine tag4refine_vperp(x,fluidvars,fluidauxvars,electrovars,intvars,flagrefine)
    class(curvmesh), intent(in) :: x
    real(wp), dimension(:,:,:,:), pointer, intent(in) :: fluidvars,fluidauxvars,electrovars
    type(gemini_work) :: intvars
    logical, intent(inout) :: flagrefine
    integer :: ix1,ix2,ix3
    real(wp) :: mlat,mlon,alt
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:,:), pointer :: rhovs1,rhoes
    real(wp), dimension(:,:,:), pointer :: rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom
    real(wp), dimension(:,:,:), pointer :: E1,E2,E3,J1,J2,J3,Phi
    real(wp) :: vperp,vpar

    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)
    call electrovar_pointers(electrovars,E1,E2,E3,J1,J2,J3,Phi)


    flagrefine=.false.
    do ix3=1,x%lx3
      do ix2=1,x%lx2
        do ix1=1,x%lx1
          mlon=x%phi(ix1,ix2,ix3)*180.0/pi
          mlat=90.0-x%theta(ix1,ix2,ix3)*180.0/pi
          alt=x%alt(ix1,ix2,ix3)
          !vperp=sqrt(vs2(ix1,ix2,ix3,1)**2 + vs3(ix1,ix2,ix3,1)**2)    ! use the major ion species
          !vpar=abs(vs1(ix1,ix2,ix3,1))
          vpar=abs(intvars%atmos%vn1(ix1,ix2,ix3))
          if (alt > 90e3 .and. alt < 350e3 .and. vpar > 25._wp) then     ! more than 50 m/s probably  means something is happening
            flagrefine=.true.
            exit
          end if
        end do
      end do
    end do
  end subroutine tag4refine_vperp


  !> Tag for refinement based on differences within some range
!  subroutine tag4refine_diff(x,fluidvars,fluidauxvars,electrovars,intvars,flagrefine)
!    class(curvmesh), intent(in) :: x
!    real(wp), dimension(:,:,:,:), pointer, intent(in) :: fluidvars,fluidauxvars,electrovars
!    type(gemini_work) :: intvars
!    logical, intent(inout) :: flagrefine
!    integer :: ix1,ix2,ix3    
!    real(wp) :: mlat,mlon,alt
!    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
!    real(wp), dimension(:,:,:,:), pointer :: rhovs1,rhoes
!    real(wp), dimension(:,:,:), pointer :: rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom
!    real(wp), dimension(:,:,:), pointer :: E1,E2,E3,J1,J2,J3,Phi
!    real(wp) :: minv,maxv,deltav
!    
!    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
!    call fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)
!    call electrovar_pointers(electrovars,E1,E2,E3,J1,J2,J3,Phi)
!
!    flagrefine=.false.
!    deltav=0.0
!    minv=0.0
!    maxv=0.0
!    do ix3=1,x%lx3
!      do ix2=1,x%lx2
!        do ix1=1,x%lx1
!          mlon=x%phi(ix1,ix2,ix3)*180.0/pi
!          mlat=90.0-x%theta(ix1,ix2,ix3)*180.0/pi
!          alt=x%alt(ix1,ix2,ix3)
!          if (sqrt( ((mlon-210)/3)**2 + (mlat-28.5)**2) < 0.35 .and. &
!               alt > 80e3 .and. alt < 300e3) then      ! only update min/max v if in region of interest
!            if (vs1(ix1,ix2,ix3,7) < minv) minv=vs1(ix1,ix2,ix3,7)
!            if (vs1(ix1,ix2,ix3,7) > maxv) maxv=vs1(ix1,ix2,ix3,7)
!          end if
!        end do
!      end do
!    end do
!    deltav=maxv-minv
!
!    if (deltav > 10.0) flagrefine=.true.
!  end subroutine tag4refine_diff

!  subroutine tag4refine_diff(x,fluidvars,fluidauxvars,electrovars,intvars,flagrefine)
!    class(curvmesh), intent(in) :: x
!    real(wp), dimension(:,:,:,:), pointer, intent(in) :: fluidvars,fluidauxvars,electrovars
!    type(gemini_work) :: intvars
!    logical, intent(inout) :: flagrefine
!    integer :: ix1,ix2,ix3    
!    real(wp) :: mlat,mlon,alt
!    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
!    real(wp), dimension(:,:,:,:), pointer :: rhovs1,rhoes
!    real(wp), dimension(:,:,:), pointer :: rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom
!    real(wp), dimension(:,:,:), pointer :: E1,E2,E3,J1,J2,J3,Phi
!    real(wp) :: minv,maxv,deltav
!    
!    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
!    call fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)
!    call electrovar_pointers(electrovars,E1,E2,E3,J1,J2,J3,Phi)
!
!    flagrefine=.false.
!    deltav=0.0
!    minv=0.0
!    maxv=0.0
!    do ix3=1,x%lx3
!      do ix2=1,x%lx2
!        do ix1=1,x%lx1
!          mlon=x%phi(ix1,ix2,ix3)*180.0/pi
!          mlat=90.0-x%theta(ix1,ix2,ix3)*180.0/pi
!          alt=x%alt(ix1,ix2,ix3)
!          if (alt > 80e3 .and. alt < 300e3) then      ! only update min/max v if in region of interest
!            if (intvars%atmos%vn1(ix1,ix2,ix3) < minv) minv=intvars%atmos%vn1(ix1,ix2,ix3)
!            if (intvars%atmos%vn1(ix1,ix2,ix3) > maxv) maxv=intvars%atmos%vn1(ix1,ix2,ix3)
!          end if
!        end do
!      end do
!    end do
!    deltav=maxv-minv
!
!    if (deltav > 5.0) flagrefine=.true.
!  end subroutine tag4refine_diff


  subroutine tag4refine_diff(x,fluidvars,fluidauxvars,electrovars,intvars,flagrefine)
    class(curvmesh), intent(in) :: x
    real(wp), dimension(:,:,:,:), pointer, intent(in) :: fluidvars,fluidauxvars,electrovars
    type(gemini_work) :: intvars
    logical, intent(inout) :: flagrefine
    integer :: ix1,ix2,ix3    
    real(wp) :: mlat,mlon,alt
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:,:), pointer :: rhovs1,rhoes
    real(wp), dimension(:,:,:), pointer :: rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom
    real(wp), dimension(:,:,:), pointer :: E1,E2,E3,J1,J2,J3,Phi
    real(wp) :: minv,maxv,deltav
    real(wp) :: minvinner,maxvinner,deltavinner
    
    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)
    call electrovar_pointers(electrovars,E1,E2,E3,J1,J2,J3,Phi)

    flagrefine=.false.
    deltav=0.0
    minv=0.0
    maxv=0.0
    do ix1=1,x%lx1
      deltavinner=0.0
      minvinner=0.0
      maxvinner=0.0
      do ix2=1,x%lx2
        do ix3=1,x%lx3
          mlon=x%phi(ix1,ix2,ix3)*180.0/pi
          mlat=90.0-x%theta(ix1,ix2,ix3)*180.0/pi
          alt=x%alt(ix1,ix2,ix3)
          if (alt > 80e3 .and. alt < 300e3) then      ! only update min/max v if in region of interest
            if (intvars%atmos%vn1(ix1,ix2,ix3) < minvinner) minvinner=intvars%atmos%vn1(ix1,ix2,ix3)
            if (intvars%atmos%vn1(ix1,ix2,ix3) > maxvinner) maxvinner=intvars%atmos%vn1(ix1,ix2,ix3)
          end if
        end do
      end do
      deltavinner=maxvinner-minvinner
      if (deltav < deltavinner) deltav=deltavinner
    end do

    if (deltav > 50.0) flagrefine=.true.
  end subroutine tag4refine_diff


  !> Tag for refinement based on differences within some range
  subroutine tag4refine_grad(x,fluidvars,fluidauxvars,electrovars,intvars,flagrefine)
    class(curvmesh), intent(in) :: x
    real(wp), dimension(:,:,:,:), pointer, intent(in) :: fluidvars,fluidauxvars,electrovars
    type(gemini_work) :: intvars
    logical, intent(inout) :: flagrefine
    integer :: ix1,ix2,ix3    
    real(wp) :: mlat,mlon,alt
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:,:), pointer :: rhovs1,rhoes
    real(wp), dimension(:,:,:), pointer :: rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom
    real(wp), dimension(:,:,:), pointer :: E1,E2,E3,J1,J2,J3,Phi
!    real(wp) :: vtest1,vtest2
    real(wp), dimension(x%lx1,x%lx2,x%lx3) :: gradvn12,gradvn13
    
    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)
    call electrovar_pointers(electrovars,E1,E2,E3,J1,J2,J3,Phi)

!    vtest1=minval(intvars%atmos%vn1)
!    vtest2=maxval(intvars%atmos%vn1)
!    if (vtest2-vtest1<15._wp) then
!      flagrefine=.false.
!    else
!      flagrefine=.true.
!    end if

    gradvn12(:,:,:)=grad3D2(intvars%atmos%vn1(:,:,:),x,1,x%lx1,1,x%lx2,1,x%lx3)
    gradvn13(:,:,:)=grad3D3(intvars%atmos%vn1(:,:,:),x,1,x%lx1,1,x%lx2,1,x%lx3)

    flagrefine=.false.
    do ix3=1,x%lx3
      do ix2=1,x%lx2
        do ix1=1,x%lx1
          mlon=x%phi(ix1,ix2,ix3)*180.0/pi
          mlat=90.0-x%theta(ix1,ix2,ix3)*180.0/pi
          alt=x%alt(ix1,ix2,ix3)
          if (alt > 90e3 .and. alt < 350e3 .and. &
                  (abs(gradvn12(ix1,ix2,ix3)) > 1.25e-3 .or. abs(gradvn13(ix1,ix2,ix3)) > 1.25e-3) ) then
            flagrefine=.true.
            exit
          end if
        end do
      end do
    end do
  end subroutine tag4refine_grad


    !> Tag for coarsening based on differences within some range
!  subroutine tag4coarsening_diff(x,fluidvars,fluidauxvars,electrovars,intvars,flagcoarsening)
!    class(curvmesh), intent(in) :: x
!    real(wp), dimension(:,:,:,:), pointer, intent(in) :: fluidvars,fluidauxvars,electrovars
!    type(gemini_work) :: intvars
!    logical, intent(inout) :: flagcoarsening
!    integer :: ix1,ix2,ix3    
!    real(wp) :: mlat,mlon,alt
!    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
!    real(wp), dimension(:,:,:,:), pointer :: rhovs1,rhoes
!    real(wp), dimension(:,:,:), pointer :: rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom
!    real(wp), dimension(:,:,:), pointer :: E1,E2,E3,J1,J2,J3,Phi
!    real(wp) :: minv,maxv,deltav
!    
!    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
!    call fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)
!    call electrovar_pointers(electrovars,E1,E2,E3,J1,J2,J3,Phi)
!
!    flagcoarsening=.false.
!    deltav=0.0
!    minv=0.0
!    maxv=0.0
!    do ix3=1,x%lx3
!      do ix2=1,x%lx2
!        do ix1=1,x%lx1
!          mlon=x%phi(ix1,ix2,ix3)*180.0/pi
!          mlat=90.0-x%theta(ix1,ix2,ix3)*180.0/pi
!          alt=x%alt(ix1,ix2,ix3)
!          if (sqrt( ((mlon-210)/3)**2 + (mlat-28.5)**2) < 0.35 .and. &
!               alt > 80e3 .and. alt < 300e3) then      ! only update min/max v if in region of interest
!            if (vs1(ix1,ix2,ix3,7) < minv) minv=vs1(ix1,ix2,ix3,7)
!            if (vs1(ix1,ix2,ix3,7) > maxv) maxv=vs1(ix1,ix2,ix3,7)
!          end if
!        end do
!      end do
!    end do
!    deltav=maxv-minv
!
!    if (deltav < 5.0) flagcoarsening=.true.
!  end subroutine tag4coarsening_diff

!  subroutine tag4coarsening_diff(x,fluidvars,fluidauxvars,electrovars,intvars,flagcoarsening)
!    class(curvmesh), intent(in) :: x
!    real(wp), dimension(:,:,:,:), pointer, intent(in) :: fluidvars,fluidauxvars,electrovars
!    type(gemini_work) :: intvars
!    logical, intent(inout) :: flagcoarsening
!    integer :: ix1,ix2,ix3    
!    real(wp) :: mlat,mlon,alt
!    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
!    real(wp), dimension(:,:,:,:), pointer :: rhovs1,rhoes
!    real(wp), dimension(:,:,:), pointer :: rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom
!    real(wp), dimension(:,:,:), pointer :: E1,E2,E3,J1,J2,J3,Phi
!    real(wp) :: minv,maxv,deltav
!    
!    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
!    call fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)
!    call electrovar_pointers(electrovars,E1,E2,E3,J1,J2,J3,Phi)
!
!    flagcoarsening=.false.
!    deltav=0.0
!    minv=0.0
!    maxv=0.0
!    do ix3=1,x%lx3
!      do ix2=1,x%lx2
!        do ix1=1,x%lx1
!          mlon=x%phi(ix1,ix2,ix3)*180.0/pi
!          mlat=90.0-x%theta(ix1,ix2,ix3)*180.0/pi
!          alt=x%alt(ix1,ix2,ix3)
!          if (alt > 80e3 .and. alt < 300e3) then      ! only update min/max v if in region of interest
!            if (intvars%atmos%vn1(ix1,ix2,ix3) < minv) minv=intvars%atmos%vn1(ix1,ix2,ix3)
!            if (intvars%atmos%vn1(ix1,ix2,ix3) > maxv) maxv=intvars%atmos%vn1(ix1,ix2,ix3)
!          end if
!        end do
!      end do
!    end do
!    deltav=maxv-minv
!
!    if (deltav < 2.5) flagcoarsening=.true.
!  end subroutine tag4coarsening_diff


    subroutine tag4coarsening_diff(x,fluidvars,fluidauxvars,electrovars,intvars,flagcoarsening)
    class(curvmesh), intent(in) :: x
    real(wp), dimension(:,:,:,:), pointer, intent(in) :: fluidvars,fluidauxvars,electrovars
    type(gemini_work) :: intvars
    logical, intent(inout) :: flagcoarsening
    integer :: ix1,ix2,ix3    
    real(wp) :: mlat,mlon,alt
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:,:), pointer :: rhovs1,rhoes
    real(wp), dimension(:,:,:), pointer :: rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom
    real(wp), dimension(:,:,:), pointer :: E1,E2,E3,J1,J2,J3,Phi
    real(wp) :: minv,maxv,deltav
    real(wp) :: minvinner,maxvinner,deltavinner
    
    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)
    call electrovar_pointers(electrovars,E1,E2,E3,J1,J2,J3,Phi)

    flagcoarsening=.false.
    deltav=0.0
    minv=0.0
    maxv=0.0
    do ix1=1,x%lx1
      deltavinner=0.0
      minvinner=0.0
      maxvinner=0.0
      do ix2=1,x%lx2
        do ix3=1,x%lx3
          mlon=x%phi(ix1,ix2,ix3)*180.0/pi
          mlat=90.0-x%theta(ix1,ix2,ix3)*180.0/pi
          alt=x%alt(ix1,ix2,ix3)
          if (alt > 80e3 .and. alt < 300e3) then      ! only update min/max v if in region of interest
            if (intvars%atmos%vn1(ix1,ix2,ix3) < minvinner) minvinner=intvars%atmos%vn1(ix1,ix2,ix3)
            if (intvars%atmos%vn1(ix1,ix2,ix3) > maxvinner) maxvinner=intvars%atmos%vn1(ix1,ix2,ix3)
          end if
        end do
      end do
      deltavinner=maxvinner-minvinner
      if (deltav < deltavinner) deltav=deltavinner
    end do

    if (deltav < 25.0) flagcoarsening=.true.
  end subroutine tag4coarsening_diff



  !> if refinement is being done it may be advantageous to have refine/interpolate done with drift and temperature
  !    it's easiest to just copy swap existing variables.
  subroutine swap_statevars(fluidvars,fluidauxvars)
    real(wp), dimension(:,:,:,:), pointer, intent(in) :: fluidvars,fluidauxvars
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:,:), pointer :: rhovs1,rhoes
    real(wp), dimension(:,:,:), pointer :: rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom
    real(wp), dimension(:,:,:,:), allocatable :: tmpvar

    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)

    !allocate(tmpvar(-1:size(ns,1)-2,-1:size(ns,2)-2,-1:size(ns,3)-2,1:size(ns,4)))
    allocate(tmpvar,mold=rhovs1)
    tmpvar(:,:,:,:)=rhovs1(:,:,:,:)
    rhovs1(:,:,:,:)=vs1(:,:,:,:)
    vs1(:,:,:,:)=tmpvar(:,:,:,:)

    tmpvar(:,:,:,:)=rhoes(:,:,:,:)
    rhoes(:,:,:,:)=Ts(:,:,:,:)
    Ts(:,:,:,:)=tmpvar(:,:,:,:)
    deallocate(tmpvar)
  end subroutine swap_statevars


!  subroutine tag4refine_dist(x,flagrefine)
!    class(curvmesh), intent(in) :: x
!    logical, intent(inout) :: flagrefine
!    real(wp) :: phi1,theta1
!    real(wp), dimension(1:x%lx1,1:x%lx2,1:x%lx3) :: xi,yi,zi
!    real(wp), parameter :: rmag=25e3
!    real(wp) :: r,xi0,yi0
!    integer :: ix1,ix2,ix3
!
!    ! retrieve grid ENU coords
!    call geog2geomag(x%glonctr,x%glatctr,phi1,theta1)
!    call ECEFspher2ENU(x%alt(1:x%lx1,1:x%lx2,1:x%lx3),x%theta(1:x%lx1,1:x%lx2,1:x%lx3), &
!                         x%phi(1:x%lx1,1:x%lx2,1:x%lx3), &
!                         theta1,phi1,xi,yi,zi)
!
!    ! find the average x,y on the end slab so we can control refinement
!    ix1=x%lx1-25
!    xi0=sum(xi(ix1,:,:))/(x%lx2*x%lx3)
!    yi0=sum(yi(ix1,:,:))/(x%lx2*x%lx3)
!
!    ! now check for distance from center to determine refinment
!    flagrefine=.false.
!    do ix3=1,lx3
!      do ix2=1,lx2
!          r=sqrt( (xi(ix1,ix2,ix3)-xi0)**2 + (yi(ix1,ix2,ix3)-yi0)**2)
!          if (r<rmag) then
!            flagrefine=.true.
!            exit
!          end if
!      end do
!    end do
!  end subroutine tag4refine_dist


  !> Call function to retrieve locations from a neutraldata3D_fclaw object
  subroutine get_locationsi_in(intvars,flagallpts,zlims,xlims,ylims,zvals,xvals,yvals,datavals)
    type(gemini_work), intent(inout) :: intvars
    logical, intent(in) :: flagallpts
    real(wp), dimension(2), intent(in) :: zlims,xlims,ylims
    real(wp), dimension(:), pointer, intent(inout) :: zvals,xvals,yvals
    real(wp), dimension(:,:), pointer, intent(inout) :: datavals
    class(neutraldata), pointer :: aperptr

    aperptr=>intvars%atmosperturb    ! apparently select case cannot handle a compound statement
    select type (aperptr)
    class is (neutraldata3D_fclaw)
      call intvars%atmosperturb%get_locationsi(flagallpts,zlims,xlims,ylims,zvals,xvals,yvals,datavals)
    class default
      print*, 'WARNING:  attempted to direct feed data (get) to object of wrong type (not neutraldata3D_fclaw)'
      zvals=>null()
      xvals=>null()
      yvals=>null()
      datavals=>null()
    end select
  end subroutine get_locationsi_in


  !> retrieve a pointer that we can directly copy data to
  subroutine get_datainow_ptr_in(intvars,datavals)
    type(gemini_work), intent(inout) :: intvars
    real(wp), dimension(:,:), pointer, intent(inout) :: datavals
    class(neutraldata), pointer :: aperptr

    aperptr=>intvars%atmosperturb    ! apparently select case cannot handle a compound statement
    select type (aperptr)
    class is (neutraldata3D_fclaw)
      datavals=>intvars%atmosperturb%get_datainow_ptr()
    class default
      print*, 'WARNING:  attempted to direct feed data (get) to object of wrong type (not neutraldata3D_fclaw)'
      datavals=>null();
    end select
  end subroutine get_datainow_ptr_in


  !> Notify object that the data have been placed in its buffer and it needs to copy-out
  subroutine set_datainow_in(intvars)
    type(gemini_work), intent(inout) :: intvars
    class(neutraldata), pointer :: aperptr

    aperptr=>intvars%atmosperturb    ! apparently select case cannot handle a compound statement
    select type (aperptr)
    class is (neutraldata3D_fclaw)
      call intvars%atmosperturb%set_datainow()
    class default
      print*, 'WARNING:  attempted to direct feed data (set) to object of wrong type (not neutraldata3D_fclaw)'
    end select
  end subroutine set_datainow_in
end module forestgemini
