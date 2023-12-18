module forestgemini_C

use forestgemini, only: interp_file2subgrid_in,grid_from_extents_in,gemini_grid_alloc,checkE1,forceZOH_all, &
            permute_fluidvars, ipermute_fluidvars, tag4refine_location, tag4refine_vperp, clean_param_after_regrid_in, &
            get_locationsi_in,set_datainow_in, get_datainow_ptr_in, swap_statevars,tag4refine_diff, tag4refine_grad, & 
            tag4coarsening_diff


implicit none (type, external)

public

contains
  ! FIXME: obviated needs to get rid of this here and in header file
  !> C wrapper for procedure to compute a grid object given extents and fullgrid reference point.  The class
  !    pointed to by xC must already have been allocated and assigned the correct fortran dynamic type.  
  subroutine grid_from_extents_C(x1lims,x2lims,x3lims,lx1wg,lx2wg,lx3wg,xtype,xC) bind(C,name='grid_from_extents_C')
    real(wp), dimension(2), intent(in) :: x1lims,x2lims,x3lims
    integer(C_INT), intent(in) :: lx1wg,lx2wg,lx3wg
    integer(C_INT), intent(inout) :: xtype
    type(c_ptr), intent(inout) :: xC
    class(curvmesh), pointer :: x

    x=>set_gridpointer_dyntype(xtype,xC)
    call grid_from_extents_in(x1lims,x2lims,x3lims,lx1wg,lx2wg,lx3wg,x)
    ! as an extra step we need to also assign a type to the grid
    xtype=detect_gridtype(x%x1,x%x2,x%x3)
  end subroutine grid_from_extents_C


  !> C wrapper to allocate grid
  subroutine gemini_grid_alloc_C(x1lims,x2lims,x3lims,lx1wg,lx2wg,lx3wg,xtype,xC) bind(C,name='gemini_grid_alloc_C')
    real(wp), dimension(2), intent(in) :: x1lims,x2lims,x3lims
    integer, intent(in) :: lx1wg,lx2wg,lx3wg
    integer, intent(inout) :: xtype
    type(c_ptr), intent(inout) :: xC
    class(curvmesh), pointer :: x

    call gemini_grid_alloc(x1lims,x2lims,x3lims,lx1wg,lx2wg,lx3wg,x,xtype,xC) 
  end subroutine gemini_grid_alloc_C


  !> C wrapper for procedure that reads full data from input file and interpolates to loca worker subgrid
  subroutine interp_file2subgrid_C(cfgC,xtype,xC,fluidvarsC,electrovarsC) bind(C,name='interp_file2subgrid_C')
    type(c_ptr), intent(in) :: cfgC
    integer(C_INT), intent(in) :: xtype
    type(c_ptr), intent(in) :: xC
    type(c_ptr), intent(inout) :: fluidvarsC
    type(c_ptr), intent(inout) :: electrovarsC
    type(gemini_cfg), pointer :: cfg
    class(curvmesh), pointer :: x
    real(wp), dimension(:,:,:,:), pointer :: fluidvars
    real(wp), dimension(:,:,:,:), pointer :: electrovars

    call c_f_pointer(cfgC,cfg)
    x=>set_gridpointer_dyntype(xtype,xC)
    call c_f_pointer(fluidvarsC,fluidvars,[(lx1+4),(lx2+4),(lx3+4),(5*lsp)])
    call c_f_pointer(electrovarsC,electrovars,[(lx1+4),(lx2+4),(lx3+4),7])
    call interp_file2subgrid_in(cfg,x,fluidvars,electrovars)
  end subroutine interp_file2subgrid_C


  !> deal with null cell solutions
  subroutine clean_param_after_regrid_C(iparm,xtype,xC,fluidvarsC,intvarsC) bind(C, name="clean_param_after_regrid_C")
    integer(C_INT), intent(in) :: iparm
    integer(C_INT), intent(in) :: xtype
    type(c_ptr), intent(in) :: xC
    type(c_ptr), intent(inout) :: fluidvarsC
    type(c_ptr), intent(in) :: intvarsC

    class(curvmesh), pointer :: x
    real(wp), dimension(:,:,:,:), pointer :: fluidvars
    type(gemini_work), pointer :: intvars   

    x=>set_gridpointer_dyntype(xtype, xC)
    call c_f_pointer(fluidvarsC,fluidvars,[(lx1+4),(lx2+4),(lx3+4),(5*lsp)])
    call c_f_pointer(intvarsC,intvars)   
    call clean_param_after_regrid_in(iparm,x,fluidvars,intvars)
  end subroutine clean_param_after_regrid_C


  !> echo print variable min/max for checking
  subroutine checkE1_C(fluidvarsC,fluidauxvarsC,electrovarsC,locID) bind(C, name="checkE1_C")
    type(c_ptr), intent(inout) :: fluidvarsC
    type(c_ptr), intent(inout) :: fluidauxvarsC
    type(c_ptr), intent(in) :: electrovarsC
    integer, intent(in) :: locID
    real(wp), dimension(:,:,:,:), pointer :: fluidvars
    real(wp), dimension(:,:,:,:), pointer :: fluidauxvars
    real(wp), dimension(:,:,:,:), pointer :: electrovars

    call c_f_pointer(fluidvarsC,fluidvars,[(lx1+4),(lx2+4),(lx3+4),(5*lsp)])
    call c_f_pointer(fluidauxvarsC,fluidauxvars,[(lx1+4),(lx2+4),(lx3+4),(2*lsp)+9])
    call c_f_pointer(electrovarsC,electrovars,[(lx1+4),(lx2+4),(lx3+4),7])
    call checkE1(fluidvars,fluidauxvars,electrovars,locID)
  end subroutine checkE1_C


  subroutine permute_fluidvars_C(fluidvarsC) bind(C, name='permute_fluidvars_C')
    type(c_ptr), intent(inout) :: fluidvarsC
    real(wp), dimension(:,:,:,:), pointer :: fluidvars

    call c_f_pointer(fluidvarsC,fluidvars,[(lx0+4),(lx2+4),(lx3+4),(5*lsp)])
    call permute_fluidvars(fluidvars)
  end subroutine permute_fluidvars_C

  subroutine ipermute_fluidvars_C(fluidvarsC) bind(C, name='ipermute_fluidvars_C')
    type(c_ptr), intent(inout) :: fluidvarsC
    real(wp), dimension(:,:,:,:), pointer :: fluidvars

    call c_f_pointer(fluidvarsC,fluidvars,[(lx0+4),(lx2+4),(lx3+4),(5*lsp)])
    call ipermute_fluidvars(fluidvars)
  end subroutine ipermute_fluidvars_C


  !> refinement type switchyard
  subroutine tag4refine_C(xtype,xC,fluidvarsC,fluidauxvarsC,electrovarsC,intvarsC,refinetype,flagrefine) &
                            bind(C, name='tag4refine_C')
    integer(C_INT), intent(in) :: xtype
    type(c_ptr), intent(in) :: xC
    type(c_ptr), intent(inout) :: fluidvarsC,fluidauxvarsC,electrovarsC,intvarsC
    integer(C_INT), intent(in) :: refinetype
    logical, intent(inout) :: flagrefine
    class(curvmesh), pointer :: x
    real(wp), dimension(:,:,:,:), pointer :: fluidvars
    real(wp), dimension(:,:,:,:), pointer :: fluidauxvars
    real(wp), dimension(:,:,:,:), pointer :: electrovars
    type(gemini_work), pointer :: intvars

    x=>set_gridpointer_dyntype(xtype, xC)
    call c_f_pointer(fluidvarsC,fluidvars,[(lx1+4),(lx2+4),(lx3+4),(5*lsp)])
    call c_f_pointer(fluidauxvarsC,fluidauxvars,[(lx1+4),(lx2+4),(lx3+4),(2*lsp)+9])
    call c_f_pointer(electrovarsC,electrovars,[(lx1+4),(lx2+4),(lx3+4),7])
    call c_f_pointer(intvarsC,intvars)

    select case (refinetype)
      case (0)    ! do a refinement based on grid location (usually used for static refines)
        call tag4refine_location(x,flagrefine)
      case (1)    ! actively refine based on perpendicular velocity magnitude (e.g. indicative of auroral or neutral forcing)
        call tag4refine_vperp(x,fluidvars,fluidauxvars,electrovars,intvars,flagrefine)
      case (2)
        call tag4refine_diff(x,fluidvars,fluidauxvars,electrovars,intvars,flagrefine)
      case (3)
        call tag4refine_grad(x,fluidvars,fluidauxvars,electrovars,intvars,flagrefine)
      case default
        error stop 'unrecognized refinment criteria from user...'
    end select
  end subroutine tag4refine_C


  !> top-level routine for determining mesh coarsening
  subroutine tag4coarsening_C(xtype,xC,fluidvarsC,fluidauxvarsC,electrovarsC,intvarsC,flagcoarsening) &
                            bind(C, name='tag4coarsening_C')
    integer(C_INT), intent(in) :: xtype
    type(c_ptr), intent(in) :: xC
    type(c_ptr), intent(inout) :: fluidvarsC,fluidauxvarsC,electrovarsC,intvarsC
    logical, intent(inout) :: flagcoarsening
    class(curvmesh), pointer :: x
    real(wp), dimension(:,:,:,:), pointer :: fluidvars
    real(wp), dimension(:,:,:,:), pointer :: fluidauxvars
    real(wp), dimension(:,:,:,:), pointer :: electrovars
    type(gemini_work), pointer :: intvars

    x=>set_gridpointer_dyntype(xtype, xC)
    call c_f_pointer(fluidvarsC,fluidvars,[(lx1+4),(lx2+4),(lx3+4),(5*lsp)])
    call c_f_pointer(fluidauxvarsC,fluidauxvars,[(lx1+4),(lx2+4),(lx3+4),(2*lsp)+9])
    call c_f_pointer(electrovarsC,electrovars,[(lx1+4),(lx2+4),(lx3+4),7])
    call c_f_pointer(intvarsC,intvars)

    call tag4coarsening_diff(x,fluidvars,fluidauxvars,electrovars,intvars,flagcoarsening)
  end subroutine tag4coarsening_C


  !> In case we need to swap state variables
  subroutine swap_statevars_C(fluidvarsC,fluidauxvarsC) bind(C, name='swap_statevars_C')
    type(c_ptr), intent(inout) :: fluidvarsC,fluidauxvarsC
    real(wp), dimension(:,:,:,:), pointer :: fluidvars
    real(wp), dimension(:,:,:,:), pointer :: fluidauxvars

    call c_f_pointer(fluidvarsC,fluidvars,[(lx1+4),(lx2+4),(lx3+4),(5*lsp)])
    call c_f_pointer(fluidauxvarsC,fluidauxvars,[(lx1+4),(lx2+4),(lx3+4),(2*lsp)+9])

    call swap_statevars(fluidvars,fluidauxvars)
  end subroutine swap_statevars_C


  !> get c pointers to magnetic coordinates of mesh sites (cell centers), note that these include ghost cells
  !  Due to the way the vtu files are created in forestclaw we likely need interface values for the cell locations
  !    As such this particular routine is likely not useful and we need separate facilities to computer and return
  !    the cell edge values.  
  subroutine get_grid_magcoords_C(xtype,xC,mlonC,mlatC,altC) bind(C, name='get_grid_magcoords_C')
    integer(C_INT), intent(in) :: xtype
    type(c_ptr), intent(in) :: xC
    type(c_ptr), intent(inout) :: mlonC,mlatC,altC
    class(curvmesh), pointer :: x

    x=>set_gridpointer_dyntype(xtype, xC)

    ! for now just use geographic since already stored; also note that these do include ghost cells!
    mlonC=c_loc(x%glon)
    mlatC=c_loc(x%glat)
    altC=c_loc(x%alt)    
  end subroutine get_grid_magcoords_C


  !> force computation (if needed) of interface cell locations (geographic coordinates) and return a C pointer
  subroutine get_grid_magcoordsi_C(xtype,xC,mloniC,mlatiC,altiC) bind(C, name='get_grid_magcoordsi_C')
    integer(C_INT), intent(in) :: xtype
    type(c_ptr), intent(in) :: xC
    type(c_ptr), intent(inout) :: mloniC,mlatiC,altiC
    class(curvmesh), pointer :: x

    x=>set_gridpointer_dyntype(xtype, xC)

    ! call for generation of interface locations if not already done
    call x%calc_geographici()

    ! for now just use geographic since already stored; also note that these do include ghost cells!
    mloniC=c_loc(x%gloni)
    mlatiC=c_loc(x%glati)
    altiC=c_loc(x%alti)    
  end subroutine get_grid_magcoordsi_C


  !> grab pointers and size of interpolation target locations
  subroutine get_locationsi_C(intvarsC,flagallpts,zlims,xlims,ylims, &
                                zvalsC,xvalsC,yvalsC,datavalsC,lpts,lparms) bind(C,name='get_locationsi_C')
    type(C_PTR), intent(inout) :: intvarsC
    logical, intent(in) :: flagallpts
    real(wp), dimension(2), intent(in) :: zlims,xlims,ylims   
    type(C_PTR), intent(inout) :: zvalsC,xvalsC,yvalsC,datavalsC
    integer(C_INT), intent(inout) :: lpts,lparms
    real(wp), dimension(:), pointer :: zvals,xvals,yvals
    real(wp), dimension(:,:), pointer :: datavals
    type(gemini_work), pointer :: intvars

    call c_f_pointer(intvarsC,intvars)
    call get_locationsi_in(intvars,flagallpts,zlims,xlims,ylims,zvals,xvals,yvals,datavals)
    zvalsC=c_loc(zvals)
    xvalsC=c_loc(xvals)
    yvalsC=c_loc(yvals)
    datavalsC=c_loc(datavals)
    lpts=size(datavals,1)
    lparms=size(datavals,2)
    !print*, 'Internal size chk:  ',lpts,lparms,shape(datavals)
  end subroutine get_locationsi_C


  !> retrieve (again?) the pointer to the data storage location for direct-feed
  subroutine get_datainow_ptr_C(intvarsC,datavalsC) bind(C,name='get_datainow_ptr_C')
    type(C_PTR), intent(inout) :: intvarsC
    type(C_PTR), intent(inout) :: datavalsC
    real(wp), dimension(:,:), pointer :: datavals
    type(gemini_work), pointer :: intvars  

    call c_f_pointer(intvarsC,intvars)   
    call get_datainow_ptr_in(intvars,datavals)
    datavalsC=c_loc(datavals)
  end subroutine get_datainow_ptr_C


  !> inform gemini object that the data are now in place and can be rotated/copied out
  subroutine set_datainow_C(intvarsC) bind(C,name='set_datainow_C')
    type(C_PTR), intent(inout) :: intvarsC
    type(gemini_work), pointer :: intvars

    call c_f_pointer(intvarsC,intvars)          
    call set_datainow_in(intvars)
  end subroutine set_datainow_C

  !> test to force ghost cells to zero-order hold extrapolated value
  subroutine forceZOH_all_C(fluidvarsC) bind(C,name='forceZOH_all_C')
    type(c_ptr), intent(inout) :: fluidvarsC
    real(wp), dimension(:,:,:,:), pointer :: fluidvars

    call c_f_pointer(fluidvarsC,fluidvars,[(lx1+4),(lx2+4),(lx3+4),(5*lsp)])
    call forceZOH_all(fluidvars)
  end subroutine forceZOH_all_C
end module forestgemini_C
