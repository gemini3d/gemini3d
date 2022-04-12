!> these are routines that are called by both the C and fortran frontends; they access module variables
!    that do not differ patch-to-patch like config variables or number of species.  
module libgemini_common

contains
  !> basic command line and grid size determination; although this takes in the split parameters it doesn't actually do any mpi
  subroutine cli_config_gridsize_C(p,lid2in,lid3in) bind(C, name="cli_config_grid_size_C")
    type(c_params), intent(in) :: p
    integer, intent(inout) :: lid2in,lid3in

    call cli_config_gridsize(p,lid2in,lid3in)
  end subroutine cli_config_gridsize


  !> return number of species from phys_consts module; this is a global variable across the entire
  !   simulation so even if multiple patches are used this will still be valid.  
  subroutine get_species_size_C(lspout) bind(C, name="get_species_size_C")
    integer, intent(inout) :: lspout

    lspout=lsp
  end subroutine get_species_size_C


  ! FIXME:  this should be in a common library, i.e. used by both fortran ane C
  !> return some data from cfg that is needed in the main program
  subroutine get_config_vars_C(flagneuBG,flagdneu,dtneuBG,dtneu) bind(C, name="get_config_vars_C")
    logical, intent(inout) :: flagneuBG
    integer, intent(inout) :: flagdneu
    real(wp), intent(inout) :: dtneuBG,dtneu

    flagneuBG=cfg%flagneuBG
    flagdneu=cfg%flagdneu
    dtneuBG=cfg%dtneuBG
    dtneu=cfg%dtneu
  end subroutine get_config_vars_C
end module libgemini_common
