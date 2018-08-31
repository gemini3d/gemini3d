module phys_consts
use, intrinsic:: iso_fortran_env,only: wp=>real64

implicit none

public

!MATHEMATICAL CONSTANTS
real(wp), parameter :: pi = 4.0_wp*atan(1.0_wp)


!PHYSICAL CONSTANTS
real(wp), parameter :: kB=1.38064852e-23_wp, &      !Boltzmann constant
                       elchrg=1.60217662e-19_wp, &  !elementary charge
                       amu=1.660539040e-27_wp, &    !atomic mass unit
                       Gconst=6.67408e-11_wp, &     !Universal gravitation constant
                       mu0=4.0_wp*pi*1e-7_wp        !vacuum permeability


!EARTH-RELATED PARAMETERS
real(wp), parameter :: Mmag=7.94e22_wp             !Earth's magnetic moment
real(wp), parameter :: Re = 6371.0e3_wp            !Earth Radius [meters]
real(wp), parameter :: Me = 5.9722e24_wp          !mass


!ION DATA (need to be doubled?)
integer, parameter :: lsp=7                                                  !number of ion/electron species
real(wp), parameter :: ms(lsp)=[16.0_wp, 30.0_wp, 28.0_wp, 32.0_wp, 14.0_wp, 1.0_wp, &
                                          5.485799090e-4_wp]*amu                !mass of each species
real(wp), parameter :: qs(lsp)=[1.0_wp,1.0_wp,1.0_wp,1.0_wp,1.0_wp,1.0_wp,-1.0_wp]*elchrg    !charge of each species
real(wp), parameter :: gammas(lsp)=[5.0_wp/3.0_wp, &
                                    7.0_wp/5.0_wp, &
                                    7.0_wp/5.0_wp, &
                                    7.0_wp/5.0_wp, &
                                    5.0_wp/3.0_wp, &
                                    5.0_wp/3.0_wp, &
                                    5.0_wp/3.0_wp]          !adiabatic index for each speces


!NEUTRAL DATA
integer, parameter :: ln=4, lnchem=6    !number of neutral densities, and number of neutrals in chem. rxns.
real(wp),parameter :: mn(ln)=[16.0_wp,28.0_wp,32.0_wp,1.0_wp]*amu      !mass of neutral species


!HOUSEKEEPING PARAMETERS
real(wp), parameter :: mindens     = 1.0e-100_wp
real(wp), parameter :: mindensnull = 1.0e-20_wp
real(wp), parameter :: mindensdiv  = 1.0e-5_wp

logical :: debug=.false.

contains

end module phys_consts
