module phys_consts

use, intrinsic:: iso_fortran_env, only: wp=>real@realbits@

implicit none (type, external)
public

integer, parameter :: comp_lvl = 3
!! compression level when writing files.
!! 0 disables compression, 1->9 increasing compression
!! only 2-D or higher variables are compressed

!> colored output text (for Unix-like systems at least).
!> It should be compatible across compiler vendors.
character(5), parameter :: &
  red = char(27) // '[31m', &
  black = char(27) // '[0m'

!MATHEMATICAL CONSTANTS
real(wp), parameter :: pi = 4*atan(1.0_wp)


!PHYSICAL CONSTANTS
real(wp), parameter :: kB=1.38064852e-23_wp, &      !Boltzmann constant [J K^-1] = [m^2 kg s^-2 K^-1]
                       elchrg=1.60217662e-19_wp, &  !elementary charge
                       amu=1.660539040e-27_wp, &    !atomic mass unit
                       Gconst=6.67408e-11_wp, &     !Universal gravitation constant
                       mu0=4*pi*1e-7_wp        !vacuum permeability


!> EARTH-RELATED PARAMETERS
real(wp), parameter :: Mmag=7.94e22_wp
  !! Earth's magnetic moment
real(wp), parameter :: Re = 6371.0e3_wp
  !! Earth Radius [meters]
real(wp), parameter :: Me = 5.9722e24_wp
  !! Earth mass


!> ION DATA (need to be doubled?)
integer, parameter :: lsp=7
  !! number of ion/electron species
real(wp), parameter :: ms(lsp)=[real(wp) :: 16, 30, 28, 32, 14, 1, 5.485799090e-4_wp]*amu
  !! mass of each species
real(wp), parameter :: qs(lsp)=[real(wp) :: 1,1,1,1,1,1,-1]*elchrg
  !! charge of each species
real(wp), parameter :: gammas(lsp)=[5.0_wp/3, &
                                    7.0_wp/5, &
                                    7.0_wp/5, &
                                    7.0_wp/5, &
                                    5.0_wp/3, &
                                    5.0_wp/3, &
                                    5.0_wp/3]
  !! adiabatic index for each speces


!> NEUTRAL DATA
integer, parameter :: ln=4, lnchem=6
  !! number of neutral densities, and number of neutrals in chem. rxns.
real(wp),parameter :: mn(ln)=[real(wp) :: 16,28,32,1]*amu
  !! mass of neutral species

!AURORAL DATA
integer, parameter :: lwave=15 !spectral auroral lines tracked in GLOWv0.982
real(wp), parameter :: wavelengths(lwave)=[real(wp) :: 3371, &
                                        4278, &
                                        5200, &
                                        5577, &
                                        6300, &
                                        7320, &
                                        10400, &
                                        3644, &
                                        7774, &
                                        8446, &
                                        3276, &
                                        1700, & !this is N2(A) LBH
                                        1356, &
                                        1493, &
                                        1304] !wavelength of each auroral line, housekeeping

!> HOUSEKEEPING PARAMETERS for conditioning densities
real(wp), parameter :: mindens     = 1.0e-100_wp
real(wp), parameter :: mindensnull = 1.0e-20_wp
real(wp), parameter :: mindensdiv  = 1.0e-5_wp

!To control the amount of console output; can be changed by user command line flag "-debug"
logical :: debug=.false.

end module phys_consts
