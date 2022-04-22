!> This module contains C/CXX wrappers for functions in libgemini_mpi.  These routines match those in libgemini_mpi.f90 and are
!!   principally meant to convert the C pointers to various data objects into fortran pointers (including in the case of the
!!   grid a class pointer (pointer to polymorphic object).  Other polymorhpic objects (neutraldata, etc.) are kept in a static
!!   derived type (intvars::gemini_work) and don't need to be passes as class pointers.  
module gemini3d_mpi_C

contains

! FIXME: put code here :)

end module gemini3d_mpi_C
