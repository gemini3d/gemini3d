! Copyright 2021 Matthew Zettergren

! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!   http://www.apache.org/licenses/LICENSE-2.0

! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

Program Gemini3D_main
!! MAIN PROGRAM FOR GEMINI3D

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use, intrinsic :: iso_c_binding, only : c_bool, c_int, c_float, c_char, c_null_char

use gemini3d, only : gemini_main
use config, only : gemini_cfg
use mpimod, only : mpi_cfg, mpibreakdown
use mpi, only : mpi_init

implicit none (type, external)

type, bind(C) :: c_params
!! this MUST match gemini3d.h and libgemini.f90 exactly including order
logical(c_bool) :: fortran_cli, debug, dryrun
character(kind=c_char) :: out_dir(1000)
!! .ini [base]
integer(c_int) :: ymd(3)
real(kind=c_float) :: UTsec0, tdur, dtout, activ(3), tcfl, Teinf
!! .ini
end type c_params


integer(c_int) :: lid2in, lid3in

character(8) :: date
character(10) :: time

integer :: ierr

type(c_params) :: p

!> INITIALIZE MESSING PASSING VARIABLES, IDS ETC.
call mpi_init(ierr)
if (ierr/=0) error stop 'gemini.bin: failed mpi_init'

p%fortran_cli = .true.
p%out_dir(1) = c_null_char

lid2in = -1
lid3in = -1
!! out_dir, lid2in, lid3in, are ignored when fortran_cli=.true.
call gemini_main(p, lid2in, lid3in)

!> SHUT DOWN MPI
ierr = mpibreakdown()

if (ierr /= 0) then
  write(stderr, *) 'GEMINI: abnormal MPI shutdown code', ierr, 'Process #', mpi_cfg%myid,' /',mpi_cfg%lid-1
  error stop
endif

call date_and_time(date,time)
print '(/,A,I0,A,I0,A)', 'GEMINI normal termination, Process # ', mpi_cfg%myid,' / ',mpi_cfg%lid-1, ' at ' // date // 'T' // time

end program
