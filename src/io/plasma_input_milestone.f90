submodule (io:plasma_input) plasma_input_milestone

use timeutils, only: date_filename,dateinc
use h5fortran, only: hdf5_file 

implicit none (type,external)   !! external procedures must be explicitly denoted thusly

contains

module procedure find_lastfile

!> Compute the last file before the target date based on a start date and cadence.  The
!  file is assumed to exist and programmer needs to check for existence outside this 
!  procedure.

integer, dimension(3) :: ymdnext
real(wp) :: UTsecnext
logical :: flagend


ymd=ymd0
UTsec=UTsec0
ymdnext=ymd0
UTsecnext=UTsec0
flagend=ymdnext(1)>=ymdtarget(1) .and. ymdnext(2)>=ymdtarget(2) .and. ymdnext(3)>=ymdtarget(3) .and. UTsecnext>=UTsectarget  ! in case the first time step is the last before target
do while ( .not.(flagend) )
  ymd=ymdnext
  UTsec=UTsecnext
  call dateinc(cadence,ymdnext,UTsecnext)
  flagend=ymdnext(1)>=ymdtarget(1) .and. ymdnext(2)>=ymdtarget(2) .and. ymdnext(3)>=ymdtarget(3) .and. UTsecnext>=UTsectarget
end do
! When the loops exits ymd,UTsec have the date of the last output file before the given target time OR the first output file in the case that the target date is before the begin date...

end procedure find_lastfile


module procedure find_milestone

!> search path having output rate cadence (s) and find the
!  last file that is a milestone.

integer, dimension(3) :: ymd
real(wp) :: UTsec
character(:), allocatable :: fn
type(hdf5_file) :: h5f
logical flagexists,flagend,flagmilestone


ymd=ymd0
UTsec=UTsec0
ymdmile=ymd0
UTsecmile=UTsec0
filemile=date_filename(path,ymd0,UTsec0)   !! This presumes the file output is a milestone
flagend=.false.
do while ( .not.(flagend) )
  !! new filename
  fn=date_filename(path,ymd,UTsec)
  fn=fn // suffix

  !! is the file there
  inquire(file=fn, exist=flagexists)
  if ( flagexists ) then
    call h5f%initialize(fn, status='old', action='r')    ! init hdf5 object instance
    flagmilestone=h5f%exists('/nsall')
    call h5f%finalize()    ! destruct hdf5 object instance so it can be reused
    if (flagmilestone) then
      ymdmile=ymd
      UTsecmile=UTsec
      filemile=fn 
    end if
  else
    flagend=.true.    !we've hit the last output file
  end if
end do

end procedure find_milestone

end submodule plasma_input_milestone
