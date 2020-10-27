submodule (io) milestone

use timeutils, only: date_filename,dateinc
use h5fortran, only: h5exist

implicit none (type,external)   !! external procedures must be explicitly denoted thusly

contains

module procedure find_milestone

!> search path having output rate cadence (s) and find the last file that is a milestone.
integer, dimension(3) :: ymd
real(wp) :: UTsec
character(:), allocatable :: fn, suffix
logical :: exists
integer :: it,lfn
real(wp) :: tsim

ymd = cfg%ymd0
UTsec = cfg%UTsec0
ymdmile = cfg%ymd0
UTsecmile = cfg%UTsec0
filemile = date_filename(cfg%outdir, ymd, UTsec)   !! This presumes the first file output is a milestone, note we don't actually test the situation wheere a first output was not produced...  User should not be restarting in that case...

it=1
tsim = 0
tmile = 0

suffix = get_suffix(cfg%indatsize)

if (suffix /= '.h5') return

milesearch : do
  !! new filename, add the 1 if it is the first
  fn = date_filename(cfg%outdir, ymd,UTsec)
  if (it==1) then
    lfn=len(fn)
    fn(lfn:lfn)='1'
  end if
  fn=fn // suffix

  inquire(file=fn, exist=exists)
  if ( .not. exists ) exit milesearch
  !! last output file

  if (h5exist(fn, '/nsall')) then
    !! this file is milestone
    ymdmile=ymd
    UTsecmile=UTsec
    filemile=fn
    tmile=tsim
  end if

  !! next time
  call dateinc(cfg%dtout, ymd,UTsec)
  it=it+1
  tsim = tsim + cfg%dtout
end do milesearch

end procedure find_milestone

!! need to add ishdf5 to


end submodule milestone
