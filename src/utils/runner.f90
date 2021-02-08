module runner

use phys_consts, only : wp
use config, only : read_configfile, gemini_cfg
use timeutils, only : date_filename,dateinc
use pathlib, only : get_suffix

implicit none (type, external)

contains

subroutine clean_output(path)

character(*), intent(in) :: path

type(gemini_cfg) :: cfg
integer, dimension(3) :: ymd
real(wp) :: UTsec
character(:), allocatable :: fn, suffix
logical :: exists

cfg%outdir = path
cfg%infile = path // '/inputs/config.nml'
inquire(file=cfg%infile, exist=exists)
if(.not.exists) error stop 'gemini3d.run: not a file: ' // cfg%infile

call read_configfile(cfg)

ymd = cfg%ymd0
UTsec = cfg%UTsec0

suffix = get_suffix(cfg%indatsize)
fn = date_filename(cfg%outdir, ymd, UTsec) // suffix

do
  !! new filename, add the 1 if it is the first
  fn = date_filename(cfg%outdir, ymd, UTsec) // suffix

  inquire(file=fn, exist=exists)
  if ( .not. exists ) exit
  !! last output file
  print *, 'delete: ', fn
  call unlink(fn)

  !! next time
  call dateinc(cfg%dtout, ymd,UTsec)
end do

end subroutine clean_output


subroutine unlink(path)
character(*), intent(in) :: path
integer :: i
logical :: e

inquire(file=path, exist=e)
if (.not.e) return

open(newunit=i, file=path, status='old')
close(i, status='delete')
end subroutine unlink

end module runner
