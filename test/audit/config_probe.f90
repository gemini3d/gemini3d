! Audit addition 2026-09-16. Apache-2.0.
program audit_config
use phys_consts, only: wp,mindens
use gemini3d_config, only: gemini_cfg,read_configfile
implicit none
type(gemini_cfg) :: cfg
character(2048) :: file
integer :: i
cfg%outdir='.'
do i=1,command_argument_count()
 call get_command_argument(i,file)
 cfg%infile=trim(file)
 call read_configfile(cfg)
 print *, "AUDIT_CONFIG",cfg%msis_version,cfg%dtneuBG,cfg%flagevibcool,cfg%flagperiodic,mindens
enddo
end program
