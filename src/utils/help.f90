module help

use, intrinsic :: iso_fortran_env, only : compiler_version
use config, only : get_compiler_vendor

implicit none (type, external)
private
public :: help_gemini_bin, help_gemini_run, help_magcalc_bin

contains

subroutine help_gemini_bin(git_revision)

character(*), intent(in) :: git_revision

print '(/,A,/)', 'GEMINI-3D: gemini.bin ' // git_revision
print '(A)', 'by Matthew Zettergren'
print '(A)', 'GLOW and auroral interfaces by Guy Grubbs'
print '(A)', 'build system and software engineering by Michael Hirsch'
print '(A)', 'Compiler vendor: '// get_compiler_vendor()
print '(A)', 'Compiler version: ' // compiler_version()
print '(/,A)', 'must specify simulation output directory. Example:'
print '(/,A,/)', '  mpiexec -np 4 build/gemini.bin /path/to/simulation_outputs'
print '(A)', '-dryrun    allows quick check of first time step'
print '(A)', '-manual_grid lx2 lx3    defines the number of MPI processes along x2 and x3.'
print '(A)', '  If -manual_grid is not specified, the MPI processes are auto-assigned along x2 and x3.'
stop 'EOF: gemini.bin'

end subroutine help_gemini_bin


subroutine help_gemini_run(git_revision)

character(*), intent(in) :: git_revision

print '(/,A,/)', 'GEMINI-3D: gemini3d.run ' // git_revision
print '(A)', 'Compiler vendor: '// get_compiler_vendor()
print '(A)', 'Compiler version: ' // compiler_version()
print '(/,A)', 'must specify simulation output directory. Example:'
print '(/,A,/)', '  build/gemini3d.run /path/to/simulation_outputs'
print '(A)', '-dryrun    allows quick check of first time step'
print '(A)', '-n   manually specify number of MPI images (default auto-calculate)'
print '(A)', '-gemexe   path to gemini.bin'
print '(A)', '-mpiexec   path to mpiexec'
stop 'EOF: gemini3d.run'

end subroutine help_gemini_run


subroutine help_magcalc_bin(git_revision)
character(*), intent(in) :: git_revision

print '(/,A,/)', 'GEMINI-3D: magcalc.bin ' // git_revision
print '(A)', 'Compiler vendor: '// get_compiler_vendor()
print '(A)', 'Compiler version: ' // compiler_version()
print '(/,A)', 'must specify input directory and fieldpoint file. Example:'
print '(/,A,/)', 'mpiexec -n 4 build/magcalc.bin test2d_fang test2d_fang/fieldpoint'
print '(A)', '-dryrun option allows quick check of first time step'
stop 'EOF: magcalc.bin'
end subroutine help_magcalc_bin


end module help
