program ex1

use filesystem

implicit none

character(:), allocatable :: h

print '(a)', "current working dir " // get_cwd()

h = get_homedir()
if(len_trim(h) == 0) error stop "homedir failed"
print '(a)', "home dir " // h

h = expanduser('~')
if(len_trim(h) == 0) error stop "expanduser failed"
print '(a)', "expanduser('~') " // h

deallocate(h)

end program
