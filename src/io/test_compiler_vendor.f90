program test_compiler_vendor

use config, only : get_compiler_vendor

implicit none (type, external)

character(:), allocatable :: vendor

vendor = get_compiler_vendor()

if (len_trim(vendor) == 0) error stop "compiler vendor not determined"

print '(A)',vendor

end program