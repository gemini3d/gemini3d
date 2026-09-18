character(:), allocatable :: r
character(*), intent(in) :: path
logical, intent(in), optional :: strict

character(kind=c_char, len=:), allocatable :: cbuf
integer(C_SIZE_T) :: N
logical(c_bool) :: s

allocate(character(max_path()) :: cbuf)
s = .false.
if(present(strict)) s = strict

N = len(cbuf)
