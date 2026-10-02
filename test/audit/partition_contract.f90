program partition_contract
use autogrid, only: max_mpi, grid_auto
implicit none
integer :: n2,n3,budget,p2,p3,expected,actual,lid2,lid3
do n2=1,33
  do n3=1,33
    do budget=1,32
      expected=1
      do p2=1,n2
        if (mod(n2,p2)/=0 .or. (n2>1 .and. n2/p2<2)) cycle
        do p3=1,n3
          if (mod(n3,p3)/=0 .or. (n3>1 .and. n3/p3<2)) cycle
          if (p2*p3<=budget) expected=max(expected,p2*p3)
        enddo
      enddo
      actual=max_mpi(n2,n3,budget)
      if (actual/=expected) then
        print *, n2,n3,budget,expected,actual
        error stop "max_mpi disagrees with exhaustive partition oracle"
      endif
      call grid_auto(n2,n3,actual,lid2,lid3)
      if (lid2*lid3/=actual .or. mod(n2,lid2)/=0 .or. mod(n3,lid3)/=0) error stop "invalid process grid"
      if (n2>1 .and. n2/lid2<2) error stop "singular local x2"
      if (n3>1 .and. n3/lid3<2) error stop "singular local x3"
    enddo
  enddo
enddo
print *, "34848 odd/even/prime/singleton/single-CPU partitions passed"
end program
