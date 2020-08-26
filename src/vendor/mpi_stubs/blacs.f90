SUBROUTINE blacs_gridinit( CNTXT, C, NPROW, NPCOL )
IMPLICIT NONE
INTEGER CNTXT, NPROW, NPCOL
CHARACTER C
error stop 'Error. BLACS_GRIDINIT should not be called.'
END SUBROUTINE blacs_gridinit

SUBROUTINE blacs_gridinfo( CNTXT, NPROW, NPCOL, MYROW, MYCOL )
IMPLICIT NONE
INTEGER CNTXT, NPROW, NPCOL, MYROW, MYCOL
error stop 'Error. BLACS_GRIDINFO should not be called.'
END SUBROUTINE blacs_gridinfo

SUBROUTINE blacs_gridexit( CNTXT )
IMPLICIT NONE
INTEGER CNTXT
error stop 'Error. BLACS_GRIDEXIT should not be called.'
END SUBROUTINE blacs_gridexit
