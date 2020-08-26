SUBROUTINE DESCINIT( DESC, M, N, MB, NB, IRSRC, ICSRC, ICTXT, LLD, INFO )
IMPLICIT NONE
INTEGER ICSRC, ICTXT, INFO, IRSRC, LLD, M, MB, N, NB
INTEGER DESC( * )
error stop 'DESCINIT should not be called.'
END SUBROUTINE DESCINIT

INTEGER FUNCTION numroc( N, NB, IPROC, ISRCPROC, NPROCS )
IMPLICIT NONE
INTEGER N, NB, IPROC, ISRCPROC, NPROCS
IF ( NPROCS .ne. 1 ) error stop 'Last parameter from NUMROC should be 1'
IF ( IPROC .ne. 0 ) error stop 'IPROC should be 0 in NUMROC.'
NUMROC = N
END FUNCTION numroc

SUBROUTINE pcpotrf( UPLO, N, A, IA, JA, DESCA, INFO )
IMPLICIT NONE
CHARACTER          UPLO
INTEGER            IA, INFO, JA, N
INTEGER            DESCA( * )
COMPLEX            A( * )
error stop 'PCPOTRF should not be called.'
END SUBROUTINE pcpotrf

SUBROUTINE pcgetrf( M, N, A, IA, JA, DESCA, IPIV, INFO )
IMPLICIT NONE
INTEGER            IA, INFO, JA, M, N
INTEGER            DESCA( * ), IPIV( * )
COMPLEX            A( * )
error stop 'PCGETRF should not be called.'
END SUBROUTINE pcgetrf

SUBROUTINE pctrtrs( UPLO, TRANS, DIAG, N, NRHS, A, IA, JA, DESCA, B, IB, JB, DESCB, INFO )
IMPLICIT NONE
CHARACTER          DIAG, TRANS, UPLO
INTEGER            IA, IB, INFO, JA, JB, N, NRHS
INTEGER            DESCA( * ), DESCB( * )
COMPLEX            A( * ), B( * )
error stop 'PCTRTRS should not be called.'
END SUBROUTINE pctrtrs

SUBROUTINE pzpotrf( UPLO, N, A, IA, JA, DESCA, INFO )
IMPLICIT NONE
CHARACTER          UPLO
INTEGER            IA, INFO, JA, N
INTEGER            DESCA( * )
COMPLEX(kind=kind(0.0D0)) ::     A( * )
error stop 'PZPOTRF should not be called.'
END SUBROUTINE pzpotrf

SUBROUTINE pzgetrf( M, N, A, IA, JA, DESCA, IPIV, INFO )
IMPLICIT NONE
INTEGER            IA, INFO, JA, M, N
INTEGER            DESCA( * ), IPIV( * )
COMPLEX(kind=kind(0.0D0)) ::     A( * )
error stop 'PZGETRF should not be called.'
END SUBROUTINE pzgetrf

SUBROUTINE pztrtrs( UPLO, TRANS, DIAG, N, NRHS, A, IA, JA, DESCA, B, IB, JB, DESCB, INFO )
IMPLICIT NONE
CHARACTER          DIAG, TRANS, UPLO
INTEGER            IA, IB, INFO, JA, JB, N, NRHS
INTEGER            DESCA( * ), DESCB( * )
COMPLEX(kind=kind(0.0D0)) ::     A( * ), B( * )
error stop 'PZTRTRS should not be called.'
END SUBROUTINE pztrtrs

SUBROUTINE pspotrf( UPLO, N, A, IA, JA, DESCA, INFO )
IMPLICIT NONE
CHARACTER          UPLO
INTEGER            IA, INFO, JA, N
INTEGER            DESCA( * )
REAL               A( * )
error stop 'PSPOTRF should not be called.'
END SUBROUTINE pspotrf

SUBROUTINE psgetrf( M, N, A, IA, JA, DESCA, IPIV, INFO )
IMPLICIT NONE
INTEGER            IA, INFO, JA, M, N
INTEGER            DESCA( * ), IPIV( * )
REAL               A( * )
error stop 'PSGETRF should not be called.'
END SUBROUTINE psgetrf

SUBROUTINE pstrtrs( UPLO, TRANS, DIAG, N, NRHS, A, IA, JA, DESCA, B, IB, JB, DESCB, INFO )
IMPLICIT NONE
CHARACTER          DIAG, TRANS, UPLO
INTEGER            IA, IB, INFO, JA, JB, N, NRHS
INTEGER            DESCA( * ), DESCB( * )
REAL               A( * ), B( * )
error stop 'PSTRTRS should not be called.'
END SUBROUTINE pstrtrs

SUBROUTINE pdpotrf( UPLO, N, A, IA, JA, DESCA, INFO )
IMPLICIT NONE
CHARACTER          UPLO
INTEGER            IA, INFO, JA, N
INTEGER            DESCA( * )
DOUBLE PRECISION   A( * )
error stop 'PDPOTRF should not be called.'
END SUBROUTINE pdpotrf

SUBROUTINE pdgetrf( M, N, A, IA, JA, DESCA, IPIV, INFO )
IMPLICIT NONE
INTEGER            IA, INFO, JA, M, N
INTEGER            DESCA( * ), IPIV( * )
DOUBLE PRECISION   A( * )
error stop 'PDGETRF should not be called.'
END SUBROUTINE pdgetrf

SUBROUTINE pdtrtrs( UPLO, TRANS, DIAG, N, NRHS, A, IA, JA, DESCA, B, IB, JB, DESCB, INFO )
IMPLICIT NONE
CHARACTER          DIAG, TRANS, UPLO
INTEGER            IA, IB, INFO, JA, JB, N, NRHS
INTEGER            DESCA( * ), DESCB( * )
DOUBLE PRECISION   A( * ), B( * )
error stop 'PDTRTRS should not be called.'
END SUBROUTINE pdtrtrs

SUBROUTINE INFOG2L( GRINDX, GCINDX, DESC, NPROW, NPCOL, MYROW, MYCOL, LRINDX, LCINDX, RSRC, CSRC )
IMPLICIT NONE
INTEGER            CSRC, GCINDX, GRINDX, LRINDX, LCINDX, MYCOL, MYROW, NPCOL, NPROW, RSRC
INTEGER            DESC( * )
error stop 'INFOG2L should not be called.'
END SUBROUTINE INFOG2L

INTEGER FUNCTION INDXG2P( INDXGLOB, NB, IPROC, ISRCPROC, NPROCS )
IMPLICIT NONE
INTEGER            INDXGLOB, IPROC, ISRCPROC, NB, NPROCS
error stop 'INFOG2L should not be called.'
END FUNCTION INDXG2P

SUBROUTINE pcscal(N, ALPHA, X, IX, JX, DESCX, INCX)
IMPLICIT NONE
INTEGER            INCX, N, IX, JX
COMPLEX            ALPHA
COMPLEX            X( * )
INTEGER            DESCX( * )
error stop 'PCSCAL should not be called.'
END SUBROUTINE pcscal

SUBROUTINE pzscal(N, ALPHA, X, IX, JX, DESCX, INCX)
IMPLICIT NONE
INTEGER            INCX, N, IX, JX
COMPLEX(kind=kind(0.0D0)) :: ALPHA, X( * )
INTEGER            DESCX( * )
error stop 'PZSCAL should not be called.'
END SUBROUTINE pzscal

SUBROUTINE pdscal(N, ALPHA, X, IX, JX, DESCX, INCX)
IMPLICIT NONE
INTEGER            INCX, N, IX, JX
DOUBLE PRECISION   ALPHA
DOUBLE PRECISION   X( * )
INTEGER            DESCX( * )
error stop 'PDSCAL should not be called.'
END SUBROUTINE pdscal

SUBROUTINE psscal(N, ALPHA, X, IX, JX, DESCX, INCX)
IMPLICIT NONE
INTEGER            INCX, N, IX, JX
REAL               ALPHA
REAL               X( * )
INTEGER            DESCX( * )
error stop 'PSSCAL should not be called.'
END SUBROUTINE psscal

SUBROUTINE pzdot ( N, DOT, X, IX, JX, DESCX, INCX, Y, IY, JY, DESCY, INCY )
IMPLICIT NONE
INTEGER N, IX, JX, IY, JY, INCX, INCY
INTEGER DESCX(*), DESCY(*)
COMPLEX(kind=kind(0.0D0)) :: X(*), Y(*)
DOUBLE PRECISION DOT
error stop 'PZDOT should not be called.'
END SUBROUTINE pzdot

SUBROUTINE pcdot ( N, DOT, X, IX, JX, DESCX, INCX, Y, IY, JY, DESCY, INCY )
IMPLICIT NONE
INTEGER N, IX, JX, IY, JY, INCX, INCY
INTEGER DESCX(*), DESCY(*)
COMPLEX X(*), Y(*)
REAL DOT
error stop 'PCDOT should not be called.'
END SUBROUTINE pcdot

SUBROUTINE pddot ( N, DOT, X, IX, JX, DESCX, INCX, Y, IY, JY, DESCY, INCY )
IMPLICIT NONE
INTEGER N, IX, JX, IY, JY, INCX, INCY
INTEGER DESCX(*), DESCY(*)
DOUBLE PRECISION X(*), Y(*), DOT
error stop 'PDDOT should not be called.'
END SUBROUTINE pddot

SUBROUTINE psdot ( N, DOT, X, IX, JX, DESCX, INCX, Y, IY, JY, DESCY, INCY )
IMPLICIT NONE
INTEGER N, IX, JX, IY, JY, INCX, INCY
INTEGER DESCX(*), DESCY(*)
REAL X(*), Y(*), DOT
error stop 'PSDOT should not be called.'
END SUBROUTINE psdot

SUBROUTINE zgebs2d( CONTXT, SCOPE, TOP, M, N, A, LDA )
IMPLICIT NONE
INTEGER CONTXT, M, N, LDA
COMPLEX(kind=kind(0.0D0)) :: A(*)
CHARACTER SCOPE, TOP
error stop 'ZGEBS2D should not be called.'
END SUBROUTINE zgebs2d

SUBROUTINE cgebs2d( CONTXT, SCOPE, TOP, M, N, A, LDA )
IMPLICIT NONE
INTEGER CONTXT, M, N, LDA
COMPLEX A(*)
CHARACTER SCOPE, TOP
error stop 'CGEBS2D should not be called.'
END SUBROUTINE cgebs2d

SUBROUTINE sgebs2d( CONTXT, SCOPE, TOP, M, N, A, LDA )
IMPLICIT NONE
INTEGER CONTXT, M, N, LDA
REAL A(*)
CHARACTER SCOPE, TOP
error stop 'SGEBS2D should not be called.'
END SUBROUTINE sgebs2d

SUBROUTINE dgebs2d( CONTXT, SCOPE, TOP, M, N, A, LDA )
IMPLICIT NONE
INTEGER CONTXT, M, N, LDA
DOUBLE PRECISION A(*)
CHARACTER SCOPE, TOP
error stop 'DGEBS2D should not be called.'
END SUBROUTINE dgebs2d

SUBROUTINE zgebr2d( CONTXT, SCOPE, TOP, M, N, A, LDA )
IMPLICIT NONE
INTEGER CONTXT, M, N, LDA
COMPLEX(kind=kind(0.0D0)) :: A(*)
CHARACTER SCOPE, TOP
error stop 'ZGEBR2D should not be called.'
END SUBROUTINE zgebr2d

SUBROUTINE cgebr2d( CONTXT, SCOPE, TOP, M, N, A, LDA )
IMPLICIT NONE
INTEGER CONTXT, M, N, LDA
COMPLEX A(*)
CHARACTER SCOPE, TOP
error stop 'CGEBR2D should not be called.'
END SUBROUTINE cgebr2d

SUBROUTINE sgebr2d( CONTXT, SCOPE, TOP, M, N, A, LDA )
IMPLICIT NONE
INTEGER CONTXT, M, N, LDA
REAL A(*)
CHARACTER SCOPE, TOP
error stop 'SGEBR2D should not be called.'
END SUBROUTINE sgebr2d

SUBROUTINE dgebr2d( CONTXT, SCOPE, TOP, M, N, A, LDA )
IMPLICIT NONE
INTEGER CONTXT, M, N, LDA
DOUBLE PRECISION A(*)
CHARACTER SCOPE, TOP
error stop 'DGEBR2D should not be called.'
END SUBROUTINE dgebr2d

SUBROUTINE pcgetrs( TRANS, N, NRHS, A, IA, JA, DESCA, IPIV, B, IB, JB, DESCB, INFO )
IMPLICIT NONE
CHARACTER          TRANS
INTEGER            IA, IB, INFO, JA, JB, N, NRHS
INTEGER            DESCA( * ), DESCB( * ), IPIV( * )
COMPLEX            A( * ), B( * )
error stop 'PCGETRS should not be called.'
END SUBROUTINE pcgetrs

SUBROUTINE pzgetrs( TRANS, N, NRHS, A, IA, JA, DESCA, IPIV, B, IB, JB, DESCB, INFO )
IMPLICIT NONE
CHARACTER          TRANS
INTEGER            IA, IB, INFO, JA, JB, N, NRHS
INTEGER            DESCA( * ), DESCB( * ), IPIV( * )
COMPLEX(kind=kind(0.0D0)) ::     A( * ), B( * )
error stop 'PZGETRS should not be called.'
END SUBROUTINE pzgetrs

SUBROUTINE psgetrs( TRANS, N, NRHS, A, IA, JA, DESCA, IPIV, B, IB, JB, DESCB, INFO )
IMPLICIT NONE
CHARACTER          TRANS
INTEGER            IA, IB, INFO, JA, JB, N, NRHS
INTEGER            DESCA( * ), DESCB( * ), IPIV( * )
REAL               A( * ), B( * )
error stop 'PSGETRS should not be called.'
END SUBROUTINE psgetrs

SUBROUTINE pdgetrs( TRANS, N, NRHS, A, IA, JA, DESCA, IPIV, B, IB, JB, DESCB, INFO )
IMPLICIT NONE
CHARACTER          TRANS
INTEGER            IA, IB, INFO, JA, JB, N, NRHS
INTEGER            DESCA( * ), DESCB( * ), IPIV( * )
DOUBLE PRECISION   A( * ), B( * )
error stop 'PDGETRS should not be called.'
END SUBROUTINE pdgetrs

SUBROUTINE pcpotrs( UPLO, N, NRHS, A, IA, JA, DESCA, B, IB, JB, DESCB, INFO )
IMPLICIT NONE
CHARACTER       UPLO
INTEGER         IA, IB, INFO, JA, JB, N, NRHS
INTEGER         DESCA( * ), DESCB( * )
COMPLEX         A( * ), B( * )
error stop 'PCPOTRS should not be called.'
END SUBROUTINE pcpotrs

SUBROUTINE pzpotrs( UPLO, N, NRHS, A, IA, JA, DESCA, B, IB, JB, DESCB, INFO )
IMPLICIT NONE
CHARACTER       UPLO
INTEGER         IA, IB, INFO, JA, JB, N, NRHS
INTEGER         DESCA( * ), DESCB( * )
COMPLEX(kind=kind(0.0D0)) ::     A( * ), B( * )
error stop 'PZPOTRS should not be called.'
END SUBROUTINE pzpotrs

SUBROUTINE pspotrs( UPLO, N, NRHS, A, IA, JA, DESCA, B, IB, JB, DESCB, INFO )
IMPLICIT NONE
CHARACTER       UPLO
INTEGER         IA, IB, INFO, JA, JB, N, NRHS
INTEGER         DESCA( * ), DESCB( * )
REAL            A( * ), B( * )
error stop 'PSPOTRS should not be called.'
END SUBROUTINE pspotrs

SUBROUTINE pdpotrs( UPLO, N, NRHS, A, IA, JA, DESCA, B, IB, JB, DESCB, INFO )
IMPLICIT NONE
CHARACTER       UPLO
INTEGER         IA, IB, INFO, JA, JB, N, NRHS
INTEGER         DESCA( * ), DESCB( * )
DOUBLE          PRECISION A( * ), B( * )
error stop 'PDPOTRS should not be called.'
END SUBROUTINE pdpotrs

SUBROUTINE pscnrm2( N, NORM2, X, IX, JX, DESCX, INCX )
IMPLICIT NONE
INTEGER N, IX, JX, INCX
INTEGER DESCX(*)
REAL NORM2
COMPLEX X( * )
error stop 'PCNRM2 should not be called.'
END SUBROUTINE pscnrm2

SUBROUTINE pdznrm2( N, NORM2, X, IX, JX, DESCX, INCX )
IMPLICIT NONE
INTEGER N, IX, JX, INCX
INTEGER DESCX(*)
DOUBLE PRECISION NORM2
COMPLEX(kind=kind(0.0D0)) :: X( * )
error stop 'PZNRM2 should not be called.'
END SUBROUTINE pdznrm2

SUBROUTINE psnrm2( N, NORM2, X, IX, JX, DESCX, INCX )
IMPLICIT NONE
INTEGER N, IX, JX, INCX
INTEGER DESCX(*)
REAL    NORM2, X( * )
error stop 'PSNRM2 should not be called.'
END SUBROUTINE psnrm2

SUBROUTINE pdnrm2( N, NORM2, X, IX, JX, DESCX, INCX )
IMPLICIT NONE
INTEGER N, IX, JX, INCX
INTEGER DESCX(*)
DOUBLE PRECISION NORM2, X( * )
error stop 'PDNRM2 should not be called.'
END SUBROUTINE pdnrm2

REAL FUNCTION pclange( NORM, M, N, A, IA,  JA, DESCA, WORK )
IMPLICIT NONE
CHARACTER    NORM
INTEGER      IA, JA, M, N
INTEGER      DESCA( * )
COMPLEX      A( * ), WORK( * )
error stop 'PCLANGE should not be called.'
END FUNCTION pclange

DOUBLE PRECISION FUNCTION pzlange( NORM, M, N, A, IA, JA, DESCA, WORK )
IMPLICIT NONE
CHARACTER    NORM
INTEGER      IA, JA, M, N
INTEGER      DESCA( * )
REAL         A( * ), WORK( * )
error stop 'PZLANGE should not be called.'
END FUNCTION pzlange

REAL FUNCTION pslange( NORM, M, N, A, IA,  JA, DESCA, WORK )
IMPLICIT NONE
CHARACTER    NORM
INTEGER      IA, JA, M, N
INTEGER      DESCA( * )
REAL         A( * ), WORK( * )
error stop 'PSLANGE should not be called.'
END FUNCTION pslange

DOUBLE PRECISION FUNCTION pdlange( NORM, M, N, A, IA,  JA, DESCA, WORK )
IMPLICIT NONE
CHARACTER    NORM
INTEGER      IA, JA, M, N
INTEGER      DESCA( * )
DOUBLE       PRECISION A( * ), WORK( * )
error stop 'PDLANGE should not be called.'
END FUNCTION pdlange

SUBROUTINE pcgecon( NORM, N,  A,  IA,  JA,  DESCA,  ANORM, RCOND,  WORK,  LWORK,  IWORK,  LIWORK, INFO )
IMPLICIT NONE
CHARACTER       NORM
INTEGER         IA, INFO, JA, LIWORK, LWORK, N
REAL            ANORM, RCOND
INTEGER         DESCA( * ), IWORK( * )
COMPLEX         A( * ), WORK( * )
error stop 'PCGECON should not be called.'
END SUBROUTINE pcgecon

SUBROUTINE pzgecon( NORM, N,  A,  IA,  JA,  DESCA,  ANORM, RCOND,  WORK,  LWORK,  IWORK,  LIWORK, INFO )
IMPLICIT NONE
CHARACTER       NORM
INTEGER         IA, INFO, JA, LIWORK, LWORK, N
DOUBLE PRECISION ANORM, RCOND
INTEGER         DESCA( * ), IWORK( * )
COMPLEX(kind=kind(0.0D0)) :: A( * ), WORK( * )
error stop 'PZGECON should not be called.'
END SUBROUTINE pzgecon

SUBROUTINE psgecon( NORM, N,  A,  IA,  JA,  DESCA,  ANORM, RCOND,  WORK,  LWORK,  IWORK,  LIWORK, INFO )
IMPLICIT NONE
CHARACTER       NORM
INTEGER         IA, INFO, JA, LIWORK, LWORK, N
REAL            ANORM, RCOND
INTEGER         DESCA( * ), IWORK( * )
REAL            A( * ), WORK( * )
error stop 'PSGECON should not be called.'
END SUBROUTINE psgecon

SUBROUTINE pdgecon( NORM, N,  A,  IA,  JA,  DESCA,  ANORM, RCOND,  WORK,  LWORK,  IWORK,  LIWORK, INFO )
IMPLICIT NONE
CHARACTER       NORM
INTEGER         IA, INFO, JA, LIWORK, LWORK, N
DOUBLE          PRECISION ANORM, RCOND
INTEGER         DESCA( * ), IWORK( * )
DOUBLE          PRECISION A( * ), WORK( * )
error stop 'PDGECON should not be called.'
END SUBROUTINE pdgecon

SUBROUTINE pcgeqpf( M,  N,  A,  IA,  JA, DESCA, IPIV, TAU, WORK, LWORK, INFO )
IMPLICIT NONE
INTEGER    IA, JA, INFO, LWORK, M, N
INTEGER    DESCA( * ), IPIV( * )
COMPLEX    A( * ), TAU( * ), WORK( * )
error stop 'PCGEQPF should not be called.'
END SUBROUTINE pcgeqpf

SUBROUTINE pzgeqpf( M,  N,  A,  IA,  JA, DESCA, IPIV, TAU, WORK, LWORK, INFO )
IMPLICIT NONE
INTEGER    IA, JA, INFO, LWORK, M, N
INTEGER    DESCA( * ), IPIV( * )
COMPLEX(kind=kind(0.0D0)) :: A( * ), TAU( * ), WORK( * )
error stop 'PZGEQPF should not be called.'
END SUBROUTINE pzgeqpf

SUBROUTINE psgeqpf( M,  N,  A,  IA,  JA, DESCA, IPIV, TAU, WORK, LWORK, INFO )
IMPLICIT NONE
INTEGER         IA, JA, INFO, LWORK, M, N
INTEGER         DESCA( * ), IPIV( * )
REAL       A( * ), TAU( * ), WORK( * )
error stop 'PSGEQPF should not be called.'
END SUBROUTINE psgeqpf

SUBROUTINE pdgeqpf( M,  N,  A,  IA,  JA, DESCA, IPIV, TAU, WORK, LWORK, INFO )
IMPLICIT NONE
INTEGER         IA, JA, INFO, LWORK, M, N
INTEGER         DESCA( * ), IPIV( * )
DOUBLE PRECISION A( * ), TAU( * ), WORK( * )
error stop 'PDGEQPF should not be called.'
END SUBROUTINE pdgeqpf

SUBROUTINE pcaxpy(N, A, X, IX, JX, DESCX, INCX, Y, IY, JY, DESCY, INCY)
IMPLICIT NONE
INTEGER N, IX, IY, JX, JY, INCX, INCY
INTEGER DESCX(*), DESCY(*)
COMPLEX A(*),X(*),Y(*)
error stop 'PCAXPY should not be called.'
END SUBROUTINE pcaxpy

SUBROUTINE pzaxpy(N, A, X, IX, JX, DESCX, INCX, Y, IY, JY, DESCY, INCY)
IMPLICIT NONE
INTEGER N, IX, IY, JX, JY, INCX, INCY
INTEGER DESCX(*), DESCY(*)
COMPLEX(kind=kind(0.0D0)) :: A(*),X(*),Y(*)
error stop 'PZAXPY should not be called.'
END SUBROUTINE pzaxpy

SUBROUTINE psaxpy(N, A, X, IX, JX, DESCX, INCX, Y, IY, JY, DESCY, INCY)
IMPLICIT NONE
INTEGER N, IX, IY, JX, JY, INCX, INCY
INTEGER DESCX(*), DESCY(*)
REAL A(*),X(*),Y(*)
error stop 'PSAXPY should not be called.'
END SUBROUTINE psaxpy

SUBROUTINE pdaxpy(N, A, X, IX, JX, DESCX, INCX, Y, IY, JY,  DESCY, INCY)
IMPLICIT NONE
INTEGER N, IX, IY, JX, JY, INCX, INCY
INTEGER DESCX(*), DESCY(*)
DOUBLE PRECISION A(*),X(*),Y(*)
error stop 'PDAXPY should not be called.'
END SUBROUTINE pdaxpy

SUBROUTINE pctrsm ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, IA, JA, DESCA, B, IB, JB, DESCB )
IMPLICIT NONE
CHARACTER          SIDE, UPLO, TRANSA, DIAG
INTEGER            M, N, IA, JA, IB, JB
COMPLEX            ALPHA
INTEGER            DESCA( * ), DESCB( * )
COMPLEX            A( * ), B( * )
error stop 'PCTRSM should not be called.'
END SUBROUTINE pctrsm

SUBROUTINE pztrsm ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, IA, JA, DESCA, B, IB, JB, DESCB )
IMPLICIT NONE
CHARACTER          SIDE, UPLO, TRANSA, DIAG
INTEGER            M, N, IA, JA, IB, JB
COMPLEX(kind=kind(0.0D0)) ::     ALPHA
INTEGER            DESCA( * ), DESCB( * )
COMPLEX(kind=kind(0.0D0)) ::     A( * ), B( * )
error stop 'PZTRSM should not be called.'
END SUBROUTINE pztrsm

SUBROUTINE pstrsm ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, IA, JA, DESCA, B, IB, JB, DESCB )
IMPLICIT NONE
CHARACTER          SIDE, UPLO, TRANSA, DIAG
INTEGER            M, N, IA, JA, IB, JB
REAL               ALPHA
INTEGER            DESCA( * ), DESCB( * )
REAL               A( * ), B( * )
error stop 'PSTRSM should not be called.'
END SUBROUTINE pstrsm

SUBROUTINE pdtrsm ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, IA, JA, DESCA, B, IB, JB, DESCB )
IMPLICIT NONE
CHARACTER          SIDE, UPLO, TRANSA, DIAG
INTEGER            M, N, IA, JA, IB, JB
DOUBLE PRECISION   ALPHA
INTEGER            DESCA( * ), DESCB( * )
DOUBLE PRECISION   A( * ), B( * )
error stop 'PDTRSM should not be called.'
END SUBROUTINE pdtrsm

SUBROUTINE pcunmqr( SIDE,  TRANS,  M,  N,  K,  A,  IA, JA, DESCA, TAU, C, IC,  JC,  DESCC, WORK, LWORK, INFO )
IMPLICIT NONE
CHARACTER SIDE, TRANS
INTEGER   IA, IC, INFO, JA, JC, K, LWORK, M, N
INTEGER   DESCA( * ), DESCC( * )
COMPLEX   A(  *  ), C( * ), TAU( * ), WORK( * )
error stop 'PCUNMQR should not be called.'
END SUBROUTINE pcunmqr

SUBROUTINE pzunmqr( SIDE,  TRANS,  M,  N,  K,  A,  IA, JA, DESCA, TAU, C, IC,  JC,  DESCC,  WORK, LWORK, INFO )
IMPLICIT NONE
CHARACTER SIDE, TRANS
INTEGER   IA, IC, INFO, JA, JC, K, LWORK, M, N
INTEGER   DESCA( * ), DESCC( * )
COMPLEX(kind=kind(0.0D0)) :: A(  *  ), C( * ), TAU( * ), WORK( * )
error stop 'PZUNMQR should not be called.'
END SUBROUTINE pzunmqr

SUBROUTINE psormqr( SIDE,  TRANS,  M,  N,  K,  A,  IA, JA, DESCA, TAU, C, IC,  JC,  DESCC, WORK, LWORK, INFO )
IMPLICIT NONE
CHARACTER SIDE, TRANS
INTEGER   IA, IC, INFO, JA, JC, K, LWORK, M, N
INTEGER   DESCA( * ), DESCC( * )
REAL      A(  *  ), C( * ), TAU( * ), WORK( * )
error stop 'PSORMQR should not be called.'
END SUBROUTINE psormqr

SUBROUTINE pdormqr( SIDE,  TRANS,  M,  N,  K,  A,  IA, JA, DESCA, TAU, C, IC,  JC,  DESCC,  WORK, LWORK, INFO )
IMPLICIT NONE
CHARACTER SIDE, TRANS
INTEGER         IA, IC, INFO, JA, JC, K, LWORK, M, N
INTEGER         DESCA( * ), DESCC( * )
DOUBLE PRECISION  A(  *  ), C( * ), TAU( * ), WORK( * )
error stop 'PDORMQR should not be called.'
END SUBROUTINE pdormqr

SUBROUTINE chk1mat( MA, MAPOS0, NA, NAPOS0, IA, JA, DESCA, DESCAPOS0, INFO )
IMPLICIT NONE
INTEGER            DESCAPOS0, IA, INFO, JA, MA, MAPOS0, NA, NAPOS0
INTEGER            DESCA( * )
error stop 'CHK1MAT should not be called.'
END SUBROUTINE chk1mat

SUBROUTINE pchk2mat( MA, MAPOS0, NA, NAPOS0, IA, JA, DESCA, DESCAPOS0, MB, MBPOS0, NB, NBPOS0, IB, JB, DESCB, DESCBPOS0, &
NEXTRA, EX, EXPOS, INFO )
IMPLICIT NONE
INTEGER            DESCAPOS0, DESCBPOS0, IA, IB, INFO, JA, JB, MA, MAPOS0, MB, MBPOS0, NA, NAPOS0, NB, NBPOS0, NEXTRA
INTEGER            DESCA( * ), DESCB( * ), EX( NEXTRA ), EXPOS( NEXTRA )
error stop 'PCHK2MAT should not be called.'
END SUBROUTINE pchk2mat

SUBROUTINE pxerbla( CONTXT, SRNAME, INFO )
IMPLICIT NONE
INTEGER CONTXT, INFO
CHARACTER SRNAME
error stop 'PXERBLA should not be called.'
END SUBROUTINE pxerbla

SUBROUTINE descset( DESC, M, N, MB, NB, IRSRC, ICSRC, ICTXT, LLD )
IMPLICIT NONE
INTEGER            ICSRC, ICTXT, IRSRC, LLD, M, MB, N, NB
INTEGER            DESC( * )
error stop 'DESCSET should not be called.'
END SUBROUTINE descset
