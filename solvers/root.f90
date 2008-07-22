FUNCTION ROOT(A,B,C,D,SGN)
!---------------------------------------------------------------------
! THIS FUNCTION GIVES BETTER ROOTS OF QUADRATICS BY AVOIDING
! CANCELLATION OF SMALLER ROOT
!---------------------------------------------------------------------
implicit none
real*8 :: root,a, b, c, d, sgn

IF (B*SGN .GE. 0.d0) THEN
  ROOT = -2.d0*C/(B+SGN*SQRT(D))
ELSE
  ROOT = (-B + SGN*SQRT(D)) / (2.d0 * A)
ENDIF
RETURN
END
