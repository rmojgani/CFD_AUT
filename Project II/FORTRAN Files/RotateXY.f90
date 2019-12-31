!#####################################
SUBROUTINE	RotateXY(x,y,AoR,xn,yn)
! Rambod Mojgani
! 90129045
! Rotate X and Y
IMPLICIT										NONE
DOUBLE PRECISION					,INTENT(IN)	:: x,y,AoR
DOUBLE PRECISION					,INTENT(OUT):: xn,yn
		xn=DCOS(AoR)*x-DSIN(AoR)*y
		yn=DSIN(AoR)*x+DCOS(AoR)*y
ENDSUBROUTINE RotateXY
!#####################################