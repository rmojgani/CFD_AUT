!#####################################
SUBROUTINE	EulerExp(dtONvoli,omegao,Fn,ds,omegan)
! Rambod Mojgani
! 90129045
! Calculates Euler Explicit when all three edges of cell have flux
IMPLICIT										NONE
DOUBLE PRECISION					,INTENT(IN)	:: dtONvoli
DOUBLE PRECISION,DIMENSION(4,3)		,INTENT(IN) :: Fn		! Answer (4 × 3)
DOUBLE PRECISION,DIMENSION(4)		,INTENT(IN) :: omegao	! Answer (4)
DOUBLE PRECISION,DIMENSION(3)		,INTENT(IN) :: ds		! Answer (3)
DOUBLE PRECISION,DIMENSION(4)		,INTENT(OUT):: omegan	! Answer (4)
	omegan(1)=omegao(1)-dtONvoli*(	Fn(1,1)*ds(1)	+	Fn(1,2)*ds(2)	+	Fn(1,3)*ds(3)	)
	omegan(2)=omegao(2)-dtONvoli*(	Fn(2,1)*ds(1)	+	Fn(2,2)*ds(2)	+	Fn(2,3)*ds(3)	)
	omegan(3)=omegao(3)-dtONvoli*(	Fn(3,1)*ds(1)	+	Fn(3,2)*ds(2)	+	Fn(3,3)*ds(3)	)
	omegan(4)=omegao(4)-dtONvoli*(	Fn(4,1)*ds(1)	+	Fn(4,2)*ds(2)	+	Fn(4,3)*ds(3)	)
ENDSUBROUTINE EulerExp
!#####################################