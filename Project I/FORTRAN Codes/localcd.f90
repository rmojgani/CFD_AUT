! ####################################################
SUBROUTINE localcd(Xi,Yi,dndx,dndy,g,cd)
! Rambod Mojgani
IMPLICIT											NONE
DOUBLE PRECISION,DIMENSION(4)			,INTENT(in) :: Xi,Yi				! Corners Coordinates
DOUBLE PRECISION,Dimension(4,4)			,INTENT(in)	:: dndx,dndy
DOUBLE PRECISION						,INTENT(in)	:: g
DOUBLE PRECISION,Dimension(4,4)			,INTENT(out):: cd
DOUBLE PRECISION,Dimension(4)						:: ddy,ddx
INTEGER												:: i,j

ddx=(/(Xi(4)-Xi(2))*5.d-1,(Xi(1)-Xi(3))*5.d-1,(Xi(2)-Xi(4))*5.d-1,(Xi(3)-Xi(1))*5.d-1/)
ddy=(/(Yi(4)-Yi(2))*5.d-1,(Yi(1)-Yi(3))*5.d-1,(Yi(2)-Yi(4))*5.d-1,(Yi(3)-Yi(1))*5.d-1/)

DO j=1,4,1
	DO i=1,4,1
		cd(i,j)=-g*(dndx(i,j)*ddy(i)-dndy(i,j)*ddx(i))
	ENDDO
ENDDO

END SUBROUTINE localcd
! ####################################################