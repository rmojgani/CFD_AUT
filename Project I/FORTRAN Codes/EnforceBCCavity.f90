! ####################################################
SUBROUTINE EnforceBCCavity(IE,JE,list,index,ux,vy,A,B)
! Rambod Mojgani 91.03.19; ED 2.0
IMPLICIT										 NONE
INTEGER								,INTENT(in)		:: IE,JE
INTEGER,DIMENSION(2*IE+2*(JE-2))	,INTENT(in)		:: list
INTEGER,DIMENSION(4)				,INTENT(in)		:: index
DOUBLE PRECISION								,INTENT(in)		:: ux,vy
DOUBLE PRECISION,DIMENSION(IE*JE*3,6*IE+11)		,INTENT(inout)	:: A
DOUBLE PRECISION,DIMENSION(IE*JE*3)				,INTENT(inout)	:: B
INTEGER												:: i,j,ii,jj,k

!-- Boundary Force on MATRIX A, banded
DO k=1,2*IE+2*(JE-2),1
	i=list(k)
	DO j=1,IE*JE*3,1
		CALL BandReady(i,j,3*IE+5,3*IE+5,ii,jj)
		IF (jj<=6*IE+11.AND.jj.ge.1) THEN
			A(3*ii-1,jj)=0.d0
			A(3*ii  ,jj)=0.d0
		ENDIF
	END DO
	CALL BandReady(3*i-1,3*i-1,3*IE+5,3*IE+5,ii,jj)
	A(ii,jj)=1.d0
	CALL BandReady( 3*i , 3*i ,3*IE+5,3*IE+5,ii,jj)
	A(ii,jj)=1.d0
END DO


DO j=1,6*IE+11,1
	A(1,j)=0.d0
END DO
CALL BandReady(1,1,3*IE+5,3*IE+5,ii,jj)
A(ii,jj)=1.d0
B( 1 )=0.d0

! -- Set Boundary Condition
! Up Row
DO k=index(1)+1,index(2),1
	i=list(k)
	B(3*i-1)=ux
	B(3*i  )=vy
END DO
! Lower Row
DO k=1,index(1),1
	i=list(k)
	B(3*i-1)=0.d0
	B(3*i  )=0.d0
END DO
! Left Row
DO k=index(2)+1,index(3),1
	i=list(k)
	B(3*i-1)=0.d0
	B(3*i  )=0.d0
END DO
! Right Row
DO k=index(3)+1,index(4),1
	i=list(k)
	B(3*i-1)=0.d0
	B(3*i  )=0.d0
END DO

END SUBROUTINE EnforceBCCavity
! ####################################################
