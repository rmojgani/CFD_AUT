! ####################################################
SUBROUTINE BCT(IE,JE,list,index,Temp,AT,BT)
! Rambod Mojgani 91.03.23; ED 3.0
IMPLICIT										 NONE
INTEGER								,INTENT(in)		:: IE,JE
INTEGER,DIMENSION(2*IE+2*(JE-2))	,INTENT(in)		:: list
INTEGER,DIMENSION(4)				,INTENT(in)		:: index
DOUBLE PRECISION,DIMENSION(4)	    ,INTENT(in)	 :: Temp !(/Up, Left , Down , Right/)
DOUBLE PRECISION,DIMENSION(IE*JE,2*IE+3)	    ,INTENT(inout)	 :: AT
DOUBLE PRECISION,DIMENSION(IE*JE)			    ,INTENT(inout)	 :: BT
INTEGER												:: i,j,ii,jj,k
!-- Boundary Force on MATRIX A, banded
!-- Boundary Force on MATRIX dd
DO k=1,2*IE+2*(JE-2),1
	i=list(k)
	DO j=1,IE*JE,1
		CALL BandReady(i,j,IE+1,IE+1,ii,jj)
		IF (jj<=2*IE+3.AND.jj>=1)	AT(ii,jj)=0.d0
		CALL BandReady(i,i,IE+1,IE+1,ii,jj)
		IF (jj<=2*IE+3.AND.jj>=1)	AT(ii,jj)=1.d0
	END DO
END DO


! -- Set Boundary Condition
 !(/Up, Left , Down , Right/)
! Lower Row
DO k=1,index(1),1
	i=list(k)
	BT(i)=Temp(3)
END DO
! Up Row
DO k=index(1)+1,index(2),1
	i=list(k)
	BT(i)=Temp(1)
END DO
! Left Row
DO k=index(2)+1,index(3),1
	i=list(k)
	BT(i)=Temp(2)
END DO
! Right Row
DO k=index(3)+1,index(4),1
	i=list(k)
	BT(i)=Temp(4)
END DO

END SUBROUTINE BCT
! ####################################################
