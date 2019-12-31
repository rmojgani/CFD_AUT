! ####################################################
SUBROUTINE generatelist(IE,JE,list,index)
! Rambod Mojgani
! Generates the list of location of boundary points number
IMPLICIT NONE
INTEGER,INTENT(in) :: IE,JE
INTEGER,DIMENSION(2*IE+2*(JE-2)),INTENT(out) :: list
INTEGER,DIMENSION(4),INTENT(out) :: index
INTEGER :: i,j,k,IEJE
IEJE=IE*JE
k=0
DO i=1,IE,1
	k=k+1
	list(k)=i
END DO
index(1)=k
DO i=IEJE-IE+1,IEJE,1
	k=k+1
	list(k)=i
END DO
index(2)=k
j=1
DO i=1,JE-2,1
	k=k+1
	list(k)=j*IE+1
	j=j+1
END DO
index(3)=k
j=2
DO i=1,JE-2,1
	k=k+1
	list(k)=j*IE
	j=j+1
END DO
index(4)=k
END SUBROUTINE generatelist
! ####################################################