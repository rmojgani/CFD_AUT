SUBROUTINE derv(Xi,Yi,volinp,dsx,dsy,l2d)
! Rambod Mojgani
! Calculates the geometric properties of grid
IMPLICIT											NONE
DOUBLE PRECISION,DIMENSION(4)			,INTENT(in) :: Xi,Yi				! Corners Coordinates
DOUBLE PRECISION,Dimension(4)			,INTENT(in)	:: volinp
DOUBLE PRECISION,Dimension(4,2)			,INTENT(out):: dsx,dsy
DOUBLE PRECISION,Dimension(4)			,INTENT(out):: l2d
DOUBLE PRECISION								    :: deta2,dksi2
INTEGER,DIMENSION(4)							    :: klist,j2list
INTEGER											    :: i,k

j2list=(/2,3,4,1/)
klist=(/4,1,2,3/)

dsx=0.d0
dsy=0.d0
l2d=0.d0
DO i=1,4,1
	k=j2list(i)
	dsx(i,1)=sum(Yi(:))*2.5d-1 - (Yi(i)+Yi(k))*5.d-1
ENDDO
DO i=1,4,1
	k=klist(i)
	dsx(i,2)=-dsx(k,1)
ENDDO
DO i=1,4,1
	k=j2list(i)
	dsy(i,1)=-( sum(Xi(:))*2.5d-1 - (Xi(i)+Xi(k))*5.d-1	)
ENDDO
DO i=1,4,1
	k=klist(i)
	dsy(i,2)=-dsy(k,1)
ENDDO
DO i=1,4,1
	deta2  = dsx(i,1)*dsx(i,1) + dsy(i,1)*dsy(i,1)
	dksi2  = volinp(i)*volinp(i) / deta2
	l2d(i) = 1.d0 / (2.d0/dksi2 + 8.d0/(3.d0*deta2) )
ENDDO
!---
ENDSUBROUTINE derv