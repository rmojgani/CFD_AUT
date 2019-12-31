! ####################################################
 SUBROUTINE upwind(Xi,Yi,u,v,cup,dsup)
! Rambod Mojgani
! 91.02.28, ED 1.1
! 91.03.02, ED 1.2 : Adding dsup
! ##---------------------------------------##
! ## Accompanied SUBROUTINS & FUNCTIONS :: ##
! ##     SUBROUTINE solve2l(a,b,c,x,y)	   ##
! ##     FUNCTION distance(x1,y1,x2,y2)	   ##
! ##---------------------------------------##
! Upwind Calculation, Findes the First Order Upwind Location and 
! Linearly Interpolates the value on the point, this is proposed 
! to avoid "False Diffusion" whilst the flow has angle to the grid.
IMPLICIT NONE
DOUBLE PRECISION,DIMENSION(4),INTENT(in)					  :: Xi,Yi				! Corners Coordinates
DOUBLE PRECISION,DIMENSION(4),INTENT(in)					  :: u,v				! Velocity Componnets on Integral Points (j)
DOUBLE PRECISION,DIMENSION(4,4),INTENT(out)					  :: cup				! CUP as in Dr. Karimian HW # 2
DOUBLE PRECISION,DIMENSION(4),INTENT(out)					  :: dsup
DOUBLE PRECISION,DIMENSION(4)								  :: Xj,Yj
DOUBLE PRECISION,DIMENSION(4,4)								  :: x,y,l1,l2,l
DOUBLE PRECISION								 			  :: m1,m2,xmid,ymid,distance
DOUBLE PRECISION,DIMENSION(4)								  :: xup,yup,iiloc
DOUBLE PRECISION,DIMENSION(4,2)								  :: xupdn,yupdn
INTEGER,DIMENSION(4,2)										  :: iloc
DOUBLE PRECISION,DIMENSION(2)								  :: a,b,c
INTEGER														  :: i,j,k
INTEGER,DIMENSION(5)										  :: nlist
nlist=(/1,2,3,4,1/)
! Finding the Cell's Center point location
xmid=(Xi(1)+Xi(2)+Xi(3)+Xi(4))*2.5d-1
ymid=(Yi(1)+Yi(2)+Yi(3)+Yi(4))*2.5d-1
! Finding the j(=1,4,1) points location
DO i=1,4,1
	Xj(i)=5.d-1*(xmid+(Xi(nlist(i))+Xi(nlist(i+1)))*5.d-1)
	Yj(i)=5.d-1*(ymid+(Yi(nlist(i))+Yi(nlist(i+1)))*5.d-1)
END DO
DO i=1,4,1
	m1=Yi(nlist(i+1))-Yi(nlist(i))
	m1=m1/(Xi(nlist(i+1))-Xi(nlist(i)))
	a(1)=-m1
	b(1)=1.d0
	c(1)=-m1*Xi(nlist(i))+Yi(nlist(i))
	DO j=1,4,1
		IF (u(j).eq.0.d0) THEN
			x(i,j)=Xj(j)
			y(i,j)=m1*x(i,j)+c(1)
		ELSE
			m2=v(j)/u(j)
			a(2)=-m2
			b(2)=1.d0
			c(2)=-m2*Xj(nlist(j))+Yj(nlist(j))
			CALL solve2l(a,b,c,x(i,j),y(i,j))
			IF ((Xi(nlist(i+1))-Xi(nlist(i))).eq.0.d0) THEN
				x(i,j)=Xi(nlist(i))
				y(i,j)=m2*x(i,j)+c(2)
			ENDIF
		ENDIF
	END DO
END DO
! --- Select The Points Which are "On the Cell"
DO j=1,4,1
	DO i=1,4,1
		l1(i,j)=distance(x(i,j),y(i,j),Xi(nlist(i)),Yi(nlist(i)))
		l2(i,j)=distance(x(i,j),y(i,j),Xi(nlist(i+1)),Yi(nlist(i+1)))
		l(i,j)=distance(Xi(nlist(i)),Yi(nlist(i)),Xi(nlist(i+1)),Yi(nlist(i+1)))
	ENDDO
ENDDO
!write(*,*) X
!write(*,*) Y
DO j=1,4,1
	k=1
	DO i=1,4,1
		IF ((l1(i,j)+l2(i,j)) == l(i,j)) THEN
			xupdn(j,k)=x(i,j)
			yupdn(j,k)=y(i,j)
			iloc(j,k)=i ! ADD Comment ####
			k=k+1
		ENDIF
	ENDDO
ENDDO
! --- Select The Upstream Point, Based on U & V Direction for each "j" point
DO j=1,4,1
	IF (u(j).lt.0.d0) THEN
		xup(j)=MAX(xupdn(j,1),xupdn(j,2))
	ELSEIF (u(j).gt.0.d0) THEN
		xup(j)=MIN(xupdn(j,1),xupdn(j,2))
	ELSEIF (v(j).lt.0.d0) THEN
		yup(j)=MAX(yupdn(j,1),yupdn(j,2))
	ELSE
		yup(j)=MIN(yupdn(j,1),yupdn(j,2))
	ENDIF
	IF (u(j) /= 0.d0) THEN
		IF (xup(j) == xupdn(j,1)) THEN
			yup(j) = yupdn(j,1)
			iiloc(j)=iloc(j,1)
		ELSEIF (xup(j) == xupdn(j,2)) THEN
			yup(j) = yupdn(j,2)
			iiloc(j)=iloc(j,2)
		ENDIF
	ELSE
		IF (yup(j) == yupdn(j,1)) THEN
			xup(j) = xupdn(j,1)
			iiloc(j)=iloc(j,1)
		ELSEIF (yup(j) == yupdn(j,2)) THEN
			xup(j) = xupdn(j,2)
			iiloc(j)=iloc(j,2)
		ENDIF	
	ENDIF
ENDDO
! Calculating The cup
DO j=1,4,1
	i=iiloc(j)
	cup(i,j)=(l(i,j)-l1(i,j))/l(i,j)
	IF (i+1==5) THEN
		cup(1,j)=(l(i,j)-l2(i,j))/l(i,j)
	ELSE
		cup(i+1,j)=(l(i,j)-l2(i,j))/l(i,j)
	END IF
ENDDO
! Calculating The dsup
DO i=1,4,1
	dsup(i)=distance(xup(i),yup(i),Xj(i),Yj(i))
END DO
END SUBROUTINE upwind
! ####################################################