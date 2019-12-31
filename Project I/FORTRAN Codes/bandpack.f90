! ##---------------------------------------##
! ##   Numerical Recipe Band Solver Pack   ##
! ##		Organized for ease of use	   ##
! ##			  Rambod Mojgani		   ##
! ##     SUBROUTINE  banbks(.......)	   ##
! ##     SUBROUTINE  bandec(.......)	   ##
! ##---------------------------------------##
! ####################################################
SUBROUTINE bandpack(n,m1,m2,np,mp,mpl,A,B)
! A X + B = 0 
IMPLICIT										 NONE
INTEGER								,INTENT(in)	 :: n,m1,m2,np,mp,mpl
DOUBLE PRECISION, Dimension(1:n,mp)	,INTENT(in)	 :: A
INTEGER, Dimension(1:n)			,INTENT(inout)	 :: B
INTEGER, Dimension(1:n)							 :: indx
DOUBLE PRECISION,Dimension(n,m1)				 :: al
DOUBLE PRECISION								 :: d
INTEGER											 :: i
CALL bandec(A,n,m1,m2,np,mp,al,mpl,indx,d)
CALL banbks(A,n,m1,m2,np,mp,al,mpl,indx,B)
ENDSUBROUTINE bandpack
! ####################################################
SUBROUTINE banbks(a,n,m1,m2,np,mp,al,mpl,indx,b)
! Numerical Recipe
IMPLICIT NONE
INTEGER m1,m2,mp,mpl,n,np,indx(n)
DOUBLE PRECISION a(np,mp),al(np,mpl),b(n)
INTEGER i,k,l,mm
DOUBLE PRECISION dum
mm=m1+m2+1
if(mm.gt.mp.or.m1.gt.mpl.or.n.gt.np) pause 'bad args in banbks'
l=m1
do  k=1,n
	i=indx(k)
	if(i.ne.k)then
		dum=b(k)
		b(k)=b(i)
		b(i)=dum
	endif
    if(l.lt.n)l=l+1
    do  i=k+1,l
		b(i)=b(i)-al(k,i-k)*b(k)
	enddo
enddo
l=1
do  i=n,1,-1
	dum=b(i)
	do  k=2,l
		dum=dum-a(i,k)*b(k+i-1)
	enddo
	b(i)=dum/a(i,1)
	if(l.lt.mm) l=l+1
enddo
ENDSUBROUTINE banbks
! ####################################################
SUBROUTINE bandec(a,n,m1,m2,np,mp,al,mpl,indx,d)
! Numerical Recipe
IMPLICIT NONE
INTEGER :: m1,m2,mp,mpl,n,np,indx(n)
DOUBLE PRECISION :: d,a(np,mp),al(np,mpl),TINY
PARAMETER (TINY=1.e-20)
INTEGER :: i,j,k,l,mm
DOUBLE PRECISION :: dum
mm=m1+m2+1
if(mm.gt.mp.or.m1.gt.mpl.or.n.gt.np) pause 'bad args in bandec'
l=m1
do  i=1,m1
	do  j=m1+2-i,mm
		a(i,j-l)=a(i,j)
	enddo
	l=l-1
	do  j=mm-l,mm
		a(i,j)=0.d0
	enddo
enddo
d=1.d0
l=m1
do  k=1,n
	dum=a(k,1)
	i=k
	if(l.lt.n)l=l+1
	do  j=k+1,l
		if(abs(a(j,1)).gt.abs(dum))then
			dum=a(j,1)
			i=j
		endif
	enddo
	indx(k)=i
	if(dum.eq.0.) a(k,1)=TINY
	if(i.ne.k)then
		d=-d
		do  j=1,mm
			dum=a(k,j)
			a(k,j)=a(i,j)
			a(i,j)=dum
		enddo
	endif
	do  i=k+1,l
		dum=a(i,1)/a(k,1)
		al(k,i-k)=dum
		do  j=2,mm
			a(i,j-1)=a(i,j)-dum*a(k,j)
		enddo
		a(i,mm)=0.d0
	enddo
enddo
ENDSUBROUTINE bandec