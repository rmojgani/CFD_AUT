! ---#####################################################
SUBROUTINE HypInt2(rref,ds1,ds2,n,k,r)
! Rambod Mojgani
IMPLICIT									NONE
DOUBLE PRECISION, DIMENSION(2), INTENT(in) :: rref
DOUBLE PRECISION, INTENT(in)			   :: ds1,ds2
INTEGER, INTENT(in)						   :: n,k
DOUBLE PRECISION, INTENT(out)			   :: r
INTEGER									   :: i
DOUBLE PRECISION						   :: a,b,d,s,u,kn
kn=1.d0*k/n
a=SQRT(ds2/ds1)
b=1.d0/(n*SQRT(ds2*ds1))
d=2.d0
DO i=1,5000,1 
	d=d-(SIN(d)-b*d)/(COS(d)-b)
END DO

u=5.d-1*(1+TANH(d*(kn-5.d-1))/TANH(d/2))
s=u/(a+(1-a)*u)
r=rref(1)+(rref(2)-rref(1))*s

ENDSUBROUTINE HypInt2