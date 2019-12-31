!#####################################
FUNCTION	dtCFL(CFL,Voli,gamma,R,omega)
! Rambod Mojgani
! 90129045
!Computed dt of specific CFL number
IMPLICIT										NONE
DOUBLE PRECISION,DIMENSION(4)					:: omega! Answer (Cell No. × 4)
DOUBLE PRECISION								:: dtCFL,CFL,Voli,gamma,R,P,u,v,density,e,T,a
	density=omega(1)
	u=omega(2)/density
	v=omega(3)/density
	e=omega(4)/density
	P=density*(gamma-1.d0)*(	e	-	5.d-1*(u*u+v*v)	)
	T=P/(density*R)
	a=DSQRT(gamma*R*T)
	dtCFL=CFL*DSQRT(Voli)/DABS(DSQRT(u*u+v*v)+a)
ENDFUNCTION
!#####################################