!#####################################
SUBROUTINE	Inlet(omegainf,omegainlet)
! Rambod Mojgani
! 90129045
! Enforce Inlet Boundary Condition (Normal)
IMPLICIT										NONE
DOUBLE PRECISION,DIMENSION(4)		,INTENT(IN) :: omegainf! Answer (Cell No. × 4)
DOUBLE PRECISION,DIMENSION(4)		,INTENT(OUT):: omegainlet! Answer (Cell No. × 4)
	omegainlet(1)=omegainf(1)
	omegainlet(2)=omegainf(2)
	omegainlet(3)=omegainf(3)
	omegainlet(4)=omegainf(4)
ENDSUBROUTINE	Inlet
!#####################################
SUBROUTINE	RiemannInlet(omegainf,omega1,omega2,gamma,R,nx,ny,omegainlet)
! Rambod Mojgani
! 90129045
! Enforce Inlet Boundary Condition (Reimann)
IMPLICIT										NONE
DOUBLE PRECISION,DIMENSION(4)		,INTENT(IN) :: omegainf,omega1,omega2	! Answer (Cell No. × 4)
DOUBLE PRECISION					,INTENT(IN) :: gamma,R,nx,ny
DOUBLE PRECISION,DIMENSION(4)		,INTENT(OUT):: omegainlet				! Answer (Cell No. × 4)
DOUBLE PRECISION								:: rhoinf,uinf,vinf,einf,pinf,rho1,u1,v1,e1,p1,rho2,u2,v2,e2,p2
DOUBLE PRECISION								:: Tinf,T1,T2,Cinf,C1,C2,Cint,uninf,utinf,uninf,utinf,uint,vint,unint
DOUBLE PRECISION								:: unB,CB,TB,rhoB,utB,nx2ny2,uB,vB,pB,eB

CALL	ExtractProp(gamma,omegainf,rhoinf,uinf,vinf,einf,pinf)
CALL	ExtractProp(gamma,omega1,rho1,u1,v1,e1,p1)
CALL	ExtractProp(gamma,omega2,rho2,u2,v2,e2,p2)

Tinf=(gamma-1.d0)/R*(einf-0.5*(uinf*uinf+vinf*vinf))
T1  =(gamma-1.d0)/R*( e1 -0.5*(u1*u1+v1*v1))
T2  =(gamma-1.d0)/R*( e2 -0.5*(u2*u2+v2*v2))

Cinf=DSQRT(gamma*R*Tinf)
C1=DSQRT(gamma*R*T1)
C2=DSQRT(gamma*R*T2)
Cint=5.d-1*(C1+C2)



uninf=DOT_PRODUCT((/uinf,vinf/),(/nx,ny/))
utinf=DOT_PRODUCT((/uinf,vinf/),(/-ny,nx/))

uint=5.d-1*(u1+u2)
vint=5.d-1*(v1+v2)
unint=DOT_PRODUCT((/uint,vint/),(/nx,ny/))

unB=5.d-1*(uninf+unint)-(Cinf-Cint)/(gamma-1.d0)
CB=-(gamma-1.d0)/4.d0*(uninf-unint)+(Cinf+Cint)/(2.d0)
TB=CB*CB/gamma/R
rhoB=rhoinf*(TB/Tinf)**(1/(gamma-1.d0))
utB=utinf

nx2ny2=nx*nx+ny*ny
uB=(nx*unB-ny*utB)/nx2ny2
vB=(ny*unB+nx*utB)/nx2ny2
pB=rhoB*R*TB
eB=(pB/(gamma-1.d0)/rhoB+5.d-1*(uB*uB+vB*vB))

omegainlet(1)=rhoB
omegainlet(2)=rhoB*uB
omegainlet(3)=rhoB*vB
omegainlet(4)=rhoB*eB

ENDSUBROUTINE	RiemannInlet