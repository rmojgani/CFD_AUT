!#####################################
SUBROUTINE	RoeMat(gamma,nelem,iL,iR,omega,nx,ny,Fn)
! Rambod Mojgani
! 90129045
! Roe Scheme 
IMPLICIT										NONE
DOUBLE PRECISION					,INTENT(IN)	:: gamma
INTEGER								,INTENT(IN)	:: nelem
INTEGER								,INTENT(IN)	:: iL,iR
DOUBLE PRECISION,DIMENSION(nelem,4)	,INTENT(IN)	:: omega	! Answer (Cell No. × 4)
DOUBLE PRECISION					,INTENT(IN)	:: nx,ny	! Normal Vectors Componnent
DOUBLE PRECISION,DIMENSION(4)		,INTENT(OUT):: Fn
INTEGER											:: i		! Counters
DOUBLE PRECISION,DIMENSION(4)					:: landa,dwbar,FL,FR,omegaR
DOUBLE PRECISION,DIMENSION(4,4)					:: Rn
DOUBLE PRECISION								:: rhoL,rhoR,uL,uR,vL,vR,unL,unR,utL,utR,pL,pR,eL,eR,hsL,hsR,sigma
DOUBLE PRECISION								:: RoeM,rhobar,ubar,vbar,hbar,cbar,unbar,vnbar,utbar

CALL	ExtractProp(gamma,omega(iL,:),rhoL,uL,vL,eL,pL)
CALL	ExtractProp(gamma,omega(iR,:),rhoR,uR,vR,eR,pR)

unL=+uL*nx+vL*ny
utL=-uL*ny+vL*nx
hsL=eL+pL/rhoL

unR=+uR*nx+vR*ny
utR=-uR*ny+vR*nx
hsR=eR+pR/rhoR

sigma=DSQRT(rhoL)/(	DSQRT(rhoL)+DSQRT(rhoR)	)
rhobar=DSQRT(rhoL*rhoR)
unbar =RoeM(unL,unR,sigma)
utbar =RoeM(utL,utR,sigma)
ubar  =RoeM(uL ,uR ,sigma)
vbar  =RoeM(vL ,vR ,sigma)
hbar  =RoeM(hsL,hsR,sigma)
cbar  =(gamma-1.d0)*(hbar-5.d-1*(ubar*ubar+vbar*vbar))
cbar  =DSQRT(cbar)

dwbar(1)=(+5.d-1/cbar/cbar)*((pR-pL)-rhobar*cbar*(unR-unL))
dwbar(2)=rhobar*(utR-utL)
dwbar(3)=(-1.d0/cbar/cbar)*((pR-pL)-cbar*cbar*(rhoR-rhoL))
dwbar(4)=(+5.d-1/cbar/cbar)*((pR-pL)+rhobar*cbar*(unR-unL))

Rn(1,1:4)=(/1.d0			,0.d0	 ,1.d0						 ,1.d0			 /)
Rn(2,1:4)=(/ubar-cbar*nx	,-ny	 ,ubar						 ,ubar+cbar*nx	 /)
Rn(3,1:4)=(/vbar-cbar*ny	,+nx	 ,vbar						 ,vbar+cbar*ny	 /)
Rn(4,1:4)=(/hbar-cbar*unbar	,utbar	 ,5.d-1*(ubar*ubar+vbar*vbar),hbar+cbar*unbar/)

landa	 =(/unbar-cbar		,unbar		 ,unbar						 ,unbar+cbar	 /)

FL=(/rhoL*unL,rhoL*unL*uL+pL*nx,rhoL*unL*vL+pL*ny,rhoL*unL*hsL/)
FR=(/rhoR*unR,rhoR*unR*uR+pR*nx,rhoR*unR*vR+pR*ny,rhoR*unR*hsR/)

Fn=0.0
DO i=1,4,1
	Fn(1)=Fn(1)-5.d-1*DABS(landa(i))*dwbar(i)*Rn(1,i)
	Fn(2)=Fn(2)-5.d-1*DABS(landa(i))*dwbar(i)*Rn(2,i)
	Fn(3)=Fn(3)-5.d-1*DABS(landa(i))*dwbar(i)*Rn(3,i)
	Fn(4)=Fn(4)-5.d-1*DABS(landa(i))*dwbar(i)*Rn(4,i)
ENDDO
Fn(1)=5.d-1*(FL(1)+FR(1))+Fn(1)
Fn(2)=5.d-1*(FL(2)+FR(2))+Fn(2)
Fn(3)=5.d-1*(FL(3)+FR(3))+Fn(3)
Fn(4)=5.d-1*(FL(4)+FR(4))+Fn(4)

ENDSUBROUTINE	RoeMat
!#####################################
FUNCTION	RoeM(L,R,sigma)
! Rambod Mojgani
! 90129045
! Roe Mean Value 
IMPLICIT										NONE
DOUBLE PRECISION								:: L,R,sigma,RoeM
	RoeM=L*sigma+R*(1.d0-sigma)
ENDFUNCTION
!#####################################