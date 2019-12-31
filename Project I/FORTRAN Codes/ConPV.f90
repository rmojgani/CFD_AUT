! ####################################################
SUBROUTINE CovPV(IE,JE,B1,B2,eps)
! Rambod Mojgani
IMPLICIT										 NONE
INTEGER								,INTENT(in)	 :: IE,JE
DOUBLE PRECISION,DIMENSION(3*IE*JE)	,INTENT(in)	 :: B1,B2
DOUBLE PRECISION,DIMENSION(3)		,INTENT(out) :: eps
DOUBLE PRECISION,DIMENSION(IE*JE)				 :: um,vm,pm
DOUBLE PRECISION,DIMENSION(IE*JE)				 :: uom,vom,pom
INTEGER											 :: i,j,IEJE
DOUBLE PRECISION								 :: epsP,epsU,epsV
IEJE=IE*JE
eps=0.d0
! Separating u,v,p just for ease of handeling the parameters
DO i=1,IEJE,1
	pom(i)=B1(3*i-2)
	uom(i)=B1(3*i-1)
	vom(i)=B1(3*i)
END DO
DO i=1,IEJE,1
	pm(i)=B2(3*i-2)
	um(i)=B2(3*i-1)
	vm(i)=B2(3*i)
END DO

j=1
DO i=1,IEJE,1
	IF (pm(i).ne.0.d0) THEN
		epsP=epsP+DABS(pm(i)-pom(i))/DABS(pm(i))
		j=j+1
	END	IF
END	DO
epsP=epsP/j

j=1
DO i=1,IEJE,1
	IF (um(i).ne.0.d0) THEN
		epsU=epsU+DABS(um(i)-uom(i))/DABS(um(i))
		j=j+1
	END	IF
END	DO
epsU=epsU/j

j=1
DO i=1,IEJE,1
	IF (vm(i).ne.0.d0) THEN
		epsV=epsV+DABS(vm(i)-vom(i))/DABS(vm(i))
		j=j+1
	END	IF
END	DO
epsV=epsV/j
eps=(/epsP,epsU,epsV/)
ENDSUBROUTINE CovPV