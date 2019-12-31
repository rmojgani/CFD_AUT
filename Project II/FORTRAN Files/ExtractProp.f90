!#####################################
SUBROUTINE	ExtractProp(gamma,omega,rho,u,v,e,p)
! Rambod Mojgani
! 90129045
! Calculated primary variables from euler standard variables
IMPLICIT										NONE
DOUBLE PRECISION					,INTENT(IN) :: gamma
DOUBLE PRECISION,DIMENSION(4)		,INTENT(IN) :: omega	! Answer (4)
DOUBLE PRECISION					,INTENT(OUT):: rho,u,v,e,p
	rho=omega(1)
	u=omega(2)/rho
	v=omega(3)/rho
	e=omega(4)/rho
	p=rho*(gamma-1.d0)*( e	- 5.d-1*(u*u+v*v)	)
ENDSUBROUTINE ExtractProp
!#####################################