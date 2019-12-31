!#####################################
SUBROUTINE	Rescheck(nelem,itime,ResCrit,omegao,omegan,IfEXIT)
! Rambod Mojgani
! 90129045
! Checks for residual magnitude and convergence criteria
IMPLICIT										NONE
INTEGER								,INTENT(IN)	:: nelem,itime
DOUBLE PRECISION					,INTENT(IN)	:: ResCrit
DOUBLE PRECISION,DIMENSION(nelem,4)	,INTENT(IN) :: omegao,omegan	! Old & New Answer (Cell No. × 4)
LOGICAL								,INTENT(OUT):: IfEXIT			! True if Converged or Diverged
DOUBLE PRECISION								:: ResRho,ResRhoU,ResRhoV,ResRhoE
INTEGER											:: i				! Counters
ResRho =0.d0
ResRhoU=0.d0
ResRhoV=0.d0
ResRhoE=0.d0
DO i=1,nelem,1
	ResRho =ResRho+DABS(omegan(i,1)-omegao(i,1))
	ResRhoU=ResRho+DABS(omegan(i,2)-omegao(i,2))
	ResRhoV=ResRho+DABS(omegan(i,3)-omegao(i,3))
	ResRhoE=ResRho+DABS(omegan(i,4)-omegao(i,4))
ENDDO

write(10,*) itime,ResRho,ResRhoU,ResRhoV,ResRhoE
write(*,*) itime,REAL(ResRho),REAL(ResRhoU),REAL(ResRhoV),REAL(ResRhoE)
IF(	(ResRho+ResRhoU+ResRhoV+ResRhoE)*2.5d-1	.LT. ResCrit	) IfEXIT=.True.
IF(	ResRho/ResRho/=1.d0	) THEN
	write(*,*) "!!!!!!!!!!!!  Divergence Detected  !!!!!!!!!!!!"
	pause
	IfEXIT=.True.
ENDIF
ENDSUBROUTINE Rescheck
!#####################################