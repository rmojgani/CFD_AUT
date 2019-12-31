program main
! Rambod Mojgani
! Cavity - ED 7.0 - 91.03.24
IMPLICIT										 NONE
! Problem Set
INTEGER, PARAMETER								 :: IE=31, JE=31
INTEGER											 :: niter=15,ntime=10
DOUBLE PRECISION								 :: dt=1.d1
DOUBLE PRECISION, PARAMETER						 :: lx=1.d0,ly=1.d0
DOUBLE PRECISION, PARAMETER						 :: thx=0.d0,thy=0.d0
DOUBLE PRECISION, PARAMETER						 :: density=1.d0,U=1.d0,Rey=100.d0
DOUBLE PRECISION, PARAMETER						 :: Tdiffusivity=1.75d-3
DOUBLE PRECISION, PARAMETER						 :: viscosity=(density*U*lx)/Rey
DOUBLE PRECISION, DIMENSION(4)					 :: Temp=(/100.d0,100.d0,500.d0,100.d0/) !(/Up, Left , Down , Right/)
CHARACTER		, PARAMETER						 :: testdiff='F' ! Do you want to test diffusion term? F for no, T for yes
CHARACTER		, PARAMETER						 :: BLrGrid ='F' ! Do you want use Boundary Layer Grid? F for no, T for yes
DOUBLE PRECISION, PARAMETER						 :: dens=10.0d0  ! Grid Wall denser parameter
! Code Parameters
INTEGER											 :: i,iloop=0,itime,iter ! Counters
INTEGER											 :: nelem=(IE-1)*(JE-1),IEJE=IE*JE,IEJE3=IE*JE*3
INTEGER, DIMENSION((IE-1)*(JE-1),4)				 :: ICORN
DOUBLE PRECISION, PARAMETER						 :: PI=ACOS(0.d0)*2
DOUBLE PRECISION, PARAMETER						 :: Ux=U*COS(thx*PI/180.d0)
DOUBLE PRECISION, PARAMETER						 :: Uy=U*SIN(thx*PI/180.d0)
DOUBLE PRECISION, DIMENSION(IE*JE)				 :: x,y
DOUBLE PRECISION, DIMENSION(IE*JE*3,6*IE+11)	 :: ABand
DOUBLE PRECISION, DIMENSION(IE*JE*3)			 :: BLagged,BOld,B,BR
DOUBLE PRECISION,DIMENSION(IE*JE,2*IE+3)		 :: AT
DOUBLE PRECISION,DIMENSION(IE*JE)				 :: BT=0.d0 ! Nulify
DOUBLE PRECISION								 :: dx,dy
DOUBLE PRECISION								 :: minerror=1.d-12
DOUBLE PRECISION,DIMENSION(3)					 :: eps
! List Generator  Parameters:
INTEGER,DIMENSION(4)							 :: index
INTEGER,DIMENSION(2*IE+2*(JE-2))				 :: list
!------------------------------------

! Information
write(*,*) 'Grid size = ','			Reynolds Number ='
write(*,*)	IE,'x',JE,					REAL(density*U*lx/viscosity)
write(*,*) 'Maximum Number of each time step iteration:',niter
write(*,*) 'Solution End time:',REAL(ntime*dt),'Seconds'
IF(BLrGrid=='T') THEN
	write(*,*) '------------------------------------'
	write(*,*) 'Note: You have chosen Boundary Layer Grid	'
	write(*,*) '	this code can perform orthogonal cavity in this case'
	write(*,*) '------------------------------------'
ENDIF
write(*,*) '------------------------------------'
write(*,*) 'Press Any Key to Confirm The above'
write(*,*) '------------------------------------'
pause
! Openin Files
OPEN(unit=10,file='convergence.plt',status='replace')
write(10,*) 'TITLE="Residual"'
write(10,*) 'VARIABLES="iloop" "Pressure" "U<sub>x</sub>" "U<sub>y</sub>"'
! Initializing
DO i=1,IEJE,1
	BOld(3*i-2)=1.d0	! P
	BOld(3*i-1)=Ux*0.d0	! U
	BOld(3*i  )=Uy*0.d0	! V
END DO

! ::: To generate the grid
IF(BLrGrid=='T') THEN
	! ::: To generate Boundary Layer Grid
	CALL gridgenBL(dens,IE,JE,thx,thy,lx,ly,x,y,ICORN)		
ELSE
	! ::: To generate Uniform Grid
	CALL gridgen(IE,JE,thx,thy,lx,ly,x,y,ICORN)
ENDIF
! ::: To generate the nodes on the boundary
CALL generatelist(IE,JE,list,index)

DO itime=1,ntime,1		! ::: The Time Marche Loop
    Blagged=BOld
    DO iter=1,niter,1	! ::: The Iteration Convergence Loop
        iloop=iloop+1
        ABand=0.d0
        BR=0.d0
        write(*,*) 'Iteration:',iter,'/',niter,'---',itime,'/',ntime

        ! ::: Assembly process of Continuity & X,Y- Momentum Equations
        CALL assemb(Blagged,BOld,ICORN,nelem,IE,JE,x,y,dt,density,viscosity,ABand,BR)	
        ! ::: Enforece Boundary Condition on Pressure-Velocity Field
		CALL EnforceBCCavity(IE,JE,list,index,Ux,Uy,ABand,BR)
		! ::: Solving Pressure-Velocity Filed
        CALL bandpack(IEJE3,3*IE+5,3*IE+5,IEJE3,6*IE+11,3*IE+5,ABand,BR)
		! ::: Convergence Check (Iteration)
		CALL CovPV(IE,JE,Blagged,BR,eps)
		write(*,*) REAL(eps)
        write(10,*) iloop,REAL(eps)
        IF(sum(eps)/3.d0.le.minerror) EXIT
        Blagged=BR
    ENDDO
	! ::: Convergence Check (Time)
	CALL CovPV(IE,JE,BOld,Blagged,eps)
	write(*,*) '--------------Time:',REAL(itime*dt),'---------------'
    write(*,*) REAL(eps)
    write(10,*) iloop+1,REAL(eps)
	IF(sum(eps)/3.d0.le.minerror) EXIT
	BOld=Blagged

ENDDO

IF(testdiff=='T') BOld=0.d0
! ::: Assembly process of Energy Equation
CALL AssemEnergy(BOld,ICORN,nelem,IE,JE,x,y,density,viscosity,Tdiffusivity,AT)
! ::: Enforece Boundary Condition on Temprature Field
CALL BCT(IE,JE,list,index,Temp,AT,BT)
! ::: Solve Energy on converged velocity-pressure domain 
CALL bandpack(IEJE,IE+1,IE+1,IEJE,2*IE+3,IE+1,AT,BT)

! ::: Writing the Results
CALL writeresult(BOld,BT,IE,JE,x,y)
CALL resultonline(BOld,BT,IE,JE,x,y)
CALL shearwrite(BOld,IE,JE,x,y,U)
CALL Nusseltwrite(BT,IE,JE,x,y,Temp(1),Temp(3))
! ::: Closing Files
close(10)
write(*,*) '------------------------------------'
write(*,*) '		END of solution'
write(*,*) '------------------------------------'
pause
ENDPROGRAM main