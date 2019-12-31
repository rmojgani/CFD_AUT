PROGRAM RoeSolver
! Rambod Mojgani
! CFD II - Spring 91
! Edition 3.3
! 90129045
! Roe Solver With Reimann Boundary Condition for Subsonic Regime
IMPLICIT										NONE
! ----------------------------------------------:: Input Values
! ------------------------------ Fluid Property
DOUBLE PRECISION								:: R=287.d0			! Gas Constant
DOUBLE PRECISION								:: gamma=1.4		! Gas Specific Heat Ratio
! ------------------------------ Flow Property
DOUBLE PRECISION								:: Minf=1.2d0		! M∞ [---]
DOUBLE PRECISION								:: AoA=0.0d0		! α	 [deg] Angle of attack, Freestream Flow anlge
DOUBLE PRECISION								:: AoR=7.0d0		! A	 [deg] Angle of Rotating the Mesh, Clockwise
DOUBLE PRECISION								:: Tinf=300.d0		! T∞ [K]
DOUBLE PRECISION								:: Pinf=101325.d0	! P∞ [Pa]
! ------------------------------ Solver Tune
DOUBLE PRECISION								:: CFL=5.0d-1		! CFL Criteria Parameter
DOUBLE PRECISION								:: ResCrit=1.d-12	! Mean Residual Stop Criteria
INTEGER											:: itimeResRep=50	! Residual Report Interval
INTEGER											:: ntime=200000		! Time Loop Limit
! ----------------------------------------------:: Other
INTEGER											:: nnode,nelem		! Number of (node		,cell)
INTEGER											:: anode,bnode		! Number of (on airfoil	, on outer boundary)
INTEGER											:: i,j,k,itime		! Counters
INTEGER											:: iL,iR			! Left & Right Cell Designator
INTEGER,DIMENSION(:)			    ,ALLOCATABLE:: ia,ib			! on Airfoil (ia) & on Outer Boundary (ib) Cell Designator Array
INTEGER,DIMENSION(:,:)				,ALLOCATABLE:: TRI,TNB			! (Triangle(Cell No, 123) and neigbour(,456)
DOUBLE PRECISION,DIMENSION(:)		,ALLOCATABLE:: x,y				! Node Location
DOUBLE PRECISION,DIMENSION(:)		,ALLOCATABLE:: xRot,yRot		! Node Location (Rotated with AoR)
DOUBLE PRECISION,DIMENSION(:,:)		,ALLOCATABLE:: omegao,omegan	! Answer (Cell No. × 4): omega(nelem,4)
DOUBLE PRECISION	,DIMENSION(:,:) ,ALLOCATABLE:: nx,ny			! Normal Vectors ( Cell No.×3 )
DOUBLE PRECISION	,DIMENSION(:)   ,ALLOCATABLE:: CellArea			! Area of each cell: CellArea(nelem)
DOUBLE PRECISION	,DIMENSION(:,:) ,ALLOCATABLE:: EdgeLen			! Edges length: EdgeLen(nelem,3)
DOUBLE PRECISION,DIMENSION(4,3)					:: Fn				! Cell Flux Array
DOUBLE PRECISION,DIMENSION(4,4)					:: Rn
DOUBLE PRECISION,DIMENSION(4)					:: landa,dwbar,FL,FR,omegainf
DOUBLE PRECISION								:: dtonvol			!Δt/Cellarea(i)
DOUBLE PRECISION								:: DEG2RAD,distance,triarea,RoeM,dtCFL,dt=1.0d-7
DOUBLE PRECISION								:: densityinf,uinf,vinf,einf,hinf,qinf
DOUBLE PRECISION								:: rhowall,uwall,vwall,ewall,Pwall
DOUBLE PRECISION								:: rho,u,v,e,p
INTEGER,DIMENSION(4)							:: ilist=(/1,2,3,1/)
LOGICAL											:: IfEXIT=.False.	! True if Converged or Diverged

OPEN(unit=10,file='convergence.plt',status='replace')
write(10,*) 'TITLE="Residual"'
write(10,*) 'VARIABLES="iloop" "<GREEK>r</GREEK>"	"<GREEK>r</GREEK>u"	"<GREEK>r</GREEK>v"	"<GREEK>r</GREEK>e"	'

! Conditioning
AoA=DEG2RAD(AoA)
uinf=Minf*sqrt(gamma*R*Tinf)*cos(AoA)
vinf=Minf*sqrt(gamma*R*Tinf)*sin(AoA)
densityinf=Pinf/(R*Tinf)
einf=(Pinf/(gamma-1.d0)/densityinf+5.d-1*(uinf*uinf+vinf*vinf))
hinf=einf+Pinf/densityinf
qinf=5.d-1*densityinf*(uinf*uinf+vinf*vinf)
!write(*,*)	(gamma-1.d0)/R*(einf-0.5*(uinf*uinf+vinf*vinf)),Tinf

! Information
write(*,*) 'AoA          = ',	REAL(AoA*3.14156d0/180.d0)
write(*,*) 'M_inf        = ',	REAL(Minf)
write(*,*) 'U_inf        = ',	REAL(uinf)
write(*,*) 'V_inf        = ',	REAL(vinf)
write(*,*) 'Density_inf  = ',	REAL(densityinf)
write(*,*) 'CFL	         = ',	REAL(CFL)

! Enter File DATA from table into program STARTED
OPEN(unit=1,file='MESH8.txt',status='old',ACTION='read')
	read(1,*),nnode
	read(1,*),nelem
	read(1,*),anode
	read(1,*),bnode
	ALLOCATE(x(nnode))
	ALLOCATE(y(nnode))
	ALLOCATE(xRot(nnode))
	ALLOCATE(yRot(nnode))
	ALLOCATE(TRI(nelem,3))
	ALLOCATE(TNB(nelem,3))
	ALLOCATE(ia(anode))
	ALLOCATE(ib(bnode))
	DO i=1,nnode,1
		read(1,*),x(i),y(i)
	ENDDO
	DO i=1,nelem,1
		read(1,*),TRI(i,1),TRI(i,2),TRI(i,3),TNB(i,1),TNB(i,2),TNB(i,3)
	ENDDO
	DO i=1,anode,1
		read(1,*),j,ia(i)
	ENDDO
	DO i=1,bnode,1
		read(1,*),j,ib(i)
	ENDDO
		
	ALLOCATE(omegao(nelem,4))
	ALLOCATE(omegan(nelem,4))
	ALLOCATE(nx(nelem,3))
	ALLOCATE(ny(nelem,3))
	ALLOCATE(CellArea(nelem))
	ALLOCATE(EdgeLen(nelem,3))
	j=0 ! Dummy
close(1)
!Entering File DATA into program is FINISHED

! Rotating the Mesh
AoR=DEG2RAD(-AoR)
DO i=1,nnode,1
	CALL RotateXY(x(i),y(i),AoR,xRot(i),yRot(i))
ENDDO
DO i=1,nnode,1
	x(i)=xRot(i)
	y(i)=yRot(i)
ENDDO

! Calculating Normal Vectors
DO i=1,nelem,1
	DO j=1,3,1
		CALL	NormalVector(x(TRI(i,j)),y(TRI(i,j)) , x(TRI(i,ilist(j+1))) ,y(TRI(i,ilist(j+1))),nx(i,j),ny(i,j))
	ENDDO
ENDDO

! Renaming Outer Boundary Neighbor Color (ColorL -10== On INLET, -20==On OUTLET)
! Numerical Test showed it postepone convergence, CFL effect yet can be tested
DO i=1,bnode,1
	IF	  (	TNB(ib(i),1)==0	) THEN
		IF (	DOT_PRODUCT((/uinf,vinf/),(/nx(ib(i),1),ny(ib(i),1)/))	.LE. 0.d0) THEN
			TNB(ib(i),1)=-10	! On INLET
		ELSE
			TNB(ib(i),1)=-20	! On OUTLET
		ENDIF
	ELSEIF(	TNB(ib(i),2)==0	) THEN
		IF (	DOT_PRODUCT((/uinf,vinf/),(/nx(ib(i),2),ny(ib(i),2)/))	.LE. 0.d0) THEN
			TNB(ib(i),2)=-10	! On INLET
		ELSE
			TNB(ib(i),2)=-20	! On OUTLET
		ENDIF
	ELSEIF(	TNB(ib(i),3)==0	) THEN
		IF (	DOT_PRODUCT((/uinf,vinf/),(/nx(ib(i),3),ny(ib(i),3)/))	.LE. 0.d0) THEN
			TNB(ib(i),3)=-10	! On INLET
		ELSE
			TNB(ib(i),3)=-20	! On OUTLET
		ENDIF
	ENDIF
ENDDO	!DO i=1,bnode,1

! Calculating Each Cell area
DO i=1,nelem,1
	CellArea(i)=triarea(	(/x(TRI(i,1)),x(TRI(i,2)),x(TRI(i,3))/) , (/y(TRI(i,1)),y(TRI(i,2)),y(TRI(i,3))/)	)
ENDDO
! Calculating Edges Length in each Cell
DO i=1,nelem,1
	DO j=1,3,1
		EdgeLen(i,j)=distance(x(TRI(i,j)),y(TRI(i,j)),x(TRI(i,ilist(j+1))),y(TRI(i,ilist(j+1))))
	ENDDO
ENDDO

! Initializing
omegainf(1)=densityinf
omegainf(2)=densityinf*uinf
omegainf(3)=densityinf*vinf
omegainf(4)=densityinf*einf

DO	i=1,nelem,1
	omegao(i,1)=omegainf(1)
	omegao(i,2)=omegainf(2)
	omegao(i,3)=omegainf(3)
	omegao(i,4)=omegainf(4)
ENDDO


DO itime=1,ntime,1										! Time March Loop
	DO i=1,nelem,1										! Loop on All Elements,Start
		iL=i											! Left Cell (Flux On Calculating) Index
		dt=dtCFL(CFL,CellArea(iL),gamma,R,omegao(iL,:)) ! Finding each Δt
		dtonvol=dt/CellArea(iL)							! Finding each Δt/Cell Area

		! On INLET	::::::::::::::::::::::::::::::::::::::::::::::::::::::::
		IF(	TNB(iL,1)==-10	 .OR.	TNB(iL,2)==-10	.OR.	TNB(iL,3)==-10)THEN			
			
			IF(Minf.GT. 1.d0)	THEN	! Equal to Boundary @ "Subsonic" inlet
				CALL	Inlet(omegainf,omegan(iL,:))
			ELSE						! Equal to Boundary @ "Supersonic" inlet
				IF(TNB(iL,1)==-10)	THEN
					CALL	RiemannInlet(omegainf,omegao(TNB(iL,2),:),omegao(TNB(iL,3),:),gamma,R,nx(iL,1),ny(iL,1),omegan(iL,:))
				ELSEIF(TNB(iL,2)==-10) THEN
					CALL	RiemannInlet(omegainf,omegao(TNB(iL,1),:),omegao(TNB(iL,3),:),gamma,R,nx(iL,2),ny(iL,2),omegan(iL,:))
				ELSEIF(TNB(iL,3)==-10)	THEN
					CALL	RiemannInlet(omegainf,omegao(TNB(iL,1),:),omegao(TNB(iL,2),:),gamma,R,nx(iL,3),ny(iL,3),omegan(iL,:))
				ENDIF
			ENDIF

		! On OUTLET	::::::::::::::::::::::::::::::::::::::::::::::::::::::::
		ELSEIF(	TNB(iL,1)==-20	)THEN		! On OUTLET
			omegan(iL,1)=5.d-1*(omegao(TNB(iL,2),1)+omegao(TNB(iL,3),1))
			omegan(iL,2)=5.d-1*(omegao(TNB(iL,2),2)+omegao(TNB(iL,3),2))
			omegan(iL,3)=5.d-1*(omegao(TNB(iL,2),3)+omegao(TNB(iL,3),3))
			omegan(iL,4)=5.d-1*(omegao(TNB(iL,2),4)+omegao(TNB(iL,3),4))

		ELSEIF(	TNB(i,2)==-20	)THEN		! On OUTLET
			omegan(iL,1)=5.d-1*(omegao(TNB(iL,1),1)+omegao(TNB(iL,3),1))
			omegan(iL,2)=5.d-1*(omegao(TNB(iL,1),2)+omegao(TNB(iL,3),2))
			omegan(iL,3)=5.d-1*(omegao(TNB(iL,1),3)+omegao(TNB(iL,3),3))
			omegan(iL,4)=5.d-1*(omegao(TNB(iL,1),4)+omegao(TNB(iL,3),4))

		ELSEIF(	TNB(iL,3)==-20	)THEN		! On OUTLET
			omegan(iL,1)=5.d-1*(omegao(TNB(iL,1),1)+omegao(TNB(iL,2),1))
			omegan(iL,2)=5.d-1*(omegao(TNB(iL,1),2)+omegao(TNB(iL,2),2))
			omegan(iL,3)=5.d-1*(omegao(TNB(iL,1),3)+omegao(TNB(iL,2),3))
			omegan(iL,4)=5.d-1*(omegao(TNB(iL,1),4)+omegao(TNB(iL,2),4))

		! On Airfoil::::::::::::::::::::::::::::::::::::::::::::::::::::::::
		ELSEIF(	TNB(i,1)== 0	)THEN		! On Airfoil

			DO j=1,3,1 ! Loop on a Cell, Start
				iR=TNB(iL,j)
				CALL RoeMat(gamma,nelem,iL,iR,omegao,nx(iL,j),ny(iL,j),Fn(:,j))
			ENDDO	! Loop on a Cell, End

			CALL	ExtractProp(gamma,omegao(iL,:),rhowall,uwall,vwall,ewall,Pwall)

			omegan(iL,1)=omegao(iL,1)-dtonvol*(								  Fn(1,2)*EdgeLen(iL,2)	+	Fn(1,3)*EdgeLen(iL,3)	)
			omegan(iL,2)=omegao(iL,2)-dtonvol*(	Pwall*nx(iL,1)*EdgeLen(iL,1)+ Fn(2,2)*EdgeLen(iL,2)	+	Fn(2,3)*EdgeLen(iL,3)	)
			omegan(iL,3)=omegao(iL,3)-dtonvol*(	Pwall*ny(iL,1)*EdgeLen(iL,1)+ Fn(3,2)*EdgeLen(iL,2)	+	Fn(3,3)*EdgeLen(iL,3)	)
			omegan(iL,4)=omegao(iL,4)-dtonvol*(								  Fn(4,2)*EdgeLen(iL,2)	+	Fn(4,3)*EdgeLen(iL,3)	)

		ELSEIF(	TNB(iL,2)== 0	)THEN		! On Airfoil

			DO j=1,3,1 ! Loop on a Cell, Start
				iR=TNB(iL,j)
				!write(*,*) i,':',iL,iR,'edge:',j
				CALL RoeMat(gamma,nelem,iL,iR,omegao,nx(iL,j),ny(iL,j),Fn(:,j))
			ENDDO	! Loop on a Cell, End

			CALL	ExtractProp(gamma,omegao(iL,:),rhowall,uwall,vwall,ewall,Pwall)

			omegan(iL,1)=omegao(iL,1)-dtonvol*(	Fn(1,1)*EdgeLen(iL,1)								+	Fn(1,3)*EdgeLen(iL,3)	)
			omegan(iL,2)=omegao(iL,2)-dtonvol*(	Fn(2,1)*EdgeLen(iL,1)+Pwall*nx(iL,2)*EdgeLen(iL,2)	+	Fn(2,3)*EdgeLen(iL,3)	)
			omegan(iL,3)=omegao(iL,3)-dtonvol*(	Fn(3,1)*EdgeLen(iL,1)+Pwall*ny(iL,2)*EdgeLen(iL,2)	+	Fn(3,3)*EdgeLen(iL,3)	)
			omegan(iL,4)=omegao(iL,4)-dtonvol*(	Fn(4,1)*EdgeLen(iL,1)								+	Fn(4,3)*EdgeLen(iL,3)	)
		
		ELSEIF(	TNB(iL,3)== 0	)THEN		! On Airfoil

			DO j=1,3,1 ! Loop on a Cell, Start
				iR=TNB(iL,j)
				!write(*,*) i,':',iL,iR,'edge:',j
				CALL RoeMat(gamma,nelem,iL,iR,omegao,nx(iL,j),ny(iL,j),Fn(:,j))
			ENDDO	! Loop on a Cell, End

			CALL	ExtractProp(gamma,omegao(iL,:),rhowall,uwall,vwall,ewall,Pwall)

			omegan(iL,1)=omegao(iL,1)-dtonvol*(	Fn(1,1)*EdgeLen(iL,1)	+	Fn(1,2)*EdgeLen(iL,2)								)
			omegan(iL,2)=omegao(iL,2)-dtonvol*(	Fn(2,1)*EdgeLen(iL,1)	+	Fn(2,2)*EdgeLen(iL,2)+Pwall*nx(iL,3)*EdgeLen(iL,3)	)
			omegan(iL,3)=omegao(iL,3)-dtonvol*(	Fn(3,1)*EdgeLen(iL,1)	+	Fn(3,2)*EdgeLen(iL,2)+Pwall*ny(iL,3)*EdgeLen(iL,3)	)
			omegan(iL,4)=omegao(iL,4)-dtonvol*(	Fn(4,1)*EdgeLen(iL,1)	+	Fn(4,2)*EdgeLen(iL,2)								)

		! In Domain ::::::::::::::::::::::::::::::::::::::::::::::::::::::::
		ELSE
			DO j=1,3,1 ! Loop on a Cell, Start
				iR=TNB(i,j)
				CALL RoeMat(gamma,nelem,iL,iR,omegao,nx(iL,j),ny(iL,j),Fn(:,j))
			ENDDO	! Loop on a Cell, End
			CALL	EulerExp(dtonvol,omegao(iL,:),Fn,EdgeLen(iL,:),omegan(iL,:))
		ENDIF
	ENDDO	!	DO i=1,nelem,1
	
	!Residual Check
	IF(MOD(itime,itimeResRep)==0)THEN ! To Lessen the calculation and writing time, every "itimeResRep" it is calculeted and reported
		CALL	Rescheck(nelem,itime,ResCrit,omegao,omegan,IfEXIT)	
		IF(IfExit) EXIT
	ENDIF
	omegao=omegan
ENDDO	!	DO itime=1,10,1
close(10)

CALL Writepack(nnode,nelem,anode,gamma,omegao,TRI,ia,x,y,pinf,qinf)

ENDPROGRAM RoeSolver