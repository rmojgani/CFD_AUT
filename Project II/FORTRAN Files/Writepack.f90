SUBROUTINE Writepack(nnode,nelem,anode,gamma,omegao,TRI,ia,x,y,pinf,qinf)
IMPLICIT										NONE
INTEGER								,INTENT(IN)	:: nnode,nelem,anode
DOUBLE PRECISION					,INTENT(IN)	:: gamma
INTEGER,DIMENSION(nelem,3)			,INTENT(IN)	:: TRI
INTEGER,DIMENSION(anode)			,INTENT(IN)	:: ia
DOUBLE PRECISION,DIMENSION(nelem,4)	,INTENT(IN)	:: omegao
DOUBLE PRECISION,DIMENSION(nnode)	,INTENT(IN) :: x,y
DOUBLE PRECISION					,INTENT(IN)	:: pinf,qinf	
DOUBLE PRECISION								:: rho,u,v,e,p,xmid,ymid
INTEGER											:: i,j	!Counters

OPEN(unit=1,file='Density.plt',status='replace')
	write(1,*) 'VARIABLES="X","Y","Density"'
	write(1,*)  'ZONE NODES=',nnode,', ELEMENTS=',nelem
	write(1,*)  ', DATAPACKING=BLOCK, ZONETYPE=FETRIANGLE'
	write(1,*)  ', varlocation=([',nelem,',',3,']=CELLCENTERED)'
	write(1,*) REAL(x)
	write(1,*) REAL(y)
	DO i=1,nelem,1
		write(1,*) REAL(omegao(i,1))
	ENDDO
	write(1,*)
	DO i=1,nelem,1
		write(1,*) TRI(i,:)
	ENDDO
close(1)

OPEN(unit=1,file='Mach.plt',status='replace')
	write(1,*) 'VARIABLES="X","Y","Mach"'
	write(1,*)  'ZONE NODES=',nnode,', ELEMENTS=',nelem
	write(1,*)  ', DATAPACKING=BLOCK, ZONETYPE=FETRIANGLE'
	write(1,*)  ', varlocation=([',nelem,',',3,']=CELLCENTERED)'
	write(1,*) REAL(x)
	write(1,*) REAL(y)
	DO i=1,nelem,1
		CALL	ExtractProp(gamma,omegao(i,:),rho,u,v,e,p)
		write(1,*) REAL(DSQRT(u*u+v*v)/DSQRT(gamma*p/rho))
	ENDDO
	write(1,*)
	DO i=1,nelem,1
		write(1,*) TRI(i,:)
	ENDDO
close(1)


OPEN(unit=1,file='u.plt',status='replace')
	write(1,*) 'VARIABLES="X","Y","u [m/s]"'
	write(1,*)  'ZONE NODES=',nnode,', ELEMENTS=',nelem
	write(1,*)  ', DATAPACKING=BLOCK, ZONETYPE=FETRIANGLE'
	write(1,*)  ', varlocation=([',nelem,',',3,']=CELLCENTERED)'
	write(1,*) REAL(x)
	write(1,*) REAL(y)
	DO i=1,nelem,1
		write(1,*) REAL(omegao(i,2)/omegao(i,1))
	ENDDO
	write(1,*)
	DO i=1,nelem,1
		write(1,*) TRI(i,:)
	ENDDO
close(1)
OPEN(unit=1,file='v.plt',status='replace')
	write(1,*) 'VARIABLES="X","Y","v [m/s]"'
	write(1,*)  'ZONE NODES=',nnode,', ELEMENTS=',nelem
	write(1,*)  ', DATAPACKING=BLOCK, ZONETYPE=FETRIANGLE'
	write(1,*)  ', varlocation=([',nelem,',',3,']=CELLCENTERED)'
	write(1,*) REAL(x)
	write(1,*) REAL(y)
	DO i=1,nelem,1
		write(1,*) REAL(omegao(i,3)/omegao(i,1))
	ENDDO
	write(1,*)
	DO i=1,nelem,1
		write(1,*) TRI(i,:)
	ENDDO
close(1)

OPEN(unit=1,file='P.plt',status='replace')
	write(1,*) 'VARIABLES="X","Y","Pressure [Pa]"'
	!write(1,*) 'VARIABLES="X","Y","Density","u","v","e"'
	write(1,*)  'ZONE NODES=',nnode,', ELEMENTS=',nelem
	write(1,*)  ', DATAPACKING=BLOCK, ZONETYPE=FETRIANGLE'
	write(1,*)  ', varlocation=([',nelem,',',3,']=CELLCENTERED)'
	!write(1,*)  ', varlocation=([',nelem,',',6,']=CELLCENTERED)'
	write(1,*) REAL(x)
	write(1,*) REAL(y)
	DO i=1,nelem,1
		write(1,*) REAL(	(gamma-1.0)*(omegao(i,4)-0.5/omegao(i,1)*(omegao(i,2)**2+omegao(i,3)**2))	)
	ENDDO
	write(1,*)
	DO i=1,nelem,1
		write(1,*) TRI(i,:)
	ENDDO
close(1)


OPEN(unit=1,file='CP-U.plt',status='replace')!upper
OPEN(unit=2,file='CP-L.plt',status='replace')!lower
OPEN(unit=3,file='CP.plt',status='replace')!lower
	write(1,*) 'TITLE="Predicted C<sub>p</sub> on Airfoil Upper Surface"'
	write(1,*) 'VARIABLES="x" "C<sub>p"'
	write(2,*) 'TITLE="Predicted C<sub>p</sub> on Airfoil Lower Surface"'
	write(2,*) 'VARIABLES="x" "C<sub>p"'
	write(3,*) 'TITLE="Predicted C<sub>p</sub> on Airfoil Surface"'
	write(3,*) 'VARIABLES="x" "C<sub>p"'
	DO i=1,anode,1
		CALL	ExtractProp(gamma,omegao(ia(i),:),rho,u,v,e,p)
		xmid=(	x(TRI(ia(i),1))+x(TRI(ia(i),2))+x(TRI(ia(i),3))	)/3.d0-4.0
		ymid=(	y(TRI(ia(i),1))+y(TRI(ia(i),2))+y(TRI(ia(i),3))	)/3.d0-4.0
		IF(ymid .GE.0.d0) write(1,*) xmid,(p-pinf)/qinf
		IF(ymid .LE.0.d0) write(2,*) xmid,(p-pinf)/qinf
						  write(3,*) xmid,(p-pinf)/qinf
	ENDDO
close(1)
close(2)
close(3)

ENDSUBROUTINE Writepack