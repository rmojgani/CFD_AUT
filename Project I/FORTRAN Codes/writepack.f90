! ####################################################
SUBROUTINE writeresult(BPV,BT,IE,JE,x,y)
IMPLICIT										 NONE
DOUBLE PRECISION,DIMENSION(3*iE*JE),INTENT(in)	 :: BPV
DOUBLE PRECISION,DIMENSION(IE*JE)	,INTENT(in)	 :: BT
INTEGER								,INTENT(in)	 :: IE,JE
DOUBLE PRECISION,DIMENSION(IE*JE)	,INTENT(in)	 :: x,y
INTEGER											 :: IEJE,i,j,IEJE3
IEJE=IE*JE
OPEN(unit=1,file='Contour.plt',status='replace')
	write(1,*) 'VARIABLES="X" "Y" "U<sub>x" "V<sub>y" "Pressure" "Temprature"'
	write(1,*) 'Zone i=',IE,',j=',JE,',k=1,f=point'
	DO i=1,IEJE,1
		write(1,*) x(i),y(i),BPV(3*i-1),BPV(3*i),BPV(3*i-2),BT(i)
	END DO
close(1)
ENDSUBROUTINE writeresult
! ####################################################
SUBROUTINE resultonline(BPV,BT,IE,JE,x,y)
IMPLICIT										 NONE
DOUBLE PRECISION,DIMENSION(3*IE*JE) ,INTENT(in)	 :: BPV
DOUBLE PRECISION,DIMENSION(IE*JE) ,INTENT(in)	 :: BT
INTEGER							   ,INTENT(in)	 :: IE,JE
DOUBLE PRECISION,DIMENSION(IE*JE) ,INTENT(in)	 :: x,y
INTEGER											 :: IEJE,IEJE3,i,imid,jmid
IEJE=IE*JE
imid=(IE+1)/2
jmid=(JE+1)/2

OPEN(unit=1,file='VCenterlineu.plt',status='replace')
write(1,*) 'TITLE="u along Vertical centerline"'
write(1,*) 'VARIABLES="u" "y"'
DO i=0,IEJE-1,IE
	write(1,*) BPV(3*(imid+i)-1),y(imid+i)
END DO
close(1)


OPEN(unit=2,file='HCenterlineV.plt',status='replace')
write(2,*) 'TITLE="v along Horiziontal centerline"'
write(2,*) 'VARIABLES="x" "v"'
DO i=(jmid-1)*IE+1,(jmid-1)*IE+IE,1
	write(2,*) x(i),BPV(3*i)
END DO
close(2)

OPEN(unit=3,file='WallPressure.plt',status='replace')
write(3,*) 'TITLE="Predicted C<sub>p</sub> along top wall"'
write(3,*) 'VARIABLES="x" "C<sub>p"'
DO i=(JE-1)*IE+1,IEJE,1
	write(3,*) x(i),BPV(3*i-2)
END DO
close(3)

ENDSUBROUTINE resultonline
! ####################################################
SUBROUTINE shearwrite(BPV,IE,JE,x,y,U)
IMPLICIT										 NONE
DOUBLE PRECISION,DIMENSION(3*IE*JE) ,INTENT(in)	 :: BPV
INTEGER								,INTENT(in)	 :: IE,JE
DOUBLE PRECISION,DIMENSION(IE*JE)	,INTENT(in)	 :: x,y
DOUBLE PRECISION					,INTENT(in)	 :: U
INTEGER											 :: IEJE,i
DOUBLE PRECISION								 :: dn
IEJE=IE*JE
dn=y((JE-2)*IE+1)-y((JE-3)*IE+1) ! Normal Distance of wall adjacent cell

OPEN(unit=1,file='dUdn.plt',status='replace')
write(1,*) 'TITLE="dU/dn"'
write(1,*) 'VARIABLES="x" "dU/dn"'
DO i=(JE-2)*IE+1,IEJE-IE,1
	write(1,*) x(i),(sqrt(BPV(3*i-1)**2+BPV(3*i)**2)-U)/dn
END DO
close(1)
ENDSUBROUTINE shearwrite
! ####################################################
SUBROUTINE Nusseltwrite(BT,IE,JE,x,y,Twall,Tinf)
IMPLICIT										 NONE
DOUBLE PRECISION,DIMENSION(IE*JE) ,INTENT(in)	 :: BT
INTEGER								,INTENT(in)	 :: IE,JE
DOUBLE PRECISION,DIMENSION(IE*JE)	,INTENT(in)	 :: x,y
DOUBLE PRECISION					,INTENT(in)	 :: Twall,Tinf
INTEGER											 :: IEJE,i,j,IEJE3
DOUBLE PRECISION								 :: dn
IEJE=IE*JE
dn=y((JE-2)*IE+1)-y((JE-3)*IE+1) ! Normal Distance of wall adjacent cell
OPEN(unit=1,file='Nusselt.plt',status='replace')
write(1,*) 'TITLE="dU/dn"'
write(1,*) 'VARIABLES="x" "dT/dn"'
DO i=(JE-2)*IE+1,IEJE-IE,1
	write(1,*) x(i),-(BT(i)-Twall)/dn/(Twall-Tinf)
END DO
close(1)
ENDSUBROUTINE Nusseltwrite