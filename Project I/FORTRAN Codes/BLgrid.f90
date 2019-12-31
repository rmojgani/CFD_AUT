! ####################################################
SUBROUTINE gridgenBL(dens,IE,JE,thx,thy,lx,ly,x,y,ICORN)
! Rambod Mojgani - 91.3.--
! Grid Generation - ED 3.0
IMPLICIT										 NONE
DOUBLE PRECISION					,INTENT(in)	 :: dens
INTEGER								,INTENT(in)	 :: IE, JE
DOUBLE PRECISION					,INTENT(in)	 :: thx,thy ! in Degree
DOUBLE PRECISION					,INTENT(in)	 :: lx,ly
DOUBLE PRECISION,		PARAMETER				 :: PI=ACOS(0.d0)*2
DOUBLE PRECISION,DIMENSION(IE*JE)	,INTENT(out) :: x,y
INTEGER,DIMENSION((IE-1)*(JE-1),4)	,INTENT(out) :: ICORN
INTEGER											 :: i,j,k,cj,ci,IEJE
DOUBLE PRECISION,DIMENSION(4)					 :: XO,YO
IEJE=IE*JE
thx=thx*PI/180
thy=thy*PI/180
! Outer Point locations:
! 2---------1
! |---------|
! |---------|
! 3---------4
! End point locations
XO=(/lx,0.d0,0.d0,lx/)
YO=(/ly,ly,0.d0,0.d0/)
i=1
j=1
k=1
DO cj=1,JE,1
	DO ci=1,IE,1
		CALL HypInt2((/XO(3),XO(4)/),lx/IE*dens,lx/IE*dens,IE-1,i-1,x(k))
		CALL HypInt2((/YO(3),YO(2)/),ly/IE*dens,ly/IE*dens,JE-1,j-1,y(k))
		i=i+1
		k=k+1
	END DO
	i=1
	j=j+1
END DO

i=1
j=1
DO cj=1,JE-1,1
	DO ci=1,IE-1,1
		ICORN(j,3)=i
		ICORN(j,4)=ICORN(j,3)+1
		ICORN(j,2)=ICORN(j,3)+IE
		ICORN(j,1)=ICORN(j,2)+1
		i=i+1
		j=j+1
	END DO
	i=i+1
END DO
END SUBROUTINE gridgenBL
! ####################################################