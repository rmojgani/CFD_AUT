! ####################################################
SUBROUTINE gridgen(IE,JE,thx,thy,lx,ly,x,y,ICORN)
! Rambod Mojgani - 90.12.--
! Grid Generation - ED 2.0
IMPLICIT										 NONE
INTEGER								,INTENT(in)	 :: IE, JE
DOUBLE PRECISION					,INTENT(in)	 :: thx,thy ! in Degree
DOUBLE PRECISION					,INTENT(in)	 :: lx,ly
DOUBLE PRECISION,		PARAMETER				 :: PI=ACOS(0.d0)*2
DOUBLE PRECISION,DIMENSION(IE*JE)	,INTENT(out) :: x,y
INTEGER,DIMENSION((IE-1)*(JE-1),4)	,INTENT(out) :: ICORN
INTEGER											 :: i,j,k,cj,ci
thx=thx*PI/180
thy=thy*PI/180
DO i=1,IE*JE
	x(i)=lx*cos(thx)*MOD(i-1,IE)/(IE-1)+ly*sin(thy)*floor((i-0.01)/IE)/(JE-1)
	y(i)=lx*sin(thx)*MOD(i-1,IE)/(IE-1)+ly*cos(thy)*floor((i-0.01)/IE)/(JE-1)
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
END SUBROUTINE gridgen
! ####################################################