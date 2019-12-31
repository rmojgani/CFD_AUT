!#####################################
FUNCTION distance(x1,y1,x2,y2)
! Rambod Mojgani
! 90129045
! Find Distance Between two points, A & B
! A=(x1,y1) & B=(x2,y2)
IMPLICIT										 NONE
DOUBLE PRECISION								:: x1,x2,y1,y2,distance
	distance=(x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)
	distance=distance**5.d-1
ENDFUNCTION distance
!#####################################
FUNCTION	triarea(x,y)
! Rambod Mojgani
! 90129045
! Use Heron Formula to calculates a triangle area
IMPLICIT										NONE
DOUBLE PRECISION,DIMENSION(3)					:: x,y
DOUBLE PRECISION								:: triarea,distance,p,a,b,c
	a=distance(x(1),y(1),x(2),y(2))
	b=distance(x(2),y(2),x(3),y(3))
	c=distance(x(3),y(3),x(1),y(1))
	p=5.d-1*(a+b+c)
	triarea=DSQRT(p*(p-a)*(p-b)*(p-c))
ENDFUNCTION
!#####################################
SUBROUTINE	NormalVector(x1,y1,x2,y2,nx,ny)
! Rambod Mojgani
! 90129045
! Calculating if normal vector
! Point(1) is the first point in CCW Roatation
! This means on CCW moving on a line, its left hand is considered as INSIDE
IMPLICIT										NONE
DOUBLE PRECISION					,INTENT(IN) :: x1,y1,x2,y2
DOUBLE PRECISION					,INTENT(OUT):: nx,ny
DOUBLE PRECISION								:: ds
		nx= (y2	- y1)
		ny=-(x2 - x1)
		ds=DSQRT(nx*nx+ny*ny)
		nx=nx / ds
		ny=ny / ds
ENDSUBROUTINE NormalVector
!#####################################