! ####################################################
SUBROUTINE solve2l(a,b,c,x,y)
! Rambod Mojgani
! Finds intersection of two lines; l1,l2
! l1 :: A(1) x + B(1) y = C(1)
! l2 :: A(2) x + B(2) y = C(2)
! NOTE :: IF one line is parallel to "Y-Axis", It does NOT work!
IMPLICIT NONE
DOUBLE PRECISION,DIMENSION(2),INTENT(in)    :: a,b,c
DOUBLE PRECISION			 ,INTENT(out)   :: x,y
DOUBLE PRECISION		    				:: delta,deltax,deltay
delta=a(1)*b(2)-a(2)*b(1)
deltax=b(2)*c(1)-b(1)*c(2)
deltay=a(1)*c(2)-c(1)*a(2)
x=deltax/delta
y=deltay/delta
END SUBROUTINE solve2l
! ####################################################