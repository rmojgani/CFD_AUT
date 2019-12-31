!#####################################
FUNCTION distance(x1,y1,x2,y2)
! Rambod Mojgani
! Finds Distance Between two points, A & B
! A=(x1,y1) & B=(x2,y2)
IMPLICIT NONE
DOUBLE PRECISION :: x1,x2,y1,y2,distance
distance=(x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)
distance=distance**5.d-1
ENDFUNCTION distance
! ####################################################