!#####################################
FUNCTION	DEG2RAD(a)
! Rambod Mojgani
! 90129045
! Convert Angle "a" from Deg. to Radians.
IMPLICIT										NONE
DOUBLE PRECISION, PARAMETER						:: PI=ACOS(0.d0)*2
DOUBLE PRECISION								:: a,DEG2RAD
	DEG2RAD=a*PI/180.d0
ENDFUNCTION
!#####################################