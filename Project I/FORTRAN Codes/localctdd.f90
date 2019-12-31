SUBROUTINE localctdd(density,dt,voli,uo,vo,ct,ddu,ddv)
! Rambod Mojgani
IMPLICIT											NONE
DOUBLE PRECISION					,INTENT(in)	 :: density,dt
DOUBLE PRECISION,Dimension(4)		,INTENT(in)	 :: voli
DOUBLE PRECISION,DIMENSION(4)		,INTENT(in)	 :: uo,vo
DOUBLE PRECISION,Dimension(4)		,INTENT(out) :: ct,ddu,ddv
INTEGER											 :: i
ct=0.d0
ddu=0.d0
ddv=0.d0
DO i=1,4,1
	ct(i) = density * voli(i) / dt
	ddu(i) = -density * voli(i) / dt * uo(i)
	ddv(i) = -density * voli(i) / dt * vo(i)
END DO
ENDSUBROUTINE localctdd