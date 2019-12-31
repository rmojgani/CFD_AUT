! ####################################################
SUBROUTINE BLN(Xi,Yi,N,dndx,dndy,JP,Jinp)
! A Bilinear Interpolation Derivative Subroutine
! Takes four corners coordinate (Xi,Yi) and the desired s,t local in
! coordinate and gives the Interpolation Weight Factors at corneres (N)
! and the N derivatives in Global Coordinatre System, also returns J volume 
! in for sub-domain point of cell (Jp) and also J volume for integral point (Jinp)
! ##---------------------------------------##
! ##	Accompanied SUBROUTINS ::		   ##
! ##     SUBROUTINE Vol(Xi,Yi,s,t,...	   ##
! ##   ,dxds,dxdt,dyds,dydt,J,N,dndt,dnds) ##
! ##---------------------------------------##
! Rambod Mojgani
! Edition 3.0 - 91.02.2
IMPLICIT											NONE
DOUBLE PRECISION,DIMENSION(4),INTENT(in)			:: Xi,Yi
DOUBLE PRECISION,DIMENSION(1:4,1:4),INTENT(out)		:: N
DOUBLE PRECISION,DIMENSION(1:4,1:4),INTENT(out)		:: dndx,dndy ! In integral points
DOUBLE PRECISION,DIMENSION(4),INTENT(out)			:: JP,Jinp
INTEGER												:: i,ii
DOUBLE PRECISION								    :: dxds,dxdt,dyds,dydt,J
DOUBLE PRECISION,DIMENSION(4)						:: dnds,dndt
DOUBLE PRECISION,DIMENSION(4)						:: slist,tlist,slistintp,tlistintp
DOUBLE PRECISION:: s,t
slistintp=(/+0.d0,-5.d-1,+0.d0,+5.d-1/)
tlistintp=(/+5.d-1,+0.d0,-5.d-1,+0.d0/)
slist=(/+5.d-1,-5.d-1,-5.d-1,+5.d-1/)
tlist=(/+5.d-1,+5.d-1,-5.d-1,-5.d-1/)
DO i=1,4,1
	s=slistintp(i)
	t=tlistintp(i)
	N(i,1:4)=(/ (1+s)*(1+t),(1-s)*(1+t),(1-s)*(1-t),(1+s)*(1-t) /)*2.5d-1
	CALL Vol(Xi,Yi,s,t,dxds,dxdt,dyds,dydt,J,N,dndt,dnds)
	dndx(i,1:4)=(1.d0/J)*(+dydt*dnds-dyds*dndt)
	dndy(i,1:4)=(1.d0/J)*(-dxdt*dnds+dxds*dndt)
END DO
! Volume on centers of incomplete control volume (JP)
DO i=1,4,1
	CALL Vol(Xi,Yi,slist(i),tlist(i),dxds,dxdt,dyds,dydt,JP(i),N(i,1:4),dndt,dnds)
END DO
! Volume on Integral Points (Jinp)
DO i=1,4,1
	CALL Vol(Xi,Yi,slistintp(i),tlistintp(i),dxds,dxdt,dyds,dydt,Jinp(i),N(i,1:4),dndt,dnds)
END DO
END SUBROUTINE BLN
! ####################################################
SUBROUTINE Vol(Xi,Yi,s,t,dxds,dxdt,dyds,dydt,J,N,dndt,dnds)
! A Bilinear Interpolation Volume Subroutine
! Is used in "BLN" Subroutine
! Takes the 4 corners coordinates of Xi & Yi in Global Coordinate System
! also takes the Interpolation Weight Factors at corneres (N)
! and returns the needed derivatives and the J volume
! Rambod Mojgani
DOUBLE PRECISION,DIMENSION(4),INTENT(in)						:: Xi,Yi
DOUBLE PRECISION,			  INTENT(in)						:: s,t
DOUBLE PRECISION,			  INTENT(out)						:: dxds,dxdt,dyds,dydt,J
DOUBLE PRECISION,DIMENSION(4),INTENT(out)						:: dndt,dnds
DOUBLE PRECISION,DIMENSION(4),INTENT(in)						:: N
dnds=(/ +(1+t)*2.5d-1, -(1+t)*2.5d-1, -(1-t)*2.5d-1, +(1-t)*2.5d-1 /)
dndt=(/ +(1+s)*2.5d-1, +(1-s)*2.5d-1, -(1-s)*2.5d-1, -(1+s)*2.5d-1 /)
dxds=dot_product(dnds,Xi)
dxdt=dot_product(dndt,Xi)
dyds=dot_product(dnds,Yi)
dydt=dot_product(dndt,Yi)
J=dxds*dydt-dyds*dxdt
END SUBROUTINE Vol
! ####################################################