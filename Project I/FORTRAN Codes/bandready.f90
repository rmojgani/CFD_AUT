! ####################################################
SUBROUTINE BandReady(i,j,N1,N2,ii,jj)
! Rambod Mojgani 91.02.05; ED 1.d0
IMPLICIT											NONE
INTEGER,INTENT(in)									:: i,j,N1,N2
INTEGER,INTENT(out)									:: ii,jj
INTEGER												:: deviation
deviation=j-i
jj=N1+1+deviation
ii=i
END SUBROUTINE BandReady
! ####################################################