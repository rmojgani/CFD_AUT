! ####################################################
SUBROUTINE localch(density,dj,uip,vip,N,dndx,dndy,dsx,dsy,chu,chv,chp)
! Rambod Mojgani
! ------------------------ chu,chv,chp (Making Continunity Equations) - بقاي جرم
IMPLICIT										 NONE
DOUBLE PRECISION					,INTENT(in)	 :: density
DOUBLE PRECISION,Dimension(4)		,INTENT(in)	 :: dj
DOUBLE PRECISION,DIMENSION(4)		,INTENT(in)  :: uip,vip
DOUBLE PRECISION,Dimension(4,4)		,INTENT(in)	 :: N,dndx,dndy
DOUBLE PRECISION,Dimension(4,2)		,INTENT(in)	 :: dsx,dsy
DOUBLE PRECISION,Dimension(4,4)		,INTENT(out) :: chu,chv,chp
INTEGER,DIMENSION(4)							 :: jlist,klist
INTEGER											 :: i,j,k,jj

! To avoid "IF" in many loops.
jlist=(/1,2,3,4/)
klist=(/4,1,2,3/)
DO jj=1,4,1
	DO i=1,4,1
		! To ready
		j=jlist(jj)
		k=klist(jj)
		! Building chu
		chu(j,i)=			density*dsx(j,1) * (  N(j,i)+(density/dj(j)) * (-vip(j)*dndy(j,i))  )
		chu(j,i)=chu(j,i) + density*dsy(j,1) * (		 (density/dj(j)) * (+vip(j)*dndx(j,i))  )
		chu(j,i)=chu(j,i) + density*dsx(j,2) * (  N(k,i)+(density/dj(k)) * (-vip(k)*dndy(k,i))  )
		chu(j,i)=chu(j,i) + density*dsy(j,2) * (		 (density/dj(k)) * (+vip(k)*dndx(k,i))  )
		! Building chv
		chv(j,i)=			density*dsx(j,1) * (		 (density/dj(j)) * (-uip(j)*dndy(j,i))  )
		chv(j,i)=chv(j,i) + density*dsy(j,1) * (  N(j,i)+(density/dj(j)) * (+uip(j)*dndx(j,i))  )
		chv(j,i)=chv(j,i) + density*dsx(j,2) * (		 (density/dj(k)) * (-uip(k)*dndy(k,i))	)
		chv(j,i)=chv(j,i) + density*dsy(j,2) * (  N(k,i)+(density/dj(k)) * (+uip(k)*dndx(k,i))	)
		! Building chp
		chp(j,i)=		  +	(density/dj(j)) * dsx(j,1) * (-dndx(j,i))
		chp(j,i)=chp(j,i) +	(density/dj(j)) * dsy(j,1) * (-dndy(j,i))
		chp(j,i)=chp(j,i) +	(density/dj(k)) * dsx(j,2) * (-dndx(k,i))
		chp(j,i)=chp(j,i) +	(density/dj(k)) * dsy(j,2) * (-dndy(k,i))
	END DO
END DO

ENDSUBROUTINE localch