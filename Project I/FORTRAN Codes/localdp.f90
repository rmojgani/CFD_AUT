SUBROUTINE localdp(Xi,Yi,N,dndx,dndy,cup,u,v,p,uip,vip,vbar,density,gamma,dsx,dsy,l2d,dsup,dj,cc,dpu,dpv,cpu,cpv)
! Rambod Mojgani
IMPLICIT										NONE
DOUBLE PRECISION,DIMENSION(4)		,INTENT(in)  :: Xi,Yi			! Corners Coordinates
DOUBLE PRECISION,Dimension(4,4)		,INTENT(in)	 :: N,dndx,dndy,cup
DOUBLE PRECISION,DIMENSION(4)		,INTENT(in)	 :: u,v,p			! In 4 conrners
DOUBLE PRECISION,DIMENSION(4)		,INTENT(in)  :: uip,vip,vbar			! Velocity Componnets on Integral Points (j)
DOUBLE PRECISION					,INTENT(in)	 :: density,gamma
DOUBLE PRECISION,Dimension(4,2)		,INTENT(in)	 :: dsx,dsy
DOUBLE PRECISION,Dimension(4)		,INTENT(in)	 :: l2d,dsup,dj
DOUBLE PRECISION,Dimension(4,4)		,INTENT(out) :: cc,dpu,dpv,cpu,cpv
DOUBLE PRECISION,Dimension(4)					 :: uhat,vhat
DOUBLE PRECISION,Dimension(4)					 :: a,b,c,mdot
INTEGER,DIMENSION(4)							 :: jlist,klist
DOUBLE PRECISION								 :: dpdx,dpdy,dudx,dudy,dvdx,dvdy ! once made, used, erased
INTEGER											 :: i,j,jj,k
jlist=(/1,2,3,4/)
klist=(/4,1,2,3/)

DO j=1,4,1
	dpdx=DOT_PRODUCT(dndx(j,1:4),p)
	dpdy=DOT_PRODUCT(dndy(j,1:4),p)
	dudx=DOT_PRODUCT(dndx(j,1:4),u)
	dudy=DOT_PRODUCT(dndy(j,1:4),u)
	dvdx=DOT_PRODUCT(dndx(j,1:4),v)
	dvdy=DOT_PRODUCT(dndy(j,1:4),v)
	uhat(j)=uip(j)+(-dpdx+density*(uip(j)*dvdy-vip(j)*dudy))/dj(j)
	vhat(j)=vip(j)+(-dpdy+density*(vip(j)*dudx-uip(j)*dvdx))/dj(j)
END DO

DO j=1,4,1
	mdot(j) = density * (uhat(j)*dsx(j,1)+vhat(j)*dsy(j,1))
	a(j) = density * vbar(j)/dsup(j)
	b(j) = gamma / l2d(j)
	c(j) = 1.d0 / ( a(j) + b(j) )
END DO
DO jj=1,4,1
	DO	i=1,4,1
		! To ready
		j=jlist(jj)
		k=klist(jj)
		! u & v Term
		cc(j,i) =			mdot(j) * c(j) * (  a(j) * cup(i,j) + b(j) * N(j,i)  )
		cc(j,i) = cc(j,i) - mdot(k) * c(k) * (  a(k) * cup(i,k) + b(k) * N(k,i)  )

		! Pressure Term
		dpu(j,i) = N(j,i)*dsx(j,1) + N(k,i)*dsx(j,2)
		dpv(j,i) = N(j,i)*dsy(j,1) + N(k,i)*dsy(j,2)
	
		cpu(j,i) = -mdot(j) * c(j) * dndx(j,i)	+	mdot(k)* c(k) * dndx(k,i)
		cpv(j,i) = -mdot(j) * c(j) * dndy(j,i)	+	mdot(k)* c(k) * dndy(k,i)
	END DO
END DO
ENDSUBROUTINE localdp