! ####################################################
SUBROUTINE AssemEnergy(BPV,ICORN,nelem,IE,JE,x,y,density,viscosity,Tdiffusivity,AT)
! Rambod Mojgani
IMPLICIT										 NONE
DOUBLE PRECISION,DIMENSION(3*(IE*JE))		   ,INTENT(in)	 :: BPV ! Pressure Veloctiy of last converged time\step
INTEGER,DIMENSION((IE-1)*(JE-1),4)			   ,INTENT(in)	 :: ICORN
INTEGER										   ,INTENT(in)	 :: nelem,IE,JE
DOUBLE PRECISION,DIMENSION(IE*JE)			   ,INTENT(in)	 :: x,y
DOUBLE PRECISION							   ,INTENT(in)	 :: density,viscosity,Tdiffusivity
DOUBLE PRECISION,DIMENSION(IE*JE,2*IE+3)	   ,INTENT(out)	 :: AT
DOUBLE PRECISION,DIMENSION(4,4)								 :: cpu,cpv,dpu,cc,cd,chu,chv,chpu,chpv,dpv,cd
DOUBLE PRECISION,DIMENSION(4)								 :: ddu,ddv
DOUBLE PRECISION,DIMENSION(4)								 :: u,v,p
DOUBLE PRECISION,DIMENSION(4)								 :: uo,vo,po
DOUBLE PRECISION,DIMENSION(IE*JE)							 :: um,vm,pm
DOUBLE PRECISION,DIMENSION(IE*JE)							 :: uom,vom,pom
INTEGER														 :: ne,nw,sw,se,elem,i,j,k,jj,kk,ib,jb
DOUBLE PRECISION,Dimension(4)								 :: Xi,Yi
DOUBLE PRECISION,Dimension(4,4)								 :: dndx,dndy,N
DOUBLE PRECISION,Dimension(4)								 :: voli,volinp
DOUBLE PRECISION,Dimension(4)								 :: uip,vip
DOUBLE PRECISION,Dimension(4,2)								 :: dsx,dsy
DOUBLE PRECISION,Dimension(4)								 :: l2d,dj,vbar
DOUBLE PRECISION,DIMENSION(4,4)								 :: cup
DOUBLE PRECISION,DIMENSION(4)								 :: dsup
! Separating u,v,p just for ease of handeling the parameters
DO i=1,IE*JE,1
	pm(i)=BPV(3*i-2)
	um(i)=BPV(3*i-1)
	vm(i)=BPV(3*i)
END DO
DO i=1,IE*JE,1
	pom(i)=BPV(3*i-2)
	uom(i)=BPV(3*i-1)
	vom(i)=BPV(3*i)
END DO

DO elem=1,nelem,1

	ne=ICORN(elem,1)
	nw=ICORN(elem,2)
	sw=ICORN(elem,3)
	se=ICORN(elem,4)

	u=(/um(ne),um(nw),um(se),um(sw)/)
	v=(/vm(ne),vm(nw),vm(se),vm(sw)/)
	p=(/pm(ne),pm(nw),pm(se),pm(sw)/)

	uo=(/uom(ne),uom(nw),uom(se),uom(sw)/)
	vo=(/vom(ne),vom(nw),vom(se),vom(sw)/)
	po=(/pom(ne),pom(nw),pom(se),pom(sw)/)

	Xi=(/x(ne),x(nw),x(sw),x(se)/)
	Yi=(/y(ne),y(nw),y(sw),y(se)/)
	CALL BLN(Xi,Yi,N,dndx,dndy,voli,volinp)

	DO j=1,4,1
		uip(j)=DOT_PRODUCT(N(1:4,j),uo)
		vip(j)=DOT_PRODUCT(N(1:4,j),vo)
	END DO

	DO j=1,4,1
		vbar(j)=( uip(j)*uip(j) + vip(j)*vip(j) )** 5.d-1 ! Velocity Magnitude in Integral Point
	END DO
	!---
	CALL upwind(Xi,Yi,uip,vip,cup,dsup)
	CALL localcd(Xi,Yi,dndx,dndy,Tdiffusivity,cd)
	CALL derv(Xi,Yi,volinp,dsx,dsy,l2d)
	DO j=1,4,1
		dj(j)=density*vbar(j)/dsup(j)+viscosity/l2d(j)
	ENDDO
	CALL localdp(Xi,Yi,N,dndx,dndy,cup,u,v,p,uip,vip,vbar,density,density*Tdiffusivity,dsx,dsy,l2d,dsup,dj,cc,dpu,dpv,cpu,cpv)
	dpu=0.d0
	dpv=0.d0
	cpu=0.d0
	dpv=0.d0
	! Bulding the Assembly Matrices in Band form; note the assmbley matrix is the compressed band form
	DO j=1,4,1
		jj=ICORN(elem,j)
		DO k=1,4,1
			kk=ICORN(elem,k)
			CALL BandReady(jj,kk,IE+1,IE+1,ib,jb)
			AT(ib,jb)=AT(ib,jb)+cc(j,k)+cd(j,k) ! As no Pressure Term in Energy Equation
		END DO
		! B(jj)=B(jj) Nothing needs to be done on B
	END DO
END DO
END SUBROUTINE AssemEnergy
! ####################################################
