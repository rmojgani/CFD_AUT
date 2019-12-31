! ####################################################
SUBROUTINE assemb(Blagged,BOld,ICORN,nelem,IE,JE,x,y,dt,density,viscosity,A,B)
! Rambod Mojgani
IMPLICIT										 NONE
DOUBLE PRECISION,DIMENSION(3*IE*JE)			   ,INTENT(in)	 :: Blagged,BOld
INTEGER,DIMENSION((IE-1)*(JE-1),4)			   ,INTENT(in)	 :: ICORN
INTEGER										   ,INTENT(in)	 :: nelem,IE,JE
DOUBLE PRECISION,DIMENSION(IE*JE)			   ,INTENT(in)	 :: x,y
DOUBLE PRECISION							   ,INTENT(in)	 :: dt,density,viscosity
DOUBLE PRECISION,DIMENSION(3*(IE*JE),3*(IE*JE)),INTENT(out)	 :: A
DOUBLE PRECISION,DIMENSION(3*(IE*JE))		   ,INTENT(out)	 :: B
DOUBLE PRECISION,DIMENSION(4,4)								 :: cpu,cpv,dpu,cc,cd,chu,chv,chp,dpv,cd
DOUBLE PRECISION,DIMENSION(4)								 :: ddu,ddv,ct
DOUBLE PRECISION,DIMENSION(4)								 :: u,v,p
DOUBLE PRECISION,DIMENSION(4)								 :: uo,vo,po
DOUBLE PRECISION,DIMENSION(IE*JE)							 :: um,vm,pm
DOUBLE PRECISION,DIMENSION(IE*JE)							 :: uom,vom,pom
INTEGER														 :: ne,nw,sw,se,elem,i,j,k,ii,jj,kk,ib,jb
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
	pm(i)=Blagged(3*i-2)
	um(i)=Blagged(3*i-1)
	vm(i)=Blagged(3*i)
END DO
DO i=1,IE*JE,1
	pom(i)=BOld(3*i-2)
	uom(i)=BOld(3*i-1)
	vom(i)=BOld(3*i)
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
	CALL localctdd(density,dt,voli,uo,vo,ct,ddu,ddv)
	CALL localcd(Xi,Yi,dndx,dndy,viscosity,cd)
	CALL derv(Xi,Yi,volinp,dsx,dsy,l2d)
	DO j=1,4,1
		dj(j)=density*vbar(j)/dsup(j)+viscosity/l2d(j)
	ENDDO
	CALL localdp(Xi,Yi,N,dndx,dndy,cup,u,v,p,uip,vip,vbar,density,viscosity,dsx,dsy,l2d,dsup,dj,cc,dpu,dpv,cpu,cpv)
	CALL localch(density,dj,uip,vip,N,dndx,dndy,dsx,dsy,chu,chv,chp)

	! Buldiing the Assembly Matrices in Band form; note the assmbley matrix is the compressed band form
	DO j=1,4,1
		jj=3*(ICORN(elem,j)-1)
		DO k=1,4,1
			kk=3*(ICORN(elem,k)-1)

			! as I Derived !
			CALL BandReady(jj+1,kk+1,3*IE+5,3*IE+5,ib,jb)
			A(ib,jb)=A(ib,jb)+chp(j,k)
			
			CALL BandReady(jj+1,kk+2,3*IE+5,3*IE+5,ib,jb)
			A(ib,jb)=A(ib,jb)+chu(j,k)

			CALL BandReady(jj+1,kk+3,3*IE+5,3*IE+5,ib,jb)
			A(ib,jb)=A(ib,jb)+chv(j,k)

			! as in Dr. Karimian Given Note
			CALL BandReady(jj+2,kk+1,3*IE+5,3*IE+5,ib,jb)
			A(ib,jb)=A(ib,jb)+dpu(j,k)+cpu(j,k)
			
			CALL BandReady(jj+2,kk+2,3*IE+5,3*IE+5,ib,jb)
			A(ib,jb)=A(ib,jb)+cd(j,k)+cc(j,k)

			CALL BandReady(jj+3,kk+1,3*IE+5,3*IE+5,ib,jb)
			A(ib,jb)=A(ib,jb)+dpv(j,k)+cpv(j,k)

			CALL BandReady(jj+3,kk+3,3*IE+5,3*IE+5,ib,jb)
			A(ib,jb)=A(ib,jb)+cd(j,k)+cc(j,k)

		END DO
		
		! as in Dr. Karimian Given Note
		CALL BandReady(jj+2,jj+2,3*IE+5,3*IE+5,ib,jb)
		A(ib,jb)=A(ib,jb)+ct(j)
		CALL BandReady(jj+3,jj+3,3*IE+5,3*IE+5,ib,jb)
		A(ib,jb)=A(ib,jb)+ct(j)

		B(jj+2)=B(jj+2)+ddu(j)
		B(jj+3)=B(jj+3)+ddv(j)
	END DO
END DO
END SUBROUTINE assemb
! ####################################################
