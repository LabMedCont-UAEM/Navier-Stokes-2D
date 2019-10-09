!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!A basic solver for 2D Navier-Stokes equations.
!The applied method is finite difference time-dependent (FDTD).
!By F. Garzon
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
program Navier2DPolar
implicit none

integer :: nr,nth,i,j,t,Tmax,imin,imax,jmin,jmax
real :: Rmin,Rmax,Lr,lth,dr,dth,pi,Re,Ga,dt,delta

real, dimension(:), allocatable :: r,th
real, dimension(:,:), allocatable :: u1,v1,u2,v2,F,G,&
				p1,p2,flr,flth

real :: dv2dr,dvudth,du2dth,duvdr,&
	d2vdr2,dvdr,d2vdth2,dudth,&
	d2udr2,dudr,d2udth2,dvdth

pi=ACOS(-1.0)

Rmin=0.01 ; Rmax=1.0 ; lth=2.0*pi
Lr=Rmax-Rmin
nr=60 ; nth=60
dr=Lr/float(nr) ; dth=lth/float(nth-1)

imin=1 ; jmin=1
imax=nr-1 ; jmax=nth-1

Ga=0.7
Re=10.0

dt=0.00001
Tmax=nint(1.0/dt)
print*, 'T=',Tmax

allocate(r(0:nr),th(0:nth))
allocate(u1(0:nr,0:nth),v1(0:nr,0:nth),&
	u2(0:nr,0:nth),v2(0:nr,0:nth),&
	p1(0:nr,0:nth),p2(0:nr,0:nth),&
	F(0:nr,0:nth),G(0:nr,0:nth),&
	flr(0:nr,0:nth),flth(0:nr,0:nth))

!-------- GRID ---------------
th(0)=0.0
th(1)=0.5*dth
DO i=(2),nth-1
th(i)=0.5*dth+(i-1)*dth
END DO
th(nth)=lth


r(0)=0.0 + Rmin
r(1)=0.5*dr+Rmin
DO j=2,nr-1
r(j)=Rmin+0.5*dr+(j-1)*dr
END DO
r(nr)=Rmax


u1=0.0 ; v1=0.0 ; p1=0.0
u2=u1 ; v2=v1 ; p2=p1

DO t=1,Tmax

CALL BOUNDCONDITIONS(u1,u2,v1,v2,nr,nth,&
			imin,imax,jmin,jmax,nth)	
DO i=imin,imax-1
DO j=jmin,jmax
	dv2dr=(0.25/dr)*((v1(i,j)+v1(i+1,j))**2 -(v1(i-1,j)+v1(i,j))**2)&
		+(0.25*Ga/dr)*((ABS(v1(i,j)+v1(i+1,j)))*(v1(i,j)-v1(i+1,j))&
		-(ABS(v1(i-1,j)+v1(i,j)))*(v1(i-1,j)-v1(i,j)))
	
	dvudth=(0.25/dth)*( (u1(i,j)+u1(i+1,j))*(v1(i,j)+v1(i,j+1)) &
		- (u1(i,j-1)+u1(i+1,j-1))*(v1(i,j-1)+v1(i,j)) ) &
		+(0.25*Ga/dth)* &
		( ABS(u1(i,j)+u1(i+1,j))*(v1(i,j)-v1(i,j+1)) &
		- ABS(u1(i,j-1)+u1(i+1,j-1))*(v1(i,j-1)-v1(i,j)) )

	d2vdr2=(v1(i+1,j)-2.0*v1(i,j)+v1(i-1,j))/(dr*dr)
	
	dvdr=(v1(i+1,j)-v1(i-1,j))/(2.0*dr)

	dudth=(u1(i,j+1)-u1(i,j-1))/(2.0*dth)

	d2vdth2=(v1(i,j+1)-2.0*v1(i,j)+v1(i,j-1))/(dth*dth)

	flr(i,j)=0.0

	F(i,j)=v1(i,j)+dt*( (1.0/Re)*(d2vdr2 + (1.0/r(i))*dvdr + (1.0/(r(i)*r(i)))*d2vdth2&
			-v1(i,j)/(r(i)*r(i)) -dudth*2.0/(r(i)*r(i)) )&
			+(u1(i,j)**2 - v1(i,j)**2)*(1./r(i)) -dv2dr -(1.0/r(i))*dvudth +flr(i,j))
ENDDO
ENDDO

DO i=imin,imax
DO j=jmin,jmax-1
	du2dth=(0.25/dth)*((u1(i,j)+u1(i,j+1))**2 &
		- (u1(i,j-1)+u1(i,j))**2) &
		+ (0.25*Ga/dth)* &
		( ABS(u1(i,j)+u1(i,j+1))*(u1(i,j)-u1(i,j+1)) &
		- ABS(u1(i,j-1)+u1(i,j))*(u1(i,j-1)-u1(i,j)) )
	
	duvdr=(0.25/dr)*((v1(i,j)+v1(i,j+1))*(u1(i,j)+u1(i+1,j))&
		 - (v1(i-1,j)+v1(i-1,j+1))*(u1(i-1,j)+u1(i,j)))&
		 + (0.25*Ga/dr)*( ABS(v1(i,j)+v1(i,j+1))*(u1(i,j)-u1(i+1,j)) &
		 - ABS(v1(i-1,j)+v1(i-1,j+1))*(u1(i-1,j)-u1(i,j)) )

	d2udr2=(u1(i+1,j)-2.0*u1(i,j)+u1(i-1,j))/(dr*dr)
	
	dudr=(u1(i+1,j)-u1(i-1,j))/(2.0*dr)
	
	dvdth=(v1(i,j+1)-v1(i,j-1))/(2.0*dth)
	
	d2udth2=(u1(i,j+1)-2.0*u1(i,j)+u1(i,j-1))/(dth*dth)

	flth(i,j)=0.0

	G(i,j)=u1(i,j)+dt*((1.0/Re)*(d2udr2 + (1.0/r(i))*dudr + (1.0/(r(i)*r(i)))*d2udth2 &
			-u1(i,j)/(r(i)*r(i))+dvdth*2.0/(r(i)*r(i)) )&
			-(2.0*u1(i,j)*v1(i,j))/r(i)-(1.0/r(i))*du2dth -duvdr + flth(i,j) )
ENDDO
ENDDO


FORALL(j=jmin:jmax)
	F(0,j)=v1(0,j) ; F(imax,j)=v1(imax,j)
END FORALL

FORALL(i=imin:imax)
	G(i,0)=u1(i,0) ; G(i,jmax)=u1(i,jmax)
END FORALL

CALL PRESSURE(p1,p2,F,G,r,nr,nth,imin,imax,jmin,jmax,dr,dth,dt)

FORALL(i=imin:imax-1,j=jmin:jmax)
	v2(i,j)=F(i,j)-(dt/dr)*(p2(i+1,j)-p2(i,j))
END FORALL

FORALL(i=imin:imax,j=jmin:jmax-1)
	u2(i,j)=G(i,j)-(1.0/r(i))*(dt/dth)*(p2(i,j+1)-p2(i,j))
END FORALL

if(mod(t,500).EQ.0) THEN
	delta=MAXVAL(ABS(u2-u1))
	write(*,*) 't=',t,'E=',delta,'v2=', v2(nr/2,nth/2),'u2=', u2(nr/2,nth/2)
end if
u1=u2
v1=v2
ENDDO

FORALL(i=imin:imax,j=jmin:jmax)
	u1(i,j)=0.5*(u2(i,j)+u2(i-1,j))
        v1(I,J)=0.5*(v2(i,j)+v2(i,j-1))
END FORALL

do i=0,nr
u1(i,nth)=(u1(i,0)+u1(i,nth-1))/2.0
v1(i,nth)=(v1(i,0)+v1(i,nth-1))/2.0
enddo

open(200,file='V_x_y.dat')
write(200,*) 'TITLE="NS2DPolar"'
write(200,*) 'VARIABLES=X,Y,U,V'
write(200,*) 'ZONE T="1", I=',imax,', J=',jmax+2
do j=jmin-1,jmax+1
do i=imin,imax
	write(200,*) r(i)*cos(th(j)), r(i)*sin(th(j)),&
			v1(i,j)*cos(th(j))-u1(i,j)*sin(th(j)),&
			v1(i,j)*sin(th(j))+u1(i,j)*cos(th(j))
enddo
enddo
close(200)

write(*,*) 'End of program'
end program

subroutine BOUNDCONDITIONS(u1,u2,v1,v2,nr,nth,&
			imin,imax,jmin,jmax,th)
implicit none

integer,value :: nr,nth,imin,jmin,imax,jmax
real, intent(inout) :: th(0:nth)
real, intent(inout) :: u1(0:nr,0:nth), v1(0:nr,0:nth),&
			u2(0:nr,0:nth), v2(0:nr,0:nth)
integer :: i,j
!r=imin-1 ,LEFT
do j=0,nth
v1(imin-1,j)=v1(imin,j)
u1(imin-1,j)=u1(imin,j)
enddo

!r=imax ,RIGHT
do j=0,nth
v1(imax+1,j)=0.0
u1(imax+1,j)=2.0*(1.0)-u1(imax,j)
enddo

!th=jmin-1 , BOT
do i=0,nr
u1(i,jmax)=u1(i,jmin)
v1(i,jmax)=v1(i,jmin)
enddo

!th=jmax+1, TOP
do i=0,nr
u1(i,jmin-1)=u1(i,jmax-1)
v1(i,jmin-1)=v1(i,jmax-1)
enddo

v2=v1
u2=u1
return
end subroutine

subroutine PRESSURE(p1,p2,F,G,r,nr,nth,imin,imax,jmin,jmax,dr,dth,dt) 
implicit none

integer, value :: nr,nth,imin,imax,jmin,jmax
real, value :: dr,dth,dt
real, value :: r(0:nr)
real, value :: p1(0:nr,0:nth),p2(0:nr,0:nth),F(0:nr,0:nth),G(0:nr,0:nth)

integer :: pt,pj,i,j
real :: pdt

pdt=0.000001
pt=100

do pj=1,pt
FORALL(i=imin:imax,j=jmin:jmax)
	p2(i,j)=p1(i,j)+pdt*( (p1(i+1,j)-2.*p1(i,j)+p1(i-1,j))/(dr*dr)&
		+(1.0/r(i))*(p1(i+1,j)-p1(i-1,j))/(2.0*dr)&
		+(1.0/(r(i)*r(i)))*(p1(i,j+1)-2.*p1(i,j)+p1(i,j-1))/(dth*dth)&
		-(1.0/dt)*( (F(i,j)-F(i-1,j))/(dr)+(1.0/r(i))*F(i,j)&
		+(1.0/r(i))*(G(i,j)-G(i,j-1))/(dth) ) )
END FORALL

do j=0,nth
p2(imax+1,j)=p2(imax,j)
p2(imin-1,j)=p2(imin,j)
enddo

do i=0,nr
p2(i,jmin-1)=p2(i,jmin)
p2(i,jmax+1)=p2(i,jmax)
p2(i,jmin)=p2(i,jmax)
enddo

p1=p2
enddo
return
end subroutine
