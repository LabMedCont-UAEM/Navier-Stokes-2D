!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!A basic solver for 2D Navier-Stokes equations.
!The applied method is finite difference time-dependent (FDTD).
!By F. Garzon
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
program navier2d
implicit none

integer :: i,j,t,nx,ny,Tmax,x1,x2,y1,y2,Nplot
real :: dx,dy,lx,ly,dt,Re,Ga,Delta,ti,tf
real, dimension(:), allocatable :: x,y
real, dimension(:,:), allocatable :: u1,u2,v1,v2,p1,p2,F,G

real, dimension(:,:), allocatable :: du2dx,duvdy,d2udx2,d2udy2
real, dimension(:,:), allocatable :: dv2dy,dvudx,d2vdx2,d2vdy2

real, dimension(:,:), allocatable :: flx,fly
call CPU_time(ti)
Nplot=100 !Step to print data arrays

!Setting domain parameters for the study of 2D Poiseuille's flow
lx=6.0 ; ly=1.0 !Domain length
nx=256 ; ny=64 !Cells in x and y direction
dx=lx/float(nx-1) ; dy=ly/float(ny-1) !Spatial step

x1=1 ; x2=nx-1 ; y1=1 ; y2=ny-1

Re=1.0 !Reynolds number
Ga=0.7 !Gamma parameter

allocate(x(0:nx),y(0:ny))
allocate(u1(0:nx,0:ny),u2(0:nx,0:ny),v1(0:nx,0:ny),v2(0:nx,0:ny))
allocate(F(0:nx,0:ny),G(0:nx,0:ny))
allocate(p1(0:nx,0:ny),p2(0:nx,0:ny))

allocate(du2dx(0:nx,0:ny),duvdy(0:nx,0:ny),d2udx2(0:nx,0:ny),d2udy2(0:nx,0:ny))
allocate(dv2dy(0:nx,0:ny),dvudx(0:nx,0:ny),d2vdx2(0:nx,0:ny),d2vdy2(0:nx,0:ny))

allocate(flx(0:nx,0:ny),fly(0:nx,0:ny))

call MESH(x,dx,nx,lx)
call MESH(y,dy,ny,ly)

dt=0.000001
Tmax=50000

u1=0.0 ; v1=0.0 ; p1=0.0
u2=u1 ; v2=v1 ; p2=p1

!=====================
!=====================
!=====================

DO t=1,Tmax

CALL BOUNDCOND(u1,u2,v1,v2,nx,ny)

FORALL(i=x1:x2-1,j=y1:y2)
	du2dx(i,j)=(0.25/dx)*( (u1(i,j)+u1(i+1,j))**2 &
		- (u1(i-1,j)+u1(i,j))**2) &
		+(0.25*Ga/dx)* &
		( ABS(u1(i,j)+u1(i+1,j))*(u1(i,j)-u1(i+1,j))&
		 - ABS(u1(i-1,j)+u1(i,j))*(u1(i-1,j)-u1(i,j)) )

	duvdy(i,j)=(0.25/dy)*( (v1(i,j)+v1(i+1,j))*(u1(i,j)+u1(i,j+1)) &
		- (v1(i,j-1)+v1(i+1,j-1))*(u1(i,j-1)+u1(i,j)) ) &
		+(0.25*Ga/dy)* &
		( ABS(v1(i,j)+v1(i+1,j))*(u1(i,j)-u1(i,j+1)) &
		- ABS(v1(i,j-1)+v1(i+1,j-1))*(u1(i,j-1)-u1(i,j)) )

	d2udx2(i,j)=(1.0/(dx*dx))*(u1(i+1,j)-2.0*u1(i,j)+u1(i-1,j))

	d2udy2(i,j)=(1.0/(dy*dy))*(u1(i,j+1)-2.0*u1(i,j)+u1(i,j-1))

	flx(i,j)=0.0 !External force in x direction

	F(i,j)=u1(i,j)+dt*( (1.0/Re)*(d2udx2(i,j) + d2udy2(i,j))&
		 - du2dx(i,j) - duvdy(i,j) &
		+ flx(i,j) )
END FORALL

FORALL(i=x1:x2,j=y1:y2-1)
	dv2dy(i,j)=(0.25/dy)*((v1(i,j)+v1(i,j+1))**2 &
		- (v1(i,j-1)+v1(i,j))**2) &
		+ (0.25*Ga/dy)* &
		( ABS(v1(i,j)+v1(i,j+1))*(v1(i,j)-v1(i,j+1)) &
		- ABS(v1(i,j-1)+v1(i,j))*(v1(i,j-1)-v1(i,j)) )

	dvudx(i,j)=(0.25/dx)*((u1(i,j)+u1(i,j+1))*(v1(i,j)+v1(i+1,j))&
		 - (u1(i-1,j)+u1(i-1,j+1))*(v1(i-1,j)+v1(i,j)))&
		 + (0.25*Ga/dx)*( ABS(u1(i,j)+u1(i,j+1))*(v1(i,j)-v1(i+1,j)) &
		 - ABS(u1(i-1,j)+u1(i-1,j+1))*(v1(i-1,j)-v1(i,j)) )

	d2vdx2(i,j)=(1.0/(dx*dx))*(v1(i+1,j)-2.0*v1(i,j)+v1(i-1,j))

	d2vdy2(i,j)=(1.0/(dy*dy))*(v1(i,j+1)-2.0*v1(i,j)+v1(i,j-1))

	fly(i,j)=0.0 !External force in y direction

	G(i,j)=v1(i,j)+dt*( (1.0/Re)*(d2vdx2(i,j)+d2vdy2(i,j))&
		 -dvudx(i,j)-dv2dy(i,j)&
		+fly(i,j) )
END FORALL

FORALL(j=y1:y2)
	F(0,j)=u1(0,j) ; F(x2,j)=u1(x2,j)
END FORALL

FORALL(i=x1:x2)
	G(i,0)=v1(i,0) ; G(i,y2)=v1(i,y2)
END FORALL

CALL PRESSURE(p1,p2,F,G,nx,ny,dx,dy,dt)

FORALL(i=x1:x2-1,j=y1:y2)
	u2(i,j)=F(i,j)-(dt/dx)*(p2(i+1,j)-p2(i,j))
END FORALL

FORALL(i=x1:x2,j=y1:y2-1)
	v2(i,j)=G(i,j)-(dt/dy)*(p2(i,j+1)-p2(i,j))
END FORALL

IF(mod(t,500)==0) THEN
	Delta=MAXVAL(ABS(u2-u1))
	write(*,*) 'T=',t,'U1=', u1(nx/2,ny/2),'V1=', v1(nx/2,ny/2), 'E=', Delta
END IF

IF(mod(t,Nplot)==0) THEN
	CALL ANIMACION(x,y,u1,v1,nx,ny,Nplot,t)
END IF

u1=u2
v1=v2
ENDDO

!=====================
!=====================
!=====================

!Finding values at the center of the numerical cell
FORALL(i=x1:x2,j=y1:y2)
	u1(i,j)=0.5*(u1(i,j)+u2(i-1,j))
        v1(I,J)=0.5*(v2(i,j)+v2(i,j-1))
END FORALL

OPEN(1,file='data.dat')

!Saving inner data
do i=x1, x2
do j=y1, y2
	write(1,*) x(i),y(j),u1(i,j),v1(i,j)
	if(i==x2/2) then
	write(2,*) x(i),y(j),u1(i,j),v1(i,j)
	end if
enddo
enddo
CLOSE(1)
call CPU_time(tf)
write(*,*) 'Time=',tf-ti 
end program


!========================
!========================
!========================
subroutine MESH(vec,step,n,LL)
implicit none
integer, value :: n
real, value :: step, LL
real, intent(inout) :: vec(0:n)

integer :: m

vec(0)=0.0
vec(1)=step*0.5
do m=2,n-1
	vec(m)=0.5*step+(m-1)*step
enddo
vec(n)=LL

return
end subroutine

subroutine BOUNDCOND(u1,u2,v1,v2,nx,ny)
implicit none

integer, value :: nx,ny
real, intent(inout) :: u1(0:nx,0:ny),u2(0:nx,0:ny),v1(0:nx,0:ny),v2(0:nx,0:ny)
integer :: i,j

!x=0,y - LEFT
!Inflow
do j=0,ny
u1(0,j)=1.0
v1(0,j)=-v1(1,j)
enddo

!x=nx,y - RIGHT
!Outflow
do j=0,ny
u1(nx-1,j)=u1(nx-2,j)
v1(nx,j)=v1(nx-1,j)
enddo


!x,y=ny - TOP
!No-slip
do i=0,nx
u1(i,ny)=-u1(i,ny-1)
v1(i,ny-1)=0.0 !No-Slip
enddo


!x,y=0 - BOT
!No-slip
do i=0,nx
u1(i,0)=-u1(i,1)
v1(i,0)=0.0
enddo

u2=u1
v2=v1
end subroutine

subroutine PRESSURE(p1,p2,F,G,nx,ny,dx,dy,dt)
implicit none

integer, value :: nx,ny
real, value :: dx,dy,dt
real, intent(inout) :: p1(0:nx,0:ny),p2(0:nx,0:ny),F(0:nx,0:ny),G(0:nx,0:ny)
	
integer :: pt,pj,i,j
real :: pdt

pdt=0.000001
pt=100

DO pj=1, pt

FORALL(i=1:nx-1,j=1:ny-1)
	p2(i,j)=p1(i,j)+pdt*( (p1(i+1,j)-2.0*p1(i,j)+p1(i-1,j))/(dx*dx) &
		+ (p1(i,j+1)-2.0*p1(i,j)+p1(i,j-1))/(dy*dy)&
		- (1.0/dt)*( (F(i,j)-F(i-1,j))/(dx) + (G(i,j)-G(i,j-1))/(dy) ) )
END FORALL

!Neumann boundary conditions for pressure
do j=0,ny
p2(0,j)=p2(1,j)
p2(nx,j)=p2(nx-1,j)
enddo

do i=0,nx
p2(i,0)=p2(i,1)
p2(i,ny)=p2(i,ny-1)
enddo

p1=p2
ENDDO
return
end subroutine

subroutine ANIMACION(x,y,u1,v1,nx,ny,Nplot,t)
implicit none

integer, value :: nx,ny, Nplot, t
real :: x(0:nx),y(0:ny),u1(0:nx,0:ny),v1(0:nx,0:ny)
integer :: i,j

integer :: NCOUNT
character :: EXT*4, DESTINY*512, NCTEXT*4, OUTFILE*512

EXT='.dat'
DESTINY='LOCAL_PATH/anim/' !Add your local anim folder path
NCOUNT=t/Nplot
WRITE(NCTEXT,'(I0)') NCOUNT

OUTFILE=trim(DESTINY)//trim(NCTEXT)//EXT

OPEN(10,file=OUTFILE)

do i=1,nx-1
do j=1,ny-1
write(10,*) x(i),y(j),u1(i,j),v1(i,j)
enddo
enddo
return
end subroutine
