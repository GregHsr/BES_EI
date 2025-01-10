program main
implicit none 
integer ,parameter ::mdim=3,nz = 1

integer ::nx,ny,ndim
!
!  PARAMETRES ICCG
!
real*8 ,dimension(:,:), allocatable :: coef 
real*8 ,dimension(:), allocatable :: rhs1,p_s,r_s,r2_s
real*8 ,dimension(:), allocatable :: q_s,s_s,x1
real*8 ,dimension(:,:), allocatable ::l_s
integer, dimension(:), allocatable :: jcoef
integer, dimension(:), allocatable :: ldiag

real*8 :: dx,dy, dt, nu, u_tmp, v_tmp, nrj_n, nrj_n1, conv, tf
real*8 :: zeta,time,Re_lu
real*8, dimension(:), allocatable :: xx
real*8, dimension(:), allocatable :: yy
real*8, dimension(:,:), allocatable :: rhs,pre,u_cent,v_cent,rot,div
real*8, dimension(:,:), allocatable :: u,v,u_dif,v_dif
real*8 pi,sum,premoy,pamoy
integer i,j,k,itmax,isto,istep,nstep
integer Nx_lu, Schema_lu, script
	
external ICCG2

script=0
if (script == 1) then
	call read_data(Re_lu, Nx_lu, Schema_lu)
else
	Re_lu = 100
	Nx_lu = 10
	Schema_lu = 1
end if

nx = Nx_lu
ny = Nx_lu
ndim = nx*ny

! Allouer les tableaux
allocate(u(0:nx+1,0:ny+1),v(0:nx+1,0:ny+1),u_dif(0:nx+1,0:ny+1),v_dif(0:nx+1,0:ny+1))
allocate(u_cent(1:nx,1:ny),v_cent(1:nx,1:ny),rot(1:nx,1:ny),div(1:nx,1:ny))
allocate(rhs(1:nx,1:ny),pre(1:nx,1:ny))
allocate(xx(1:nx),yy(1:ny))

allocate(coef(1:ndim,1:mdim),rhs1(1:ndim),p_s(1:ndim),r_s(1:ndim),r2_s(1:ndim))
allocate(q_s(1:ndim),s_s(1:ndim),x1(1:ndim),l_s(1:ndim,1:5),jcoef(1:mdim),ldiag(1:5))

pi=4.*atan(1.)
time=0.
tf=0.
nrj_n=0.
nrj_n1=0.
zeta=1.e-8
itmax=300
	
dx=1./float(nx)
dy=1./float(ny)
	
do i=1,nx
	xx(i)=(i-0.5)*dx
enddo
	
do j=1,ny
	yy(j)=(j-0.5)*dy
enddo

!  Calcul des vitesses initiales
call initialize_un(u,v,nx,ny)

! Initilisation de l'algorithme 
istep=0
isto=1000
nstep=40000
dt = 0.
write(*,*) "Re = ",Re_lu
nu = 1/Re_lu
write(*,*) "nu = ",nu

u_cent=0.d0
v_cent=0.d0   
rot=0.
div=0.

do istep=0,nstep
	!   TIMESTEP
	call timestep(u,v,nx,ny,dx,dy,dt,nu)
	! write(*,*) "dt = ",dt
	time=time+dt

	!	ADVECTION
	call diffusion(u,v,u_dif,v_dif,dx,dy,dt,nu,nx,ny)

	call vitesse_upwind(u,v,u_dif,v_dif,nx,ny,dt,dx,dy)

	call boundaries(u,v,nx,ny)

	! Calcul du second membre de l'equation de pression
	do i=1,nx
       	do j=1,ny
        		rhs(i,j) = ((u(i,j) - u(i-1,j))/dx + (v(i,j) - v(i,j-1))/dy)/dt    
    	end do
   	end do
    	!write(*,'(4F10.2)') rhs

	!	GENERATION DE LA MATRICE des coef de l 'equation de pression

	call matgen_cavite(coef,jcoef,nx,ny,ndim,mdim,dx,dy) 

	sum=0.0
	do j=1,ny
	   do i=1,nx
	      sum=sum+rhs(i,j)
	   enddo
	enddo

	sum=sum/float(nx*ny)
	k=1
	do j=1,ny
	   do i=1,nx
	      rhs(i,j)=rhs(i,j)-sum
	      rhs1(k) =-rhs(i,j)
	      x1(k)=0.
	      k=k+1
	    enddo
	enddo

	!	RESOLUTION PRESSION

	call ICCG2(coef,jcoef,l_s,Ldiag,rhs1,x1, &
	  	       ndim,mdim,zeta,p_s,r_s,r2_s,q_s,s_s,itmax)
	k=1
	do j=1,ny 
	   do i=1,nx 
	      pre(i,j)=x1(k)
	      k=k+1
	    enddo
	enddo

	do i=1,nx-1
		do j=1,ny
			u_tmp = u(i,j) - dt*(pre(i+1,j)-pre(i,j))/dx
			u(i,j) = u_tmp
		end do 
	end do
	do i=1,nx
		do j=1,ny-1
			v_tmp = v(i,j) - dt*(pre(i,j+1)-pre(i,j))/dy
			v(i,j) = v_tmp
		end do 
	end do

	
	do i=1,nx
		do j=1,ny
			u_cent(i,j)=(u(i-1,j)+u(i,j))/2
			v_cent(i,j)=(v(i,j-1)+v(i,j))/2
		end do
	end do
	
	! Vérification de la convergence à l'aide de l'énergie cinétique
	nrj_n1=0.
	do i=1,nx
		do j=1,ny
			nrj_n1=nrj_n1+0.5*(u(i,j)**2+v(i,j)**2)
		end do
	end do
	nrj_n1=nrj_n1/(nx*ny)

	conv=abs(nrj_n-nrj_n1)/dt
	nrj_n=nrj_n1

	write(*,*) "Convergence = ",conv

	if (conv.lt.1.e-4 .and. istep>500) then
		write(*,*) "Convergence atteinte en ",istep," iterations"
		write(*,*) "Time step stockage", isto
		exit
	end if

	!  Calcul du rotationnel et de la divergence
	call divergence(u,v,dx,dy,nx,ny,div)
	call rotationnel(u,v,dx,dy,nx,ny,rot)

	call write_result_ensight(xx,yy,u_cent,v_cent,rot,div,pre,nx,ny,nz,istep,isto,nstep)

end do

write(*,*) "Temps final = ",time

! Export u on ny/2 and v on nx/2
open(1,file='v.dat')
write(1,*) 'x v'
do i=1,nx
	write(1,*) xx(i),v_cent(i,ny/2)
end do
close(1)
open(2,file='u.dat')
write(2,*) ' y u'
do j=1,ny
	write(2,*) yy(j), u_cent(nx/2,j)
end do
close(2)

deallocate(u,v,u_dif,v_dif,u_cent,v_cent,rot,div,rhs,pre,xx,yy)
deallocate(coef,rhs1,p_s,r_s,r2_s,q_s,s_s,x1,l_s,jcoef,ldiag)

end 
