program main
implicit none 
integer ,parameter ::nx=16,ny=16,nz=1
!
!  PARAMETRES ICCG
!
integer ,parameter ::ndim = nx*ny,mdim=3
real*8 ,dimension(1:ndim,1:mdim) :: coef 
real*8 ,dimension(1:ndim) :: rhs1,p_s,r_s,r2_s
real*8 ,dimension(1:ndim) :: q_s,s_s,x1
real*8 ,dimension(1:ndim,1:5)::l_s
integer, dimension(1:mdim):: jcoef
integer, dimension(1:5):: ldiag

real*8 :: dx,dy, dt, nu, u_tmp, v_tmp
real*8 :: zeta,time
real*8, dimension(1:nx) :: xx
real*8, dimension(1:ny) :: yy
real*8, dimension(1:nx,1:ny) :: rhs,pre,u_cent,v_cent,rot,div
real*8, dimension(0:nx+1,0:ny+1) :: u,v
real*8 pi,sum,premoy,pamoy
integer i,j,k,itmax,isto,istep,nstep
	
external ICCG2
	
pi=4.*atan(1.)
time=0.
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
call initialize(u,v,nx,ny)

! Initilisation de l'algorithme 
istep=2
isto=1
nstep=20
dt = 0.1
nu = 0.01

u_cent=0.d0
v_cent=0.d0   
rot=0.
div=0.

do istep=0,nstep
	!   TIMESTEP
	call timestep(u,v,nx,ny,dx,dy,dt, nu)
	write(*,*) "dt = ",dt

	!	ADVECTION
	call vitesse_tilt(u,v,u_cent,v_cent,nx,ny)
	call vitesse_upwind(u,v,u_cent,v_cent,nx,ny,dt,dx,dy)

	!        Calcul du second membre de l'equation de pression
	call calcul_rhs(u,v,rhs,dx,dy,dt,nx,ny)

	! 	  Calcul du rotationnel et de la divergence
	call divergence(u,v,dx,dy,nx,ny,div)
	call rotationnel(u,v,dx,dy,nx,ny,rot)

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
			u_tmp = u(i,j)/dt - pre(i,j)
			u(i,j) = u_tmp
		end do 
	end do
	do i=1,nx
		do j=1,ny-1
			v_tmp = v(i,j)/dt - pre(i,j)
			v(i,j) = v_tmp
		end do 
	end do

	call vitesse_tilt(u,v,u_cent,v_cent,nx,ny)
	call write_result_ensight(xx,yy,u_cent,v_cent,rot,div,pre,nx,ny,nz,istep,isto,nstep)

end do

end 
