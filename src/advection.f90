subroutine initialize(un,vn,nx,ny) 
    implicit none
    integer, intent(in) :: nx,ny
    real*8, dimension(0:nx+1,0:ny+1), intent(out) :: un,vn 
    integer :: i,j
    do i=0,nx+1
        do j=0,ny+1
            un(i,j)=0.0
            vn(i,j)=0.0
        end do
    end do
end subroutine initialize


subroutine initialize_un(un,vn,nx,ny)
    implicit none
    integer, intent(in) :: nx,ny
    real*8, dimension(0:nx+1,0:ny+1), intent(out) :: un,vn 
    integer :: i,j
    do i=0,nx+1
        do j=0,ny+1
            un(i,j)=0.0
            vn(i,j)=0.0
        end do
    end do
    ! /!\ ATTENTION CONDITIONS LIMITES A IMPLEMENTER
    do i=1,nx
   	 ! Condition Bas
    	un(i,0) = - un(i,1)
        vn(i,0) = 0
        
        ! Condition Haut
        un(i,ny+1) = 2 - un(i,ny)
        vn(i,ny) = 0
    end do        
                
    do j=1,ny
	    ! Condition Gauche
	    un(0,j) = 0
        vn(0,j) = - vn(1,j)
        
        ! Condition Droite
	    un(nx,j) = 0
        vn(nx+1,j) = - vn(nx,j)      
    end do
   
    do i=0,ny+1
    !write(*,'(6F10.2)') un(0:nx+1, i)
    end do
    
end subroutine initialize_un


subroutine boundaries(un,vn,nx,ny)
    implicit none
    integer, intent(in) :: nx,ny
    real*8, dimension(0:nx+1,0:ny+1), intent(inout) :: un,vn 
    integer :: i,j
    ! /!\ ATTENTION CONDITIONS LIMITES A IMPLEMENTER
    do i=1,nx
   	 ! Condition Bas
    	un(i,0) = - un(i,1)
        vn(i,0) = 0
        
        ! Condition Haut
        un(i,ny+1) = 2 - un(i,ny)
        vn(i,ny) = 0
    end do        
                
    do j=1,ny
	    ! Condition Gauche
	    un(0,j) = 0
        vn(0,j) = - vn(1,j)
        
        ! Condition Droite
	    un(nx,j) = 0
        vn(nx+1,j) = - vn(nx,j)      
    end do
   
    do i=0,ny+1
    !write(*,'(6F10.2)') un(0:nx+1, i)
    end do
    
end subroutine boundaries



subroutine timestep(un,vn,nx,ny,dx,dy,dt,nu)
    ! Schéma Upwind et Centré
    implicit none
    integer, intent(in) :: nx,ny
    real*8, dimension(0:nx+1,0:ny+1), intent(in) :: un,vn
    real*8, intent(in) :: dx,dy,nu
    real*8, intent(out) :: dt
    real*8 :: umax, vmax
    integer :: i,j

    umax = 0.0
    vmax = 0.0

    do i = 0, nx+1
        do j = 0, ny+1
            if (abs(un(i,j)) > umax) then
                umax = abs(un(i,j))
            end if
            if (abs(vn(i,j)) > vmax) then
                vmax = abs(vn(i,j))
            end if
        end do
    end do

    dt = 1.0 / (2.0 * (umax / dx + vmax / dy) + 2.0 * nu * (1.0 / (dx * dx) + 1.0 / (dy * dy)))
    
end subroutine timestep
            


subroutine vitesse_upwind(un,vn,u_dif,v_dif,nx,ny, dt, dx, dy)
    implicit none
    integer, intent(in) :: nx,ny
    real*8, dimension(0:nx+1,0:ny+1), intent(inout) :: un,vn
    real*8, dimension(0:nx+1,0:ny+1), intent(in) :: u_dif,v_dif
    real*8, intent(in) :: dt,dx,dy
    real*8, dimension(0:nx+1,0:ny+1) :: u_temp, v_temp
    real*8 :: u_west, u_est, u_nord, u_sud
    real*8 :: v_west, v_est, v_nord, v_sud
    real*8 :: u_tilt_nord, u_tilt_sud, u_tilt_west, u_tilt_est
    real*8 :: v_tilt_nord, v_tilt_sud, v_tilt_west, v_tilt_est
    integer :: i,j

    ! Calcul U
    do i=1,nx-1
        do j=1,ny
            u_tilt_est = 0.5*(un(i,j)+un(i+1,j))
            u_tilt_west = 0.5*(un(i-1,j)+un(i,j))
            v_tilt_nord = 0.5*(vn(i,j)+vn(i+1,j))
            v_tilt_sud = 0.5*(vn(i,j-1)+vn(i+1,j-1))
            
            
            if (u_tilt_est >= 0) then
                u_est = un(i,j)
            else
                u_est = un(i+1,j)
            end if

            if (u_tilt_west >= 0) then
                u_west = un(i-1,j)
            else
                u_west = un(i,j)
            end if

            if (v_tilt_nord >= 0) then
                u_nord = un(i,j)
            else
                u_nord = un(i,j+1)
            end if

            if (v_tilt_sud >= 0) then
                u_sud = un(i,j-1)
            else
                u_sud = un(i,j)
            end if
            
            ! Calcul des vitesses au temps n+1
            u_temp(i,j) = un(i,j) - dt/dx*(u_tilt_est*u_est - u_tilt_west*u_west)&
                         - dt/dy*(v_tilt_nord*u_nord - v_tilt_sud*u_sud)
                         
            un(i,j) = u_temp(i,j) + u_dif(i,j)
                         
        end do
    end do  
    

    ! Calcul V
    do i=1,nx
        do j=1,ny-1
            u_tilt_est = 0.5*(un(i,j)+un(i,j+1))
            u_tilt_west = 0.5*(un(i-1,j)+un(i-1,j+1))
            v_tilt_nord = 0.5*(vn(i,j)+vn(i,j+1))
            v_tilt_sud = 0.5*(vn(i,j-1)+vn(i,j))
            
            if (u_tilt_est >= 0) then
                v_est = vn(i,j)
            else
                v_est = vn(i+1,j)
            end if

            if (u_tilt_west >= 0) then
                v_west = vn(i-1,j)
            else
                v_west = vn(i,j)
            end if

            if (v_tilt_nord >= 0) then
                v_nord = vn(i,j)
            else
                v_nord = vn(i,j+1)
            end if

            if (v_tilt_sud >= 0) then
                v_sud = vn(i,j-1)
            else
                v_sud = vn(i,j)
            end if

            ! Calcul des vitesses au temps n+1
            v_temp(i,j) = vn(i,j) - dt/dx*(u_tilt_est*v_est - u_tilt_west*v_west)&
                         - dt/dy*(v_tilt_nord*v_nord - v_tilt_sud*v_sud)
                         
            vn(i,j) = v_temp(i,j) + v_dif(i,j)    
                         
        end do
    end do    
end subroutine vitesse_upwind    

subroutine diffusion(un,vn,u_dif,v_dif,dx,dy,dt,nu,nx,ny)
    implicit none 
    integer, intent(in) :: nx,ny
    real*8, dimension(0:nx+1,0:ny+1), intent(in) :: un,vn
    real*8, intent(in) :: dx,dy,nu,dt
    real*8, dimension(0:nx+1,0:ny+1), intent(out) :: u_dif, v_dif
    integer :: i,j    
    u_dif=0
    v_dif=0
    do i=1,nx
        do j=1,ny
            u_dif(i,j) = ((nu/(dx*dx))*(un(i+1,j)-2*un(i,j)+un(i-1,j))+(nu/(dy*dy))*(un(i,j+1)-2*un(i,j)+un(i,j-1)))*dt
            v_dif(i,j) = ((nu/(dx*dx))*(vn(i+1,j)-2*vn(i,j)+vn(i-1,j))+(nu/(dy*dy))*(vn(i,j+1)-2*vn(i,j)+vn(i,j-1)))*dt
        end do
    end do

end subroutine diffusion


subroutine calcul_rhs(un,vn,rhs,dx,dy,dt,nx,ny)
    implicit none
    integer, intent(in) :: nx,ny
    real*8, dimension(0:nx+1,0:ny+1), intent(in) :: un,vn
    real*8, intent(in) :: dx,dy,dt
    real*8, dimension(0:nx+1,0:ny+1), intent(out) :: rhs
    real*8, dimension(1:nx,1:ny) :: div
    integer :: i,j
    call divergence(un,vn,dx,dy,nx,ny,div)
    do i=1,nx
        do j=1,ny
            rhs(i,j) = ((un(i,j) - un(i-1,j))/dx + (vn(i,j) - vn(i,j-1))/dy)/dt    
        end do
    end do
end subroutine calcul_rhs


subroutine divergence(un,vn,dx,dy,nx,ny,div)
    implicit none
    integer, intent(in) :: nx,ny
    real*8, dimension(0:nx+1,0:ny+1), intent(in) :: un,vn
    real*8, intent(in) :: dx,dy
    real*8, dimension(1:nx,1:ny), intent(out) :: div
    integer :: i,j

    do i=1,nx
        do j=1,ny
            div(i,j) = (un(i+1,j) - un(i,j))/dx + (vn(i,j+1) - vn(i,j))/dy
        end do
    end do

end subroutine divergence


subroutine rotationnel(un,vn,dx,dy,nx,ny,rot)
    implicit none
    integer, intent(in) :: nx,ny
    real*8, dimension(0:nx+1,0:ny+1), intent(in) :: un,vn
    real*8, intent(in) :: dx,dy
    real*8, dimension(1:nx,1:ny), intent(out) :: rot
    integer :: i,j

    do i=1,nx
        do j=1,ny
            rot(i,j) = (vn(i+1,j) - vn(i,j))/dx - (un(i,j+1) - un(i,j))/dy
        end do
    end do

end subroutine rotationnel
