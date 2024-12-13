subroutine initialize(un,vn,nx,ny) 
    implicit none
    real*8, dimension(0:nx+1,0:ny+1), intent(out) :: un,vn 
    integer, intent(in) :: nx,ny
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
    real*8, dimension(0:nx+1,0:ny+1), intent(out) :: un,vn 
    integer, intent(in) :: nx,ny
    integer :: i,j
    do i=0,nx+1
        do j=0,ny+1
            un(i,j)=0.0
            vn(i,j)=0.0
        end do
    end do
    ! /!\ ATTENTION CONDITIONS LIMITES A IMPLEMENTER
    do i=0,nx+1
        do j=0,ny+1  
            ! Condition Bas
            if (j == 0) then
                un(i,0) = un(i,1)
                vn(i,0) = 0
            end if
            
            ! Condition Gauche
            if (i == 0) then
                un(0,j) = 0
                vn(0,j) = - vn(1,j)
            end if

            ! Condition Droite
            if (i == nx + 1) then
                un(nx,j) = 0
                vn(nx+1,j) = - vn(nx,j)
            end if

            ! Condition Haut
            if (j == ny + 1) then
                un(i,ny+1) = 2 - un(i,ny)
                vn(i,ny) = 0
            end if
        end do
    end do
    do i=0,ny+1
    write(*,'(6F10.2)') un(0:nx+1, i)
    end do
end subroutine initialize_un


subroutine timestep(un,vn,nx,ny,dx,dy,dt,nu)
    ! Schéma Upwind et Centré
    implicit none
    integer, intent(in) :: nx,ny
    real*8, dimension(0:nx+1,0:ny+1), intent(in) :: un,vn
    real*8, intent(in) :: dx,dy,nu
    real*8, intent(out) :: dt
    real*8 :: umax, vmax

    umax = maxval(abs(un(:,:)))
    vmax = maxval(abs(vn(:,:)))
    dt = 1.0/(2*(umax/dx + vmax/dy) + 2*nu*(1/(dx*dx) + 1/(dy*dy)))
end subroutine timestep
            

subroutine vitesse_tilt(un,vn,u_tilt,v_tilt,nx,ny)
    implicit none
    real*8, dimension(0:nx+1,0:ny+1), intent(in) :: un,vn
    integer, intent(in) :: nx,ny
    real*8, dimension(0:nx+1,0:ny+1), intent(out) :: u_tilt,v_tilt
    integer :: i,j
    do i=1,nx
        do j=1,ny
            u_tilt(i,j) = 0.5*(un(i,j)+un(i+1,j))
            v_tilt(i,j) = 0.5*(vn(i,j)+vn(i,j+1))
        end do
    end do
end subroutine vitesse_tilt


subroutine vitesse_upwind(un,vn,u_tilt,v_tilt,u_dif,v_dif,nx,ny, dt, dx, dy)
    implicit none
    integer, intent(in) :: nx,ny
    real*8, dimension(0:nx+1,0:ny+1), intent(inout) :: un,vn
    real*8, dimension(0:nx+1,0:ny+1), intent(in) :: u_tilt,v_tilt,u_dif,v_dif
    real*8, intent(in) :: dt,dx,dy
    real*8, dimension(0:nx+1,0:ny+1) :: u_temp, v_temp
    real*8 :: u_ouest, u_est, u_nord, u_sud
    real*8 :: v_ouest, v_est, v_nord, v_sud
    integer :: i,j

    call initialize_un(u_temp,v_temp,nx,ny)
    
    ! Calcul U
    do i=1,nx-1
        do j=1,ny
            ! Calcul des U au centre
            if (u_tilt(i,j) >= 0) then
                u_est = un(i,j)
            else
                u_est = un(i+1,j)
            end if

            if (u_tilt(i-1,j) >= 0) then
                u_ouest = un(i-1,j)
            else
                u_ouest = un(i,j)
            end if

            if (v_tilt(i,j) >= 0) then
                u_nord = un(i,j)
            else
                u_nord = un(i,j+1)
            end if

            if (v_tilt(i,j-1) >= 0) then
                u_sud = un(i,j-1)
            else
                u_sud = un(i,j)
            end if

            ! Calcul des V au centre
            if (u_tilt(i,j) >= 0) then
                v_est = vn(i,j)
            else
                v_est = vn(i+1,j)
            end if

            if (u_tilt(i-1,j) >= 0) then
                v_ouest = vn(i-1,j)
            else
                v_ouest = vn(i,j)
            end if

            if (v_tilt(i,j) >= 0) then
                v_nord = vn(i,j)
            else
                v_nord = vn(i,j+1)
            end if

            if (v_tilt(i,j-1) >= 0) then
                v_sud = vn(i,j-1)
            else
                v_sud = vn(i,j)
            end if

            ! Calcul des vitesses au temps n+1
            u_temp(i,j) = un(i,j) - dt/dx*(u_tilt(i,j)*u_est - u_tilt(i-1,j)*u_ouest)&
                         - dt/dy*(v_tilt(i,j)*u_nord - v_tilt(i,j-1)*u_sud)
        end do
    end do  

    ! Calcul V
    do i=1,nx
        do j=1,ny-1
            ! Calcul des U au centre
            if (u_tilt(i,j) >= 0) then
                u_est = un(i,j)
            else
                u_est = un(i+1,j)
            end if

            if (u_tilt(i-1,j) >= 0) then
                u_ouest = un(i-1,j)
            else
                u_ouest = un(i,j)
            end if

            if (v_tilt(i,j) >= 0) then
                u_nord = un(i,j)
            else
                u_nord = un(i,j+1)
            end if

            if (v_tilt(i,j-1) >= 0) then
                u_sud = un(i,j-1)
            else
                u_sud = un(i,j)
            end if

            ! Calcul des V au centre
            if (u_tilt(i,j) >= 0) then
                v_est = vn(i,j)
            else
                v_est = vn(i+1,j)
            end if

            if (u_tilt(i-1,j) >= 0) then
                v_ouest = vn(i-1,j)
            else
                v_ouest = vn(i,j)
            end if

            if (v_tilt(i,j) >= 0) then
                v_nord = vn(i,j)
            else
                v_nord = vn(i,j+1)
            end if

            if (v_tilt(i,j-1) >= 0) then
                v_sud = vn(i,j-1)
            else
                v_sud = vn(i,j)
            end if

            ! Calcul des vitesses au temps n+1
            v_temp(i,j) = vn(i,j) - dt/dx*(u_tilt(i,j)*v_est - u_tilt(i-1,j)*v_ouest)&
                         - dt/dy*(v_tilt(i,j)*v_nord - v_tilt(i,j-1)*v_sud)
        end do
    end do

    do i=0,nx+1
        do j=0,ny+1
            un(i,j) = u_temp(i,j) + u_dif(i,j)
            vn(i,j) = v_temp(i,j) + v_dif(i,j)
        end do 
    end do
end subroutine vitesse_upwind    

subroutine diffusion(un,vn,u_dif,v_dif,dx,dy,nu,nx,ny)
    implicit none 
    real*8, dimension(0:nx+1,0:ny+1), intent(in) :: un,vn
    real*8, intent(in) :: dx,dy,nu
    integer, intent(in) :: nx,ny
    real*8, dimension(0:nx+1,0:ny+1), intent(out) :: u_dif, v_dif
    integer :: i,j    
    call initialize(u_dif,v_dif,nx,ny)
    do i=1,nx
        do j=1,ny
            u_dif(i,j) = (nu/(dx*dx))*(un(i+1,j)-2*un(i,j)+un(i-1,j))+(nu/(dy*dy))*(un(i,j+1)-2*un(i,j)+un(i,j-1))
            v_dif(i,j) = (nu/(dx*dx))*(vn(i+1,j)-2*vn(i,j)+vn(i-1,j))+(nu/(dy*dy))*(vn(i,j+1)-2*vn(i,j)+vn(i,j-1))
        end do
    end do

end subroutine diffusion


subroutine calcul_rhs(un,vn,dx,dy,dt,nx,ny,rhs)
    implicit none
    real*8, dimension(0:nx+1,0:ny+1), intent(in) :: un,vn
    real*8, intent(in) :: dx,dy,dt
    integer, intent(in) :: nx,ny
    real*8, dimension(0:nx+1,0:ny+1), intent(out) :: rhs
    integer :: i,j
    do i=1,nx
        do j=1,ny
            rhs(i,j) = ((un(i,j) - un(i-1,j))/dx + (vn(i,j) - vn(i,j-1))/dy)/dt    
        end do
    end do
end subroutine calcul_rhs


subroutine divergence(un,vn,dx,dy,nx,ny,div)
    implicit none
    real*8, dimension(0:nx+1,0:ny+1), intent(in) :: un,vn
    real*8, intent(in) :: dx,dy
    integer, intent(in) :: nx,ny
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
    real*8, dimension(0:nx+1,0:ny+1), intent(in) :: un,vn
    real*8, intent(in) :: dx,dy
    integer, intent(in) :: nx,ny
    real*8, dimension(1:nx,1:ny), intent(out) :: rot
    integer :: i,j

    do i=1,nx
        do j=1,ny
            rot(i,j) = (vn(i+1,j) - vn(i,j))/dx - (un(i,j+1) - un(i,j))/dy
        end do
    end do

end subroutine rotationnel