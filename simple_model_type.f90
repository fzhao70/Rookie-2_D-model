module utils 
    use type_def

    implicit none

    contains

    subroutine message_display(str_in)
    !
    !>\ Author : Fanghe Zhao
    !>\ Date : 10/2018
    !>\ Details:
    !>  A warpper for print message in fix format 
    !
        implicit none
        character (len = 20), intent(in) :: str_in
        write(6, *) "==========================="
        write(6, *) str_in
        write(6, *) "==========================="

        return
    end subroutine message_display

    function distance(lon1, lat1, lon2, lat2)
    !
    !>\ Author : Fanghe Zhao
    !>\ Date : 10/2018
    !>\ Details:
    !> This function is for compute distance between two points 
    !
        implicit none

        real, parameter :: R_e = 6371 
        real, intent(in) :: lon1, lat1, lon2, lat2
        real :: distance
        real :: dlon, dlat
        real :: a, b
        dlon = lon2 - lon1
        dlat = lat2 - lat1
        a = SIN(dlat / 2)**2 + COS(lat1) * COS(lat2) * SIN(dlon / 2)**2
        b = 2 * ATAN2(SQRT(a), SQRT(1 - a))
        distance = R_e * b * 1000

        return
    end function distance

    function centered_diff(input, i, j, lon, lat, nx, ny, xyflags)
    !
    !>\ Author : Fanghe Zhao
    !>\ Date : 10/2018
    !>\ Details:
    !> This function is for compute centered_diff
    !
        implicit none

        integer, intent(in)  :: i, j, nx, ny, xyflags
        real, dimension(ny), intent(in) :: lat
        real, dimension(nx), intent(in) :: lon
        real, dimension(ny, nx), intent(in)  :: input
        real :: ds
        real :: centered_diff

        IF (xyflags .eq. 1) THEN
            ! Y direction
            ds = distance(lon(i), lat(j + 1), lon(i), lat(j - 1)) 
            centered_diff = (input(j + 1, i) - input(j - 1, i)) / (ds)  
        ELSEIF (xyflags .eq. 0) THEN
            ! X direction
            ds = distance(lon(i + 1), lat(j), lon(i - 1), lat(j)) 
            centered_diff = (input(j, i + 1) - input(j, i - 1)) / (ds)  
        ELSE
            call message_display("Wrong Arg input xyflags must be 0 for x , 1 for y")
        END IF

        return
    end function centered_diff

    function forward_diff(input, i, j, lon, lat, nx, ny, xyflags)
    !
    !>\ Author : Fanghe Zhao
    !>\ Date : 10/2018
    !>\ Details:
    !> This function is for compute forward_diff
    !
        implicit none

        integer, intent(in)  :: i, j, nx, ny, xyflags
        real, dimension(ny), intent(in) :: lat
        real, dimension(nx), intent(in) :: lon
        real, dimension(ny, nx), intent(in)  :: input
        real :: ds
        real :: forward_diff

        IF (xyflags .eq. 1) THEN
            ! Y direction
            ds = distance(lon(i), lat(j), lon(i), lat(j + 1)) 
            forward_diff = (input(j + 1, i) - input(j, i)) / (ds)  
        ELSEIF (xyflags .eq. 0) THEN
            ! X direction
            ds = distance(lon(i), lat(j), lon(i + 1), lat(j)) 
            forward_diff = (input(j, i + 1) - input(j, i)) / (ds)  
        ELSE
            call message_display("Wrong Arg input xyflags must be 0 for x , 1 for y")
        END IF

        return
    end function forward_diff

    function backward_diff(input, i, j, lon, lat, nx, ny, xyflags)
    !
    !>\ Author : Fanghe Zhao
    !>\ Date : 10/2018
    !>\ Details:
    !> This function is for compute backward_diff
    !
        implicit none

        integer, intent(in)  :: i, j, nx, ny, xyflags
        real, dimension(ny), intent(in) :: lat
        real, dimension(nx), intent(in) :: lon
        real, dimension(ny, nx), intent(in)  :: input
        real :: ds
        real :: backward_diff

        IF (xyflags .eq. 1) THEN
            ! Y direction
            ds = distance(lon(i), lat(j - 1), lon(i), lat(j)) 
            backward_diff = (input(j, i) - input(j - 1, i)) / (ds)  
        ELSEIF (xyflags .eq. 0) THEN
            ! X direction
            ds = distance(lon(i), lat(j), lon(i - 1), lat(j)) 
            backward_diff = (input(j, i) - input(j, i - 1)) / (ds)  
        ELSE
            call message_display("Wrong Arg input xyflags must be 0 for x , 1 for y")
        END IF

        return
    end function backward_diff

    function horizontal_diff(source, lon, lat, nx, ny, xyflags) result(diff)
    !
    !>\ Author : Fanghe Zhao
    !>\ Date : 10/2018
    !>\ Details:
    !> This function is for compute backward_diff
    !
        implicit none
        integer, intent(in)  :: nx, ny, xyflags
        integer :: i, j
        real, dimension(ny), intent(in) :: lat
        real, dimension(nx), intent(in) :: lon
        real, dimension(ny, nx), intent(in)  :: source
        real, dimension(ny, nx)  :: diff

        !$OMP PARALLEL DO
        do i = 1, nx
            do j = 1, ny
                if (i .eq. 0 .and. j .eq. 0) then
                    diff(j, i) = forward_diff(source, i, j, lon, lat, nx, ny, xyflags)

                else if (i .eq. nx .and. j .eq. ny) then
                    diff(j, i) = backward_diff(source, i, j, lon, lat, nx, ny, xyflags)

                else if (i .eq. nx .and. j .eq. 0) then
                    if (xyflags .eq. 0) then
                        diff(j, i) = backward_diff(source, i, j, lon, lat, nx, ny, xyflags)
                    elseif (xyflags .eq. 1) then
                        diff(j, i) = forward_diff(source, i, j, lon, lat, nx, ny, xyflags)
                    end if

                else if (i .eq. 0 .and. j .eq. ny) then
                    if (xyflags .eq. 0) then
                        diff(j, i) = forward_diff(source, i, j, lon, lat, nx, ny, xyflags)
                    elseif (xyflags .eq. 1) then
                        diff(j, i) = backward_diff(source, i, j, lon, lat, nx, ny, xyflags)
                    end if

                else
                        diff(j, i) = centered_diff(source, i, j, lon, lat, nx, ny, xyflags)
                end if
            end do
        end do
        !$OMP END PARALLEL DO

        return
    end function horizontal_diff

end module

module solver
    use type_def
    use utils
    implicit none
    integer, parameter :: X_AXIS = 0, Y_AXIS = 1
    integer, parameter :: g = 9.8, H = 8000
    real, parameter :: omega = 7.292e-5 

    contains

    subroutine Matsuno(u, v, Z, u_next, v_next, Z_next, lon, lat, nx, ny)
        implicit none
        integer, intent(in)  :: nx, ny
        integer :: i, j
        real :: dt, f
        real, dimension(ny), intent(in) :: lat
        real, dimension(nx), intent(in) :: lon
        real, dimension(ny, nx), intent(in)  :: Z, u, v 
        real, dimension(ny, nx) :: u_next, v_next, Z_next 
        real, dimension(ny, nx) :: u_star, v_star, Z_star
        real, dimension(ny, nx) :: u_diff_x, u_diff_y, Z_diff_x, Z_diff_y, v_diff_x, v_diff_y
        real, dimension(ny, nx) :: u_star_diff_x, u_star_diff_y, Z_star_diff_x, Z_star_diff_y, v_star_diff_x, v_star_diff_y

        u_diff_x = horizontal_diff(u, lon, lat, nx, ny, X_AXIS) 
        u_diff_y = horizontal_diff(u, lon, lat, nx, ny, Y_AXIS) 
        v_diff_x = horizontal_diff(v, lon, lat, nx, ny, X_AXIS) 
        v_diff_y = horizontal_diff(v, lon, lat, nx, ny, Y_AXIS) 
        Z_diff_x = horizontal_diff(Z, lon, lat, nx, ny, X_AXIS) 
        Z_diff_y = horizontal_diff(Z, lon, lat, nx, ny, Y_AXIS) 

        u_star = u + dt * (-u *  u_diff_x - v * u_diff_y - g * Z_diff_x + f * v) 
        v_star = v + dt * (-u *  v_diff_x - v * v_diff_y - g * Z_diff_x - f * u) 
        Z_star = Z + dt * (-u *  Z_diff_x - v * Z_diff_y - H * (u_diff_x + v_diff_y))

        u_star_diff_x = horizontal_diff(u_star, lon, lat, nx, ny, X_AXIS)
        u_star_diff_y = horizontal_diff(u_star, lon, lat, nx, ny, Y_AXIS)
        v_star_diff_x = horizontal_diff(v_star, lon, lat, nx, ny, X_AXIS)
        v_star_diff_y = horizontal_diff(v_star, lon, lat, nx, ny, Y_AXIS)
        Z_star_diff_x = horizontal_diff(Z_star, lon, lat, nx, ny, X_AXIS)
        Z_star_diff_y = horizontal_diff(Z_star, lon, lat, nx, ny, Y_AXIS)

        u_next = u + dt * (-u_star * u_star_diff_x - v * u_star_diff_y - g * Z_diff_x + f * v)
        v_next = v + dt * (-u_star * v_star_diff_x - v * v_star_diff_y - g * Z_diff_x + f * v)
        Z_next = Z + dt * (-u *  Z_star_diff_x - v * Z_star_diff_y - H * (u_diff_x + v_diff_y))

        return
    end subroutine Matsuno

    subroutine Leapfrog(u, v, Z, u_pre, v_pre, Z_pre, u_next, v_next, Z_next, lon, lat, nx, ny)
        implicit none
        integer, intent(in)  :: nx, ny
        integer :: i, j
        real :: dt, f
        real, dimension(ny), intent(in) :: lat
        real, dimension(nx), intent(in) :: lon
        real, dimension(ny, nx), intent(in) :: Z, u, v 
        real, dimension(ny, nx), intent(in) :: u_pre, v_pre, Z_pre
        real, dimension(ny, nx) :: u_next, v_next, Z_next 
        real, dimension(ny, nx) :: u_diff_x, u_diff_y, Z_diff_x, Z_diff_y, v_diff_x, v_diff_y

        u_diff_x = horizontal_diff(u, lon, lat, nx, ny, X_AXIS) 
        u_diff_y = horizontal_diff(u, lon, lat, nx, ny, Y_AXIS) 
        v_diff_x = horizontal_diff(v, lon, lat, nx, ny, X_AXIS) 
        v_diff_y = horizontal_diff(v, lon, lat, nx, ny, Y_AXIS) 
        Z_diff_x = horizontal_diff(Z, lon, lat, nx, ny, X_AXIS) 
        Z_diff_y = horizontal_diff(Z, lon, lat, nx, ny, Y_AXIS) 

        u_next = u_pre + 2 * dt * (-u * u_diff_x - v * u_diff_y - g * Z_diff_x + f * v)
        v_next = v_pre + 2 * dt * (-u * v_diff_x - v * v_diff_y - g * Z_diff_x - f * u)
        Z_next = Z_pre + 2 * dt * (-u * Z_diff_x - v * Z_diff_y - H * (u_diff_x + v_diff_y))

        return
    end subroutine 

end module solver

module io
    use netcdf
    use utils 
    implicit none

    character (len = *), parameter :: IN_FILE_NAME = "hgt_location_selected.nc"
    character (len = *), parameter :: OUT_FILE_NAME = "geo_wind.nc"
    integer :: ncid, hgtid, lonid, latid, uid, vid
    integer :: y_dimid, x_dimid
    integer(kind=1) :: ierr

    contains

end module io

module type_def
    implicit none

    integer, parameter :: nx = 33,ny = 23
    integer :: nt = 1000

    type element
    !meteorological elements
        real, dimension(ny, nx) :: Z = 0, u = 0, v = 0
    end type element

    type grid
    !meteorological elements
        real, dimension(ny) :: lat
        real, dimension(nx) :: lon
    end type grid

end module

program simple_model
    !
    !>\ Author : Fanghe Zhao
    !>\ Date : 10/2018
    !>\ Details:
    !>  This Program is for compute geostrophic wind use ECMWF data
    !
    use utils 
    use solver
    use io
    use type_def

    implicit none

    ! Var Defination
    integer :: i, j
    real :: dx = 1, dy = 1, dt = 1
    real :: f

    type(grid), target :: coor
    type(element), target :: ele_pre, ele, ele_next

    real, dimension(:, :), pointer :: Z, u, v
    real, dimension(:, :), pointer :: Z_next, u_next, v_next
    real, dimension(:, :), pointer :: Z_pre, u_pre, v_pre
    real, dimension(:), pointer :: lat
    real, dimension(:), pointer :: lon

    ! pointer assignment
    Z => ele%Z
    u => ele%u
    v => ele%v
    Z_next => ele_next%Z
    u_next => ele_next%u
    v_next => ele_next%v
    Z_pre => ele_pre%Z
    u_pre => ele_pre%u
    v_pre => ele_pre%v
    lon => coor%lon
    lat => coor%lat

    do i = 1, nt
        if (i .eq. 1) then
            call Matsuno(u, v, Z, u_next, v_next, Z_next, lon, lat, nx, ny)
            !update_var
            ele = ele_next
        else 
            call Leapfrog(u, v, Z, u_pre, v_pre, Z_pre, u_next, v_next, Z_next, lon, lat, nx, ny)
            !update_var
            ele_pre = ele
            ele = ele_next
        end if
    end do

end program simple_model