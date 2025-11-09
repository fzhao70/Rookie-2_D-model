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

    subroutine Matsuno(u, v, Z, u_next, v_next, Z_next, lon, lat, nx, ny, dt_in)
        implicit none
        integer, intent(in)  :: nx, ny
        integer :: i, j
        real, intent(in) :: dt_in
        real :: dt, f, lat_mean
        real, dimension(ny), intent(in) :: lat
        real, dimension(nx), intent(in) :: lon
        real, dimension(ny, nx), intent(in)  :: Z, u, v
        real, dimension(ny, nx) :: u_next, v_next, Z_next
        real, dimension(ny, nx) :: u_star, v_star, Z_star
        real, dimension(ny, nx) :: u_diff_x, u_diff_y, Z_diff_x, Z_diff_y, v_diff_x, v_diff_y
        real, dimension(ny, nx) :: u_star_diff_x, u_star_diff_y, Z_star_diff_x, Z_star_diff_y, v_star_diff_x, v_star_diff_y

        dt = dt_in
        lat_mean = SUM(lat) / ny
        f = 2 * omega * SIN(lat_mean)

        u_diff_x = horizontal_diff(u, lon, lat, nx, ny, X_AXIS)
        u_diff_y = horizontal_diff(u, lon, lat, nx, ny, Y_AXIS)
        v_diff_x = horizontal_diff(v, lon, lat, nx, ny, X_AXIS)
        v_diff_y = horizontal_diff(v, lon, lat, nx, ny, Y_AXIS)
        Z_diff_x = horizontal_diff(Z, lon, lat, nx, ny, X_AXIS)
        Z_diff_y = horizontal_diff(Z, lon, lat, nx, ny, Y_AXIS)

        u_star = u + dt * (-u *  u_diff_x - v * u_diff_y - g * Z_diff_x + f * v)
        v_star = v + dt * (-u *  v_diff_x - v * v_diff_y - g * Z_diff_y - f * u)
        Z_star = Z + dt * (-u *  Z_diff_x - v * Z_diff_y - H * (u_diff_x + v_diff_y))

        u_star_diff_x = horizontal_diff(u_star, lon, lat, nx, ny, X_AXIS)
        u_star_diff_y = horizontal_diff(u_star, lon, lat, nx, ny, Y_AXIS)
        v_star_diff_x = horizontal_diff(v_star, lon, lat, nx, ny, X_AXIS)
        v_star_diff_y = horizontal_diff(v_star, lon, lat, nx, ny, Y_AXIS)
        Z_star_diff_x = horizontal_diff(Z_star, lon, lat, nx, ny, X_AXIS)
        Z_star_diff_y = horizontal_diff(Z_star, lon, lat, nx, ny, Y_AXIS)

        u_next = u + dt * (-u_star * u_star_diff_x - v_star * u_star_diff_y - g * Z_star_diff_x + f * v_star)
        v_next = v + dt * (-u_star * v_star_diff_x - v_star * v_star_diff_y - g * Z_star_diff_y - f * u_star)
        Z_next = Z + dt * (-u_star * Z_star_diff_x - v_star * Z_star_diff_y - H * (u_star_diff_x + v_star_diff_y))

        return
    end subroutine Matsuno

    subroutine Leapfrog(u, v, Z, u_pre, v_pre, Z_pre, u_next, v_next, Z_next, lon, lat, nx, ny, dt_in)
        implicit none
        integer, intent(in)  :: nx, ny
        integer :: i, j
        real, intent(in) :: dt_in
        real :: dt, f, lat_mean
        real, dimension(ny), intent(in) :: lat
        real, dimension(nx), intent(in) :: lon
        real, dimension(ny, nx), intent(in) :: Z, u, v
        real, dimension(ny, nx), intent(in) :: u_pre, v_pre, Z_pre
        real, dimension(ny, nx) :: u_next, v_next, Z_next
        real, dimension(ny, nx) :: u_diff_x, u_diff_y, Z_diff_x, Z_diff_y, v_diff_x, v_diff_y

        dt = dt_in
        lat_mean = SUM(lat) / ny
        f = 2 * omega * SIN(lat_mean)

        u_diff_x = horizontal_diff(u, lon, lat, nx, ny, X_AXIS)
        u_diff_y = horizontal_diff(u, lon, lat, nx, ny, Y_AXIS)
        v_diff_x = horizontal_diff(v, lon, lat, nx, ny, X_AXIS)
        v_diff_y = horizontal_diff(v, lon, lat, nx, ny, Y_AXIS)
        Z_diff_x = horizontal_diff(Z, lon, lat, nx, ny, X_AXIS)
        Z_diff_y = horizontal_diff(Z, lon, lat, nx, ny, Y_AXIS)

        u_next = u_pre + 2 * dt * (-u * u_diff_x - v * u_diff_y - g * Z_diff_x + f * v)
        v_next = v_pre + 2 * dt * (-u * v_diff_x - v * v_diff_y - g * Z_diff_y - f * u)
        Z_next = Z_pre + 2 * dt * (-u * Z_diff_x - v * Z_diff_y - H * (u_diff_x + v_diff_y))

        return
    end subroutine Leapfrog 

end module solver

module io
    use netcdf
    use utils
    use type_def
    implicit none

    character (len = *), parameter :: IN_FILE_NAME = "hgt_location_selected.nc"
    character (len = *), parameter :: OUT_FILE_NAME = "simple_model_output.nc"
    integer :: ncid, hgtid, lonid, latid, uid, vid, zid
    integer :: y_dimid, x_dimid, t_dimid
    integer(kind=1) :: ierr

    contains

    subroutine read_input(lon, lat, hgt)
    !
    !>\ Author : Fanghe Zhao
    !>\ Date : 10/2018
    !>\ Details:
    !>  Read input NetCDF file with geopotential height and coordinates
    !
        implicit none
        real, dimension(ny), intent(out) :: lat
        real, dimension(nx), intent(out) :: lon
        real, dimension(ny, nx), intent(out) :: hgt

        include 'netcdf.inc'

        call message_display("Reading Input File")
        ncid = 100
        ierr = nf90_open(path = IN_FILE_NAME, mode = NF90_NOWRITE, ncid = ncid)
        if (ierr /= NF90_NOERR) then
            call message_display("Error opening input file")
            stop
        end if

        ierr = nf_inq_varid(ncid, 'hgt', hgtid)
        ierr = nf_inq_varid(ncid, 'lon', lonid)
        ierr = nf_inq_varid(ncid, 'lat', latid)

        ierr = nf_get_var_real(ncid, lonid, lon)
        ierr = nf_get_var_real(ncid, latid, lat)
        ierr = nf_get_var_real(ncid, hgtid, hgt)

        ierr = nf_close(ncid)
        call message_display("Input File Read")

        return
    end subroutine read_input

    subroutine write_output(lon, lat, u, v, Z, time_step)
    !
    !>\ Author : Fanghe Zhao
    !>\ Date : 10/2018
    !>\ Details:
    !>  Write output NetCDF file with model results
    !
        implicit none
        integer, intent(in) :: time_step
        integer, dimension(2) :: dimids
        real, dimension(ny), intent(in) :: lat
        real, dimension(nx), intent(in) :: lon
        real, dimension(ny, nx), intent(in) :: u, v, Z

        include 'netcdf.inc'

        call message_display("Writing Output File")

        ierr = nf90_create(OUT_FILE_NAME, NF90_CLOBBER, ncid)
        ierr = nf90_def_dim(ncid, "lon", nx, x_dimid)
        ierr = nf90_def_dim(ncid, "lat", ny, y_dimid)
        dimids = (/ y_dimid, x_dimid /)

        ierr = nf90_def_var(ncid, "u", NF90_FLOAT, dimids, uid)
        ierr = nf90_def_var(ncid, "v", NF90_FLOAT, dimids, vid)
        ierr = nf90_def_var(ncid, "Z", NF90_FLOAT, dimids, zid)
        ierr = nf90_def_var(ncid, "lon", NF90_FLOAT, x_dimid, lonid)
        ierr = nf90_def_var(ncid, "lat", NF90_FLOAT, y_dimid, latid)

        ierr = nf90_enddef(ncid)

        ierr = nf90_put_var(ncid, uid, u)
        ierr = nf90_put_var(ncid, vid, v)
        ierr = nf90_put_var(ncid, zid, Z)
        ierr = nf90_put_var(ncid, lonid, lon)
        ierr = nf90_put_var(ncid, latid, lat)

        ierr = nf90_close(ncid)
        call message_display("Output File Written")

        return
    end subroutine write_output

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
    !>  This Program is a simple 2D atmospheric model using shallow water equations
    !
    use utils
    use solver
    use io
    use type_def

    implicit none

    ! Var Defination
    integer :: i, j
    real :: dx = 1, dy = 1, dt = 100.0
    real :: f
    real, dimension(ny, nx) :: hgt

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

    ! Read input data
    call read_input(lon, lat, hgt)

    ! Initialize height field from geopotential height
    Z = hgt

    ! Initialize u and v from geostrophic balance
    call message_display("Initializing winds")
    do i = 1, nx
        do j = 1, ny
            f = 2 * omega * SIN(lat(j))
            if (i .eq. 1 .and. j .eq. 1) then
                dx = distance(lon(i), lat(j), lon(i + 1), lat(j))
                dy = distance(lon(i), lat(j), lon(i), lat(j + 1))
                u(j, i) = -1 * (g / f) * (hgt(j + 1, i) - hgt(j, i)) / dy
                v(j, i) = (g / f) * (hgt(j, i + 1) - hgt(j, i)) / dx
            else if (i .eq. nx .and. j .eq. ny) then
                dx = distance(lon(i - 1), lat(j), lon(i), lat(j))
                dy = distance(lon(i), lat(j - 1), lon(i), lat(j))
                u(j, i) = -1 * (g / f) * (hgt(j, i) - hgt(j - 1, i)) / dy
                v(j, i) = (g / f) * (hgt(j, i) - hgt(j, i - 1)) / dx
            else if (i .eq. nx .and. j .eq. 1) then
                dx = distance(lon(i - 1), lat(j), lon(i), lat(j))
                dy = distance(lon(i), lat(j), lon(i), lat(j + 1))
                u(j, i) = -1 * (g / f) * (hgt(j + 1, i) - hgt(j, i)) / dy
                v(j, i) = (g / f) * (hgt(j, i) - hgt(j, i - 1)) / dx
            else if (i .eq. 1 .and. j .eq. ny) then
                dx = distance(lon(i), lat(j), lon(i + 1), lat(j))
                dy = distance(lon(i), lat(j - 1), lon(i), lat(j))
                u(j, i) = -1 * (g / f) * (hgt(j, i) - hgt(j - 1, i)) / dy
                v(j, i) = (g / f) * (hgt(j, i + 1) - hgt(j, i)) / dx
            else
                dx = distance(lon(i - 1), lat(j), lon(i + 1), lat(j))
                dy = distance(lon(i), lat(j - 1), lon(i), lat(j + 1))
                u(j, i) = -1 * (g / f) * (hgt(j + 1, i) - hgt(j - 1, i)) / dy
                v(j, i) = (g / f) * (hgt(j, i + 1) - hgt(j, i - 1)) / dx
            end if
        end do
    end do

    ! Time integration
    call message_display("Starting time integration")
    do i = 1, nt
        if (i .eq. 1) then
            call Matsuno(u, v, Z, u_next, v_next, Z_next, lon, lat, nx, ny, dt)
            ! Update variables
            ele_pre = ele
            ele = ele_next
        else
            call Leapfrog(u, v, Z, u_pre, v_pre, Z_pre, u_next, v_next, Z_next, lon, lat, nx, ny, dt)
            ! Update variables
            ele_pre = ele
            ele = ele_next
        end if

        ! Print progress every 100 time steps
        if (MOD(i, 100) .eq. 0) then
            write(6, *) "Time step:", i, "/", nt
        end if
    end do
    call message_display("Time integration complete")

    ! Write output
    call write_output(lon, lat, u, v, Z, nt)

end program simple_model