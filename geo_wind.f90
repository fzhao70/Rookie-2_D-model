!------------------------------------------------------------------------------
! USTC-AEMOL, geostrophic_wind_4_teaching
!------------------------------------------------------------------------------
!
!> @author
!> Fanghe Zhao}
!> zfh1997 at mail.ustc.edu.cn}
!
! DESCRIPTION: 
!>  A simple model just for solving geostrophic_wind
!
! LICENSE: MIT
!
! REVISION HISTORY:
! 02 10 2018 - Initial Version
!------------------------------------------------------------------------------
module utils 
    implicit none

    contains

    function distance(lon1, lat1, lon2, lat2)
    !
    ! This function is for compute distance between two points 
    !
    ! Author : Fanghe Zhao
        implicit none

        real, parameter :: R_e = 6371 
        real :: lon1, lat1, lon2, lat2
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

    subroutine message_display(str_in)
    !
    ! This function is for print message in fix format 
    !
    ! Author : Fanghe Zhao
        implicit none
        character (len = 20), intent(in) :: str_in
        write(6, *) "==========================="
        write(6, *) str_in
        write(6, *) "==========================="
    end subroutine message_display

end module

program geostrophic_wind
    !
    ! This Program is for compute geostrophic wind use ECMWF data
    !
    ! Author : Fanghe Zhao
    use netcdf
    use utils 
    implicit none
    ! Var Defination
    integer, parameter :: nx = 33,ny = 23
    real, parameter :: g = 9.8 
    real, parameter :: omega = 7.292e-5 
    character (len = *), parameter :: IN_FILE_NAME = "hgt_location_selected.nc"
    character (len = *), parameter :: OUT_FILE_NAME = "geo_wind.nc"
    integer(kind=1) :: ierr
    integer :: i, j
    integer :: ncid, hgtid, lonid, latid, uid, vid
    integer :: y_dimid, x_dimid
    real :: dx, dy ,f
    integer,dimension(2) :: dimids
    real, dimension(ny, nx) :: hgt = 0, u = 0, v = 0
    real, dimension(ny) :: lat
    real, dimension(nx) :: lon

    include 'netcdf.inc'
    
    ! Read input netcdf file
    call message_display("Read Input File")
    ncid = 100 
    ierr = nf90_open(path = IN_FILE_NAME, mode = NF90_NOWRITE ,ncid = ncid)
    ierr = nf_inq_varid(ncid, 'hgt', hgtid)
    ierr = nf_inq_varid(ncid, 'lon', lonid)
    ierr = nf_inq_varid(ncid, 'lat', latid)

    ierr = nf_get_var_real(ncid, lonid, lon)
    ierr = nf_get_var_real(ncid, latid, lat)
    ierr = nf_get_var_real(ncid, hgtid, hgt)

    ierr = nf_close(ncid)
    call message_display("Read Input File End")

    ! Compute Geo-wind
    call message_display("Compute Geo-wind")
    do i = 1, nx
        do j = 1, ny
            f = (2 * omega * SIN(lat(j)))
            if (i .eq. 0 .and. j .eq. 0) then
                dx = distance(lon(i), lat(j), lon(i + 1), lat(j)) 
                dy = distance(lon(i), lat(j), lon(i), lat(j + 1)) 
                u(j, i) = -1 * (g / f) * (hgt(j + 1, i) - hgt(j, i)) / (dy) 
                v(j, i) = (g / f) * (hgt(j, i + 1) - hgt(j, i)) / (dx) 
            else if (i .eq. nx .and. j .eq. ny) then
                dx = distance(lon(i - 1), lat(j), lon(i), lat(j)) 
                dy = distance(lon(i), lat(j - 1), lon(i), lat(j)) 
                u(j, i) = -1 * (g / f) * (hgt(j, i) - hgt(j - 1, i)) / (dy)  
                v(j, i) = (g / f) * (hgt(j, i) - hgt(j, i - 1)) / (dx) 
            else if (i .eq. nx .and. j .eq. 0) then
                dx = distance(lon(i - 1), lat(j), lon(i), lat(j)) 
                dy = distance(lon(i), lat(j), lon(i), lat(j + 1)) 
                u(j, i) = -1 * (g / f) * (hgt(j + 1, i) - hgt(j, i)) / (dy)  
                v(j, i) = (g / f) * (hgt(j, i) - hgt(j, i - 1)) / (dx) 
            else if (i .eq. 0 .and. j .eq. ny) then
                dx = distance(lon(i), lat(j), lon(i + 1), lat(j)) 
                dy = distance(lon(i), lat(j - 1), lon(i), lat(j)) 
                u(j, i) = -1 * (g / f) * (hgt(j, i) - hgt(j - 1, i)) / (dy)  
                v(j, i) = (g / f) * (hgt(j, i + 1) - hgt(j, i)) / (dx) 
            else
                dx = distance(lon(i - 1), lat(j), lon(i + 1), lat(j)) 
                dy = distance(lon(i), lat(j - 1), lon(i), lat(j + 1)) 
                u(j, i) = -1 * (g / f) * (hgt(j + 1, i) - hgt(j - 1, i)) / (dy)  
                v(j, i) = (g / f) * (hgt(j, i + 1) - hgt(j, i - 1)) / (dx) 
            end if
        end do
    end do
    call message_display("Compute Geo-wind End")

    ! Output File
    call message_display("Write Output File")
    ierr = nf90_create(OUT_FILE_NAME, NF90_CLOBBER, ncid)
    ierr = nf90_def_dim(ncid, "lon", nx, x_dimid)
    ierr = nf90_def_dim(ncid, "lat", ny, y_dimid)
    dimids =  (/ y_dimid, x_dimid /)
    ierr = nf90_def_var(ncid, "u", NF90_FLOAT, dimids, uid)
    ierr = nf90_def_var(ncid, "v", NF90_FLOAT, dimids, vid)
    ierr = nf90_def_var(ncid, "lon", NF90_INT, x_dimid, lonid)
    ierr = nf90_def_var(ncid, "lat", NF90_INT, y_dimid, latid)
    ierr = nf90_enddef(ncid)
    ierr = nf90_put_var(ncid, uid, u)
    ierr = nf90_put_var(ncid, vid, v)
    ierr = nf90_put_var(ncid, lonid, lon)
    ierr = nf90_put_var(ncid, latid, lat)
    ierr = nf90_close(ncid)
    call message_display("Write Output File End")

end program geostrophic_wind
