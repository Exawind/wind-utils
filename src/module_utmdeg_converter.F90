!------------------------------------------------------------------------------
!
! MODULE: module_utm2deg_bc
!
!> @author
!> Marc T. Henry de Frahan, National Renewable Energy Laboratory
!
!> @brief
!> Convert arrays of UTM coordinates into Lat/Lon arrays (and vice-versa)
!
!> @details

!> This is a Fortran version of the <a
!> href="https://www.mathworks.com/matlabcentral/fileexchange/10914-utm2deg">utm2deg</a>
!> and <a
!> href="https://www.mathworks.com/matlabcentral/fileexchange/10915-deg2utm">deg2utm</a>
!> matlab functions.  All credit goes to Rafael Palacios and Gabriel Ruiz Martinez.
!
! 
!------------------------------------------------------------------------------


MODULE module_utmdeg_converter
  use module_constants
  use module_str2int
  implicit none

  double precision, parameter :: sa = 6378137.000000, &
       sb = 6356752.314245, &
       e2 = ( sqrt(sa**2 - sb**2) )  / sb, &
       e22 = e2**2, &
       c = (sa/sb) * sa, &
       alpha = 0.75 * e22, &
       beta  = 5.d0/3.d0 * (alpha**2), &
       gamma = 35.d0/27.d0 * (alpha**3)

CONTAINS

  
  !--------------------------------------------------------------------------------
  !
  !> @brief Convert UTM coordinates to Lat/Lon
  !
  !> @param[in] n       number of points
  !> @param[in] x       UTM x coordinate
  !> @param[in] y       UTM y coordinate
  !> @param[in] utmzone UTM time zone
  !> @param[out] lat    Latitude   +ddd.ddddd  WGS84
  !> @param[out] lon    Longitude  +ddd.ddddd  WGS84
  !
  !> @details A Fortran version of href="https://www.mathworks.com/matlabcentral/fileexchange/10914-utm2deg">utm2deg</a>
  !>
  !> Example 1:
  !>   integer, dimension(6) :: test_x = (/458731,  407653,  239027,  230253,  343898,  362850/)
  !>   integer, dimension(6) :: test_y = (/4462881, 5126290, 4163083, 3171843, 4302285, 2772478/)
  !>   character(len=4), dimension(6) :: test_utmzone = (/'30 T', '32 T', '11 S', '28 R', '15 S', '51 R'/)
  !>   real, dimension(6) :: test_lat, test_lon
  !>   call utm2deg(6,test_x,test_y,test_utmzone,test_lat,test_lon)
  !>   write(*,*)"lat:"  test_lat
  !>   write(*,*)"lon:" test_lon
  !> Result:
  !>   lat: 40.315430   46.283902   37.577834   28.645647   38.855552   25.061780
  !>   lon: -3.485713    7.801235 -119.955246  -17.759537  -94.799019  121.640266
  !
  !--------------------------------------------------------------------------------
  SUBROUTINE utm2deg (n, x, y, utmzone, lat, lon)

    !================================================================================
    ! Initialize
    integer,      intent(in)  :: n
    integer,      dimension(n), intent(in)  :: x, y
    character(4), dimension(n), intent(in)  :: utmzone
    real,         dimension(n), intent(out) :: lat, lon

    ! internal
    integer :: i
    integer,          dimension(n) :: xc, yc, zonenum, stat
    logical,          dimension(n) :: mask
    character(1),     dimension(n) :: zonestr, hemis
    double precision, dimension(n) :: S, latnum,  v,  a, &
         a1,  a2, j2, j4, j6, Bm, &
         b, Epsi, Eps, nab, senoheps, Delt, TaO

    !================================================================================
    ! Conversion
    
    zonestr = utmzone(:)(4:4)
    mask(:)  = zonestr(:) .gt. 'M'
    hemis(:) = merge("N","S",mask(:))
    
    call str2int(utmzone(:)(1:2), zonenum, stat)

    xc(:) = x(:) - 500000
    mask(:)  = (hemis(:) .eq. "S") .or. (hemis(:) .eq. "s")
    yc(:) = merge( y(:) - 10000000, y(:), mask(:))
    
    S(:) = ( ( zonenum(:) * 6.d0 ) - 183.d0 )
    latnum(:) = dble(yc(:)) / ( 6366197.724 * 0.9996 )
    v(:) = ( c / sqrt( 1.d0 +  e22 * cos(latnum(:))**2 ) ) * 0.9996
    a(:) = xc(:)/v(:)
    a1(:) = sin(2.d0*latnum(:))
    a2(:) = a1(:) * (cos(latnum(:))**2)
    j2(:) = latnum(:) + 0.5*a1(:)
    j4(:) = 0.25 * ( 3.d0 * j2(:) + a2(:) )
    j6(:) = ( 5.d0 * j4(:)  +  a2(:) * ( cos(latnum(:))**2) ) / 3.0
    Bm(:) = 0.9996 * c * ( latnum(:) - alpha*j2(:) + beta*j4(:) - gamma*j6(:) )
    b(:) = ( yc(:) - Bm(:) ) / v(:)
    Epsi(:) =  0.5 * ( e22 * (a(:)**2) ) * ( cos(latnum(:))**2)
    Eps(:) = a(:) * ( 1 -  Epsi(:)/3.d0  )
    nab(:) = b(:) * ( 1 - Epsi(:) ) + latnum(:)
    senoheps(:) = 0.5 * ( exp(Eps(:)) - exp(-Eps(:)) )
    Delt(:) = atan(senoheps(:) / cos(nab(:)) )
    TaO(:) = atan(cos(Delt(:)) * tan(nab(:)))

    lat = ( latnum(:) + ( 1.d0 + e22* (cos(latnum(:))**2) - ( 3.d0 / 2.d0 ) * e22 * sin(latnum(:)) * cos(latnum(:)) * ( TaO(:) - latnum(:) ) ) &
         * ( TaO(:) - latnum(:) ) ) *  (180.d0 / pi)
    lon = Delt(:) *(180.d0 / pi )  + S

  END SUBROUTINE utm2deg

  !--------------------------------------------------------------------------------
  !
  !> @brief Convert Lat/Lon coordinates to UTM
  !
  !> @param[in] n       number of points
  !> @param[in] x       UTM x coordinate
  !> @param[in] y       UTM y coordinate
  !> @param[in] utmzone UTM time zone
  !> @param[out] lat    Latitude   +ddd.ddddd  WGS84
  !> @param[out] lon    Longitude  +ddd.ddddd  WGS84
  !
  !> @details A Fortran version of href="https://www.mathworks.com/matlabcentral/fileexchange/10915-deg2utm">deg2utm</a>
  !>
  !> Example 1:
  
  !>   integer, dimension(6) :: test_x = (/458731,  407653,  239027,  230253,  343898,  362850/)
  !>   integer, dimension(6) :: test_y = (/4462881, 5126290, 4163083, 3171843, 4302285, 2772478/)
  !>   character(len=4), dimension(6) :: test_utmzone = (/'30 T', '32 T', '11 S', '28 R', '15 S', '51 R'/)
  !>   real, dimension(6) :: test_lat, test_lon
  !>   call utm2deg(6,test_x,test_y,test_utmzone,test_lat,test_lon)
  !>   call deg2utm(6,test_x,test_y,test_utmzone,test_lat,test_lon)
  !>   write(*,*)"x:" test_x
  !>   write(*,*)"y:" test_y
  !>   write(*,*)"utmzone:" test_utmzone
  !> Result:
  !>   x: 458731,  407653,  239027,  230253,  343898,  362850
  !>   y: 4462881, 5126290, 4163083, 3171843, 4302285, 2772478
  !>   utmzone: '30 T', '32 T', '11 S', '28 R', '15 S', '51 R'
  !
  !--------------------------------------------------------------------------------
  SUBROUTINE deg2utm (n, lat, lon, x, y, utmzone)

    !================================================================================
    ! Initialize
    integer,      intent(in)  :: n
    real,         dimension(n), intent(in) :: lat, lon
    integer,      dimension(n), intent(out)  :: x, y
    character(4), dimension(n), intent(out)  :: utmzone

    ! internal
    integer :: i
    integer,          dimension(n) :: zonenum
    logical,          dimension(n) :: mask
    character(1),     dimension(n) :: zonestr
    double precision, dimension(n) :: S, latnum,  v,  a, &
         a1,  a2, j2, j4, j6, Bm, &
         epsilon, nu,  ta, deltaS

    
    !================================================================================
    ! Conversion
    latnum(:) = lat(:) * pi/180.d0

    zonenum(:) = floor( lon(:) / 6.d0  + 31.d0)
    S(:) =   zonenum(:) * 6.d0 - 183.d0
    deltaS(:) = (lon(:) -  S(:)) * ( pi / 180.d0 )

    do i = 1, n
       if (lat(i)<-72) then
          zonestr(i)='C'
       else if (lat(i)<-64) then
          zonestr(i)='D'
       else if (lat(i)<-56) then
          zonestr(i)='E'
       else if (lat(i)<-48) then
          zonestr(i)='F'
       else if (lat(i)<-40) then
          zonestr(i)='G'
       else if (lat(i)<-32) then
          zonestr(i)='H'
       else if (lat(i)<-24) then
          zonestr(i)='J'
       else if (lat(i)<-16) then
          zonestr(i)='K'
       else if (lat(i)<-8) then
          zonestr(i)='L'
       else if (lat(i)<0) then
          zonestr(i)='M'
       else if (lat(i)<8) then
          zonestr(i)='N'
       else if (lat(i)<16) then
          zonestr(i)='P'
       else if (lat(i)<24) then
          zonestr(i)='Q'
       else if (lat(i)<32) then
          zonestr(i)='R'
       else if (lat(i)<40) then
          zonestr(i)='S'
       else if (lat(i)<48) then
          zonestr(i)='T'
       else if (lat(i)<56) then
          zonestr(i)='U'
       else if (lat(i)<64) then
          zonestr(i)='V'
       else if (lat(i)<72) then
          zonestr(i)='W'
       else
          zonestr(i)='X'
       endif
    enddo

    a(:) = cos(latnum(:)) * sin(deltaS(:))
    epsilon(:) = 0.5 * log( (1 +  a(:)) / (1 - a(:)) )
    nu(:) = atan( tan(latnum(:)) / cos(deltaS(:)) ) - latnum(:)
    v(:) = c / ( sqrt( 1.d0 + ( e22 * cos(latnum(:))**2 ) ) ) * 0.9996
    ta(:) = 0.5 * e22 * (epsilon(:) * cos(latnum(:)))**2
    a1(:) = sin(2.d0 * latnum(:))
    a2(:) = a1(:) * cos(latnum(:))**2
    j2(:) = latnum(:) + 0.5*a1(:)
    j4(:) = 0.25 * (3.d0 * j2(:) + a2(:))
    j6(:) = ( 5.d0 * j4(:) + a2(:) * cos(latnum(:))**2 ) / 3.d0
    Bm(:) = 0.9996 * c * ( latnum(:) - alpha * j2(:) + beta * j4(:) - gamma * j6(:) )

    ! UTM coordinates
    x(:) =  nint(epsilon(:) * v(:) * ( 1 + ( ta(:) / 3.d0 ) ) + 500000)
    y(:) =  nint(nu(:) * v(:) * ( 1 + ta(:) ) + Bm(:))
    mask(:)  = y(:) .lt. 0
    y(:) = merge(9999999+y(:), y(:) , mask(:))
    do i = 1, n
       write(utmzone(i),'(I2,A2)')zonenum(i),zonestr(i)
    enddo
    
  END SUBROUTINE deg2utm
  
END MODULE module_utmdeg_converter
