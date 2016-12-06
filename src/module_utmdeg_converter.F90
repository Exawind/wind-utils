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
  use module_str2int
  implicit none
  

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

    ! Initialize
    integer,      intent(in)  :: n
    integer,      dimension(n), intent(in)  :: x, y
    character(4), dimension(n), intent(in)  :: utmzone
    real,         dimension(n), intent(out) :: lat, lon

    ! internal
    integer :: i
    integer,      dimension(n) :: xc, yc, zonenum, stat
    logical,      dimension(n) :: mask
    character(1), dimension(n) :: zonestr, hemis
    real,         dimension(n) :: S, latnum, v, a, a1, a2, j2, j4, j6, &
         alpha, beta, gamma, Bm, b, Epsi, Eps, nab, &
         senoheps, Delt, TaO
    real, parameter :: sa = 6378137.000000, sb = 6356752.314245
    real, parameter :: e2 = ( sqrt(sa**2 - sb**2) )  / sb, &
         e22 = e2**2, &
         c = (sa**2) / sb
    double precision, parameter :: pi=4.D0*datan(1.D0)


    ! Get hemisphere
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
    alpha(:) = 0.75 * e22
    beta(:)  = 5.d0/3.d0 * (alpha(:)**2)
    gamma(:) = 35.d0/27.d0 * (alpha(:)**3)
    Bm(:) = 0.9996 * c * ( latnum(:) - alpha(:)*j2(:) + beta(:)*j4(:) - gamma(:)*j6(:) )
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


  
END MODULE module_utmdeg_converter
