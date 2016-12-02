!------------------------------------------------------------------------------
!
! MODULE: module_openfoam_bc
!
!> @author
!> J. Michalakes and M. Churchfield, National Renewable Energy Laboratory
!
!> @brief
!> Reads openfoam data
!
!> @date 01/12/2016 J. Michalakes and M. Churchfield
!> - Initial version from WRFTOOOF
! 
!------------------------------------------------------------------------------


MODULE module_openfoam_bc
  USE module_dm
  IMPLICIT NONE
  TYPE of_pt_t
     DOUBLE PRECISION lat,lon,lz
     INTEGER i,j
  END TYPE of_pt_t
  TYPE bdy_t
     INTEGER npoints
     TYPE (of_pt_t), ALLOCATABLE :: point(:)
  END TYPE bdy_t
  INTEGER, PARAMETER :: nbdys = 7
  TYPE (bdy_t) bdy(nbdys)
  INTEGER, PARAMETER :: BDY_XS = 1, BDY_XE = 2, &
       BDY_YS = 3, BDY_YE = 4, &
       BDY_ZS = 5, BDY_ZE = 6, &
       INTERIOR = 7 ! not a boundary per se but we can handle this too

CONTAINS

  
  !--------------------------------------------------------------------------------
  !
  !> @brief Output an error message and stop the code
  !> @param message messsage to output
  !
  !--------------------------------------------------------------------------------
  SUBROUTINE wrf_error_fatal (message)
    CHARACTER*(*) :: message
    WRITE(0,*)TRIM(message)
    STOP 999
  END SUBROUTINE wrf_error_fatal

  
  !--------------------------------------------------------------------------------
  !
  !> @brief Determine OpenFoam coordinates for the BCs
  !> @param whichbdy  body index 
  !> @param fname     OpenFoam BC file data
  !
  !--------------------------------------------------------------------------------
  SUBROUTINE read_openfoam_bdy_coords( whichbdy , fname )
    INTEGER, INTENT(IN) :: whichbdy
    CHARACTER*(*), INTENT(IN) :: fname
    !local
    LOGICAL , EXTERNAL :: wrf_dm_on_monitor
    CHARACTER*256 message,latstr,lonstr
    INTEGER ipoint,ibdy, ierr, i
    !exec

    ierr = 0
    IF ( wrf_dm_on_monitor()) THEN
       OPEN ( 75, file=TRIM(fname), form="formatted", status="old", err=2222 )
       ! count up the number of lines in the file
       bdy(whichbdy)%npoints = 0
       DO WHILE ( .TRUE. )
          READ(75,*,END=2210)
          bdy(whichbdy)%npoints = bdy(whichbdy)%npoints + 1
       ENDDO
2210   CONTINUE
       CLOSE (75)
       ALLOCATE(bdy(whichbdy)%point(bdy(whichbdy)%npoints))
       ! now read them in
       OPEN ( 75, file=TRIM(fname), form="formatted", status="old", err=2222 )
       DO ipoint = 1,bdy(whichbdy)%npoints
          READ(75,*,ERR=2222)latstr, lonstr, bdy(whichbdy)%point(ipoint)%lz
          i = INDEX(latstr,'N')
          if ( i .NE. 0 ) latstr(i:i) = ' '
          i = INDEX(latstr,'S')
          if ( i .NE. 0 ) latstr(i:i) = '-'
          READ(latstr,*)bdy(whichbdy)%point(ipoint)%lat
          i = INDEX(lonstr,'E')
          if ( i .NE. 0 ) lonstr(i:i) = ' '
          i = INDEX(lonstr,'W')
          if ( i .NE. 0 ) lonstr(i:i) = '-'
          READ(lonstr,*)bdy(whichbdy)%point(ipoint)%lon
       ENDDO
       GOTO 2220
2222   CONTINUE
       ierr = 1
2220   CONTINUE
    ENDIF
    CALL wrf_dm_bcast_integer(ierr,1)
    IF ( ierr .NE. 0 ) THEN
       WRITE(message,*)'read_openfoam_bdy_coords: some error reading in bdy coords from ',TRIM(fname)
    ENDIF
    CALL wrf_dm_bcast_integer(bdy(whichbdy)%npoints,1)
    DO ipoint = 1,bdy(whichbdy)%npoints
       CALL wrf_dm_bcast_double(bdy(whichbdy)%point(ipoint)%lat)
       CALL wrf_dm_bcast_double(bdy(whichbdy)%point(ipoint)%lon)
       CALL wrf_dm_bcast_double(bdy(whichbdy)%point(ipoint)%lz)
    ENDDO
    RETURN
  END SUBROUTINE read_openfoam_bdy_coords

  
  !--------------------------------------------------------------------------------
  !
  !> @brief Precompute the OpenFoam coordinate points
  !> @param ibdy      body index
  !> @param xlat      latitude
  !> @param xlon      longitude
  !> @param ids       start i index
  !> @param ide       end i index   
  !> @param jds       start j index 
  !> @param jde       end j index   
  !> @param kds       start k index 
  !> @param kde       end k index   
  !> @param ips       start i index
  !> @param ipe       end i index   
  !> @param jps       start j index 
  !> @param jpe       end j index   
  !> @param ims       start i index
  !> @param ime       end i index   
  !> @param jms       start j index 
  !> @param jme       end j index
  !
  !--------------------------------------------------------------------------------
  SUBROUTINE precompute_openfoam_points(ibdy,xlat,xlong,ids,ide,jds,jde,ips,ipe,jps,jpe,ims,ime,jms,jme )
    INTEGER, INTENT(IN) :: ibdy,ids,ide,jds,jde,ips,ipe,jps,jpe,ims,ime,jms,jme
    REAL, DIMENSION(ids:ide-1,jds:jde-1), INTENT(IN) :: xlat, xlong
    ! local
    INTEGER i,j,ipoint,idummy
    REAL of_lat, of_lon
    REAL dsw,dse,dnw,dne,lim,dmin
    CHARACTER*256 bdy_cache_name, message
    LOGICAL incache

    IF ( ALLOCATED(bdy(ibdy)%point) ) THEN
       incache=.TRUE.
       WRITE(bdy_cache_name,'("bdy_cache_",I1)')ibdy
       OPEN(75,file=TRIM(bdy_cache_name),form="formatted",status="old",ERR=9911)
       GOTO 9910
9911   CONTINUE
       OPEN(75,file=TRIM(bdy_cache_name),form="formatted",status="new",ERR=9911)
       incache=.FALSE.
9910   CONTINUE
       DO ipoint = 1,bdy(ibdy)%npoints
          IF ( incache ) THEN
             READ(75,*)idummy,bdy(ibdy)%point(ipoint)%i,bdy(ibdy)%point(ipoint)%j
             IF ( idummy .NE. ipoint ) THEN
                WRITE(message,*)'problem reading: ',TRIM(bdy_cache_name),': ',idummy,' ne ',ipoint
                CALL wrf_error_fatal(message)
             ENDIF
          ELSE
             of_lat = bdy(ibdy)%point(ipoint)%lat
             of_lon = bdy(ibdy)%point(ipoint)%lon
             dmin = 999999.9
             DO j = jps,min(jpe,jde-2)
                DO i = ips,min(ipe,ide-2)
                   ! ignore special case where of point lies outside the grid of cell centers
                   ! should not put OF grid that close to a WRF boundary
                   ! also note the cavalier way we ignore curvature and assume the
                   ! grid cells are perfectly square and that lat and lon are Cartesian
                   dsw = sqrt((of_lat-xlat(i ,j ))*(of_lat-xlat(i ,j )) + (of_lon-xlong(i ,j ))*(of_lon-xlong(i ,j )))
                   !!absolute closest
                   !IF ( dsw .LT. dmin ) THEN
                   !alternate scheme, pick the point that is closest to the sw of the openfoam point
                   IF ( dsw .LT. dmin .AND. of_lat .GE. xlat(i,j) .AND. of_lon .GE. xlong(i,j) ) THEN
                      bdy(ibdy)%point(ipoint)%i = i
                      bdy(ibdy)%point(ipoint)%j = j
                      dmin = dsw
                   ENDIF
                ENDDO
             ENDDO
             WRITE(75,*)ipoint,bdy(ibdy)%point(ipoint)%i,bdy(ibdy)%point(ipoint)%j
          ENDIF
       ENDDO
       CLOSE(75)
    ENDIF
  END SUBROUTINE precompute_openfoam_points

  
  !--------------------------------------------------------------------------------
  !
  !> @brief Calculate the rotation angle of the WRF w.r.t. true
  !> @param xlat      latitude
  !> @param dx        mesh spacing
  !> @param ids       start i index
  !> @param ide       end i index   
  !> @param jds       start j index 
  !> @param jde       end j index   
  !> @param kds       start k index 
  !> @param kde       end k index   
  !> @param ips       start i index
  !> @param ipe       end i index   
  !> @param jps       start j index 
  !> @param jpe       end j index   
  !> @param ims       start i index
  !> @param ime       end i index   
  !> @param jms       start j index 
  !> @param jme       end j index
  !
  !--------------------------------------------------------------------------------
  REAL FUNCTION rotation_angle( xlat,dx,ids,ide,jds,jde,ips,ipe,jps,jpe,ims,ime,jms,jme )
    INTEGER, INTENT(IN) :: ids,ide,jds,jde,ips,ipe,jps,jpe,ims,ime,jms,jme
    REAL, DIMENSION(ims:ime-1,jms:jme-1), INTENT(IN) :: xlat
    REAL, INTENT(IN) :: dx
    !local
    REAL cen_lat_west, cen_lat_east, dlat, dist, domlen
    cen_lat_west = -9999.
    cen_lat_east = -9999.
    IF ( jps .LE. (jde-jds)/2 .AND. (jde-jds)/2 .LT. jpe ) THEN
       IF ( ips .EQ. ids ) cen_lat_west = xlat(ips,(jde-jds)/2)
       IF ( ipe .EQ. ide ) cen_lat_east = xlat(ipe-1,(jde-jds)/2)
    ENDIF
    cen_lat_west = wrf_dm_max_real( cen_lat_west )
    cen_lat_east = wrf_dm_max_real( cen_lat_east )
    dlat = (cen_lat_west-cen_lat_east)/360.
    dist = (dlat * 40000000)
    domlen = ( ide-ids )*dx
    rotation_angle = asin( dist / domlen )
  END FUNCTION rotation_angle

END MODULE module_openfoam_bc
