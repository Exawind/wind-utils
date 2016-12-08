!------------------------------------------------------------------------------
!
! MODULE: module_exodus_bc
!
!> @author
!> M. T. Henry de Frahan, National Renewable Energy Laboratory
!
!> @brief
!> Reads Exodus data
!
!> @date 01/12/2016 J. Michalakes and M. Churchfield
!> - Initial version from WRFTOOOF
! 
!------------------------------------------------------------------------------


module module_exodus
  use module_dm
  use module_ncderrcheck
  use module_utmdeg_converter
  implicit none

  type bdy_t
     integer num_nodes
     integer num_nod_var_dimid
     integer ofid, time_whole_id, coordx_id, coordy_id, coordz_id, &
          vals_nod_var1_id, vals_nod_var2_id, vals_nod_var3_id, &
          vals_nod_var4_id, vals_nod_var5_id, vals_nod_var6_id, &
          vals_nod_var7_id     
     integer,      dimension(:),   allocatable :: exo_wrf_i, exo_wrf_j
     real,         dimension(:),   allocatable :: coordx, coordy, coordz, lat, lon
     real,         dimension(:,:), allocatable ::vals_nod_var1, vals_nod_var2, vals_nod_var3, &
          vals_nod_var4, vals_nod_var5, vals_nod_var6, vals_nod_var7     

     character(4), dimension(:), allocatable :: utmzone
  end type bdy_t

  integer, parameter :: nbdys = 7
  type (bdy_t) bdy(nbdys)
  integer, parameter :: BDY_XS = 1, BDY_XE = 2, &
       BDY_YS = 3, BDY_YE = 4, &
       BDY_ZS = 5, BDY_ZE = 6, &
       INTERIOR = 7 ! not a boundary per se but we can handle this too

contains


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
  !> @brief Read in Exodus coordinates from a mesh
  !> @param bdynum    body index 
  !> @param ofname    Exodus BC file data
  !> @param ntimes    Number of time steps
  !> @param times     Time steps
  !
  !--------------------------------------------------------------------------------
  subroutine prep_exodus( bdynum, ofname, ntimes, times)

    ! initialize
    integer,           intent(in) :: bdynum, ntimes
    character*(*),     intent(in) :: ofname
    character(len=19), intent(in) :: times(100)

    ! internal
    include 'netcdf.inc'
    integer stat, rec

    ! Store time information
    integer,dimension(8) :: date_values
    character(12)  :: date_str
    character(2) :: dd
    character(2) :: mm
    character(4) :: yyyy

    ! Variable ids and values
    integer num_nodes_dimid, time_step_dimid, len_name_dimid
    integer info_records_id, eb_names_id,  &
         vals_nod_var1_dims(2), vals_nod_var2_dims(2), &
         vals_nod_var3_dims(2), vals_nod_var4_dims(2), &
         vals_nod_var5_dims(2), vals_nod_var6_dims(2), &
         vals_nod_var7_dims(2), &
         name_nod_var_id, name_nod_var_dims(2)

    character(len=255) :: info_records
    integer :: info_records_start(2), info_records_count(2)

    character(len=255) :: eb_names

    character(len=33), dimension(7) :: name_nod_var
    integer :: name_nod_var_start(2), name_nod_var_count(2)

    double precision time_whole(ntimes)
    integer :: time_whole_start(2), time_whole_count(2)
    
    WRITE(0,*)'opening output file ',TRIM(ofname)
    stat = NF_OPEN(ofname, NF_WRITE, bdy(bdynum)%ofid )
    call ncderrcheck( __LINE__ ,stat )

    !================================================================================
    ! Get system information
    call date_and_time(VALUES=date_values)
    write(  dd,'(i2)') date_values(3)
    write(  mm,'(i2)') date_values(2)
    write(yyyy,'(i4)') date_values(1)
    write(date_str,*),dd,'/',mm,'/',yyyy

    !================================================================================
    ! Get mesh information
    stat = NF_INQ_DIMID(bdy(bdynum)%ofid,"num_nodes", num_nodes_dimid)
    call ncderrcheck( __LINE__,stat)
    stat = NF_INQ_DIMLEN(bdy(bdynum)%ofid, num_nodes_dimid, bdy(bdynum)%num_nodes)
    call ncderrcheck( __LINE__,stat)

    stat = NF_INQ_DIMID(bdy(bdynum)%ofid,"time_step", time_step_dimid)
    call ncderrcheck( __LINE__,stat)

    stat = NF_INQ_DIMID(bdy(bdynum)%ofid,"len_name", len_name_dimid)
    call ncderrcheck( __LINE__,stat)

    stat = NF_INQ_VARID(bdy(bdynum)%ofid,"info_records", info_records_id)
    call ncderrcheck( __LINE__,stat)

    stat = NF_INQ_VARID(bdy(bdynum)%ofid,"eb_names", eb_names_id)
    call ncderrcheck( __LINE__,stat)

    stat = NF_INQ_VARID(bdy(bdynum)%ofid,"time_whole", bdy(bdynum)%time_whole_id)
    call ncderrcheck( __LINE__,stat)

    stat = NF_INQ_VARID(bdy(bdynum)%ofid,"coordx", bdy(bdynum)%coordx_id)
    call ncderrcheck( __LINE__,stat)

    stat = NF_INQ_VARID(bdy(bdynum)%ofid,"coordy", bdy(bdynum)%coordy_id)
    call ncderrcheck( __LINE__,stat)

    stat = NF_INQ_VARID(bdy(bdynum)%ofid,"coordz", bdy(bdynum)%coordz_id)
    call ncderrcheck( __LINE__,stat)

    !================================================================================
    ! Define new dimensions and variables

    ! put in define mode
    stat = NF_REDEF(bdy(bdynum)%ofid)
    call ncderrcheck( __LINE__,stat)

    ! Write out num_nod_var dimension
    stat = nf_def_dim(bdy(bdynum)%ofid, "num_nod_var", 7 , bdy(bdynum)%num_nod_var_dimid) ! CHANGE
    call ncderrcheck( __LINE__,stat)

    ! vals_nod_var1 variable
    vals_nod_var1_dims(1) = num_nodes_dimid
    vals_nod_var1_dims(2) = time_step_dimid
    stat = nf_def_var(bdy(bdynum)%ofid, "vals_nod_var1", nf_double, 2, vals_nod_var1_dims , bdy(bdynum)%vals_nod_var1_id)
    call ncderrcheck( __LINE__,stat)

    ! vals_nod_var2 variable
    vals_nod_var2_dims(1) = num_nodes_dimid
    vals_nod_var2_dims(2) = time_step_dimid
    stat = nf_def_var(bdy(bdynum)%ofid, "vals_nod_var2", nf_double, 2, vals_nod_var2_dims , bdy(bdynum)%vals_nod_var2_id)
    call ncderrcheck( __LINE__,stat)

    ! vals_nod_var3 variable
    vals_nod_var3_dims(1) = num_nodes_dimid
    vals_nod_var3_dims(2) = time_step_dimid
    stat = nf_def_var(bdy(bdynum)%ofid, "vals_nod_var3", nf_double, 2, vals_nod_var3_dims , bdy(bdynum)%vals_nod_var3_id)
    call ncderrcheck( __LINE__,stat)

    ! vals_nod_var4 variable
    vals_nod_var4_dims(1) = num_nodes_dimid
    vals_nod_var4_dims(2) = time_step_dimid
    stat = nf_def_var(bdy(bdynum)%ofid, "vals_nod_var4", nf_double, 2, vals_nod_var4_dims , bdy(bdynum)%vals_nod_var4_id)
    call ncderrcheck( __LINE__,stat)

    ! vals_nod_var5 variable
    vals_nod_var5_dims(1) = num_nodes_dimid
    vals_nod_var5_dims(2) = time_step_dimid
    stat = nf_def_var(bdy(bdynum)%ofid, "vals_nod_var5", nf_double, 2, vals_nod_var5_dims , bdy(bdynum)%vals_nod_var5_id)
    call ncderrcheck( __LINE__,stat)

    ! vals_nod_var6 variable
    vals_nod_var6_dims(1) = num_nodes_dimid
    vals_nod_var6_dims(2) = time_step_dimid
    stat = nf_def_var(bdy(bdynum)%ofid, "vals_nod_var6", nf_double, 2, vals_nod_var6_dims , bdy(bdynum)%vals_nod_var6_id)
    call ncderrcheck( __LINE__,stat)

    ! vals_nod_var7 variable
    vals_nod_var7_dims(1) = num_nodes_dimid
    vals_nod_var7_dims(2) = time_step_dimid
    stat = nf_def_var(bdy(bdynum)%ofid, "vals_nod_var7", nf_double, 2, vals_nod_var7_dims , bdy(bdynum)%vals_nod_var7_id)
    call ncderrcheck( __LINE__,stat)

    ! name_nod_var variable
    name_nod_var_dims(1) = len_name_dimid
    name_nod_var_dims(2) = bdy(bdynum)%num_nod_var_dimid
    stat = nf_def_var(bdy(bdynum)%ofid, "name_nod_var", nf_char, 2, name_nod_var_dims , name_nod_var_id)
    call ncderrcheck( __LINE__,stat)

    ! End define mode
    stat = NF_ENDDEF(bdy(bdynum)%ofid)
    call ncderrcheck( __LINE__,stat)

    !================================================================================
    ! Fill in some of the variables

    ! time_whole variable
    do rec = 1, ntimes
       time_whole(rec) = 4.0 + dble(rec)
    end do
    time_whole_count(1) = ntimes ! number of doubles to write
    time_whole_count(2) = 1      ! only write one record
    time_whole_start(1) = 1      ! start at beginning of variable
    time_whole_start(2) = 1      ! record number to write                        
    stat = nf_put_vara_double(bdy(bdynum)%ofid, bdy(bdynum)%time_whole_id, time_whole_start, time_whole_count, time_whole)
    call ncderrcheck( __LINE__ ,stat )
   
    ! info_records variable
    info_records = "Made with WRFTONALU on"//trim(date_str)
    info_records_start(1) = 1                              ! start at beginning of variable
    info_records_start(2) = 1                              ! record number to write
    info_records_count(1) = len(trim(info_records))        ! number of chars to write
    info_records_count(2) = 1                              ! only write one record
    stat = nf_put_vara_text(bdy(bdynum)%ofid, info_records_id, info_records_start, info_records_count, info_records)
    call ncderrcheck( __LINE__ ,stat )

    ! eb_names
    eb_names =  "block_101"
    stat = nf_put_var_text(bdy(bdynum)%ofid, eb_names_id, trim(eb_names))
    call ncderrcheck( __LINE__ ,stat )

    ! name_nod_var
    name_nod_var(1) = 'cont_velocity_bc_x'
    name_nod_var(2) = 'cont_velocity_bc_y'
    name_nod_var(3) = 'cont_velocity_bc_z'
    name_nod_var(4) = 'temperature_bc'
    name_nod_var(5) = 'velocity_bc_x'
    name_nod_var(6) = 'velocity_bc_y'
    name_nod_var(7) = 'velocity_bc_z'  
    name_nod_var_start(1) = 1                              ! start at beginning of variable
    name_nod_var_start(2) = 1                              ! record number to write
    name_nod_var_count(1) = LEN(trim(name_nod_var(1)))     ! number of chars to write
    name_nod_var_count(2) = 1                              ! only write one record
    stat = nf_put_vara_text(bdy(bdynum)%ofid, name_nod_var_id, name_nod_var_start, name_nod_var_count, name_nod_var(1))
    call ncderrcheck( __LINE__ ,stat )
    name_nod_var_start(1) = 1
    name_nod_var_start(2) = 2
    name_nod_var_count(1) = LEN(trim(name_nod_var(2)))
    name_nod_var_count(2) = 1
    stat = nf_put_vara_text(bdy(bdynum)%ofid, name_nod_var_id, name_nod_var_start, name_nod_var_count, name_nod_var(2))
    call ncderrcheck( __LINE__ ,stat )
    name_nod_var_start(1) = 1
    name_nod_var_start(2) = 3
    name_nod_var_count(1) = LEN(trim(name_nod_var(3)))
    name_nod_var_count(2) = 1
    stat = nf_put_vara_text(bdy(bdynum)%ofid, name_nod_var_id, name_nod_var_start, name_nod_var_count, name_nod_var(3))
    call ncderrcheck( __LINE__ ,stat )
    name_nod_var_start(1) = 1
    name_nod_var_start(2) = 4
    name_nod_var_count(1) = LEN(trim(name_nod_var(4)))
    name_nod_var_count(2) = 1
    stat = nf_put_vara_text(bdy(bdynum)%ofid, name_nod_var_id, name_nod_var_start, name_nod_var_count, name_nod_var(4))
    call ncderrcheck( __LINE__ ,stat )
    name_nod_var_start(1) = 1
    name_nod_var_start(2) = 5
    name_nod_var_count(1) = LEN(trim(name_nod_var(5)))
    name_nod_var_count(2) = 1
    stat = nf_put_vara_text(bdy(bdynum)%ofid, name_nod_var_id, name_nod_var_start, name_nod_var_count, name_nod_var(5))
    call ncderrcheck( __LINE__ ,stat )
    name_nod_var_start(1) = 1
    name_nod_var_start(2) = 6
    name_nod_var_count(1) = LEN(trim(name_nod_var(6)))
    name_nod_var_count(2) = 1
    stat = nf_put_vara_text(bdy(bdynum)%ofid, name_nod_var_id, name_nod_var_start, name_nod_var_count, name_nod_var(6))
    call ncderrcheck( __LINE__ ,stat )
    name_nod_var_start(1) = 1
    name_nod_var_start(2) = 7
    name_nod_var_count(1) = LEN(trim(name_nod_var(7)))
    name_nod_var_count(2) = 1
    stat = nf_put_vara_text(bdy(bdynum)%ofid, name_nod_var_id, name_nod_var_start, name_nod_var_count, name_nod_var(7))
    call ncderrcheck( __LINE__ ,stat )


    !================================================================================
    ! Allocate variables while we are at it
    allocate( bdy(bdynum)%coordx(bdy(bdynum)%num_nodes))
    allocate( bdy(bdynum)%coordy(bdy(bdynum)%num_nodes))
    allocate( bdy(bdynum)%coordz(bdy(bdynum)%num_nodes))
    allocate( bdy(bdynum)%utmzone(bdy(bdynum)%num_nodes))
    allocate( bdy(bdynum)%lat(bdy(bdynum)%num_nodes))
    allocate( bdy(bdynum)%lon(bdy(bdynum)%num_nodes))

    allocate( bdy(bdynum)%vals_nod_var1(bdy(bdynum)%num_nodes, ntimes))
    allocate( bdy(bdynum)%vals_nod_var2(bdy(bdynum)%num_nodes, ntimes))
    allocate( bdy(bdynum)%vals_nod_var3(bdy(bdynum)%num_nodes, ntimes))
    allocate( bdy(bdynum)%vals_nod_var4(bdy(bdynum)%num_nodes, ntimes))
    allocate( bdy(bdynum)%vals_nod_var5(bdy(bdynum)%num_nodes, ntimes))
    allocate( bdy(bdynum)%vals_nod_var6(bdy(bdynum)%num_nodes, ntimes))
    allocate( bdy(bdynum)%vals_nod_var7(bdy(bdynum)%num_nodes, ntimes))

    bdy(bdynum)%coordx  = 0.0
    bdy(bdynum)%coordy  = 0.0
    bdy(bdynum)%coordz  = 0.0
    bdy(bdynum)%utmzone = ''
    bdy(bdynum)%lat     = 0.0
    bdy(bdynum)%lon     = 0.0

    bdy(bdynum)%vals_nod_var1 = 0.0
    bdy(bdynum)%vals_nod_var2 = 0.0
    bdy(bdynum)%vals_nod_var3 = 0.0
    bdy(bdynum)%vals_nod_var4 = 0.0
    bdy(bdynum)%vals_nod_var5 = 0.0
    bdy(bdynum)%vals_nod_var6 = 0.0
    bdy(bdynum)%vals_nod_var7 = 0.0    
    
  end subroutine  prep_exodus


  !--------------------------------------------------------------------------------
  !
  !> @brief Close netCDF file
  !> @param bdynum  body index 
  !
  !--------------------------------------------------------------------------------
  subroutine close_exodus( bdynum )

    ! initialize
    integer,       intent(in) :: bdynum

    ! internal
    include 'netcdf.inc'
    integer stat

    write(*,*)'closing output file'
    
    ! Close the file
    stat = NF_CLOSE(bdy(bdynum)%ofid);
    call ncderrcheck( __LINE__ ,stat )

  end subroutine close_exodus


  !--------------------------------------------------------------------------------
  !
  !> @brief Read in Exodus coordinates from a mesh
  !> @param bdynum        body index
  !> @param lat_offset    latitude offset
  !> @param lon_offset    longitude offset
  !
  !--------------------------------------------------------------------------------
  subroutine read_exodus_bdy_coords( bdynum, lat_offset, lon_offset)

    ! initialize
    integer, intent(in) :: bdynum
    real,    intent(in) :: lat_offset(1), lon_offset(1)

    ! internal
    include 'netcdf.inc'
    integer stat
    integer :: utmx_offset(1), utmy_offset(1)
    character(4) :: utmzone_offset(1)

    !================================================================================
    ! Read mesh coordinates
    stat = nf_get_var_real(bdy(bdynum)%ofid,bdy(bdynum)%coordx_id,bdy(bdynum)%coordx)
    call ncderrcheck( __LINE__ ,stat )
    stat = nf_get_var_real(bdy(bdynum)%ofid,bdy(bdynum)%coordy_id,bdy(bdynum)%coordy)
    call ncderrcheck( __LINE__ ,stat )
    stat = nf_get_var_real(bdy(bdynum)%ofid,bdy(bdynum)%coordz_id,bdy(bdynum)%coordz)
    call ncderrcheck( __LINE__ ,stat )

    !================================================================================
    ! Transform these to lat and lon with the offset
    call deg2utm(1, &
         lat_offset, lon_offset, &
         utmx_offset, utmy_offset, utmzone_offset)

    bdy(bdynum)%coordx(:)  = bdy(bdynum)%coordx(:) + utmx_offset(1)
    bdy(bdynum)%coordy(:)  = bdy(bdynum)%coordy(:) + utmy_offset(1)
    bdy(bdynum)%utmzone(:) = utmzone_offset(1)

    call utm2deg(bdy(bdynum)%num_nodes, &
         nint(bdy(bdynum)%coordx), nint(bdy(bdynum)%coordy), bdy(bdynum)%utmzone, &
         bdy(bdynum)%lat, bdy(bdynum)%lon)

  end subroutine read_exodus_bdy_coords


  !--------------------------------------------------------------------------------
  !
  !> @brief Find the closest matching point in WRF
  !> @param bdynum  body index
  !> @param xlat      WRF latitude
  !> @param xlong     WRF longitude
  !> @param ids       start i index
  !> @param ide       end i index   
  !> @param jds       start j index 
  !> @param jde       end j index   
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
  subroutine relate_exodus_wrf( bdynum,xlat,xlong,ids,ide,jds,jde,ips,ipe,jps,jpe,ims,ime,jms,jme)

    ! initialize
    integer, intent(in) :: bdynum,ids,ide,jds,jde,ips,ipe,jps,jpe,ims,ime,jms,jme
    real, dimension(ids:ide-1,jds:jde-1), intent(in) :: xlat, xlong

    ! internal
    integer i,j,ipoint
    REAL exo_lat, exo_lon
    REAL dsw,dse,dnw,dne,lim,dmin

    allocate( bdy(bdynum)%exo_wrf_i(bdy(bdynum)%num_nodes))
    allocate( bdy(bdynum)%exo_wrf_j(bdy(bdynum)%num_nodes))

    ! Loop on all the nodes
    do ipoint = 1, bdy(bdynum)%num_nodes

       exo_lat = bdy(bdynum)%lat(ipoint)
       exo_lon = bdy(bdynum)%lon(ipoint)
       dmin = 999999.9

       ! loop on WRF data
       do j = jps,min(jpe,jde-2)
          do i = ips,min(ipe,ide-2)
             ! ignore special case where of point lies outside the grid of cell centers
             ! should not put EXO grid that close to a WRF boundary
             ! also note the cavalier way we ignore curvature and assume the
             ! grid cells are perfectly square and that lat and lon are Cartesian
             dsw = sqrt((exo_lat-xlat(i ,j ))*(exo_lat-xlat(i ,j )) + (exo_lon-xlong(i ,j ))*(exo_lon-xlong(i ,j )))
             !!absolute closest
             !IF ( dsw .LT. dmin ) THEN
             !alternate scheme, pick the point that is closest to the sw of the exodus point
             if ( dsw .lt. dmin .and. exo_lat .ge. xlat(i,j) .and. exo_lon .ge. xlong(i,j) ) then
                bdy(bdynum)%exo_wrf_i(ipoint) = i
                bdy(bdynum)%exo_wrf_j(ipoint) = j
                dmin = dsw
             endif
          enddo
       enddo
    enddo

  end subroutine relate_exodus_wrf


  !--------------------------------------------------------------------------------
  !
  !> @brief Write out the necessary variables to the new Exodus file
  !> @param bdynum        body index
  !
  !--------------------------------------------------------------------------------
  subroutine write_vars_exodus( bdynum )

    ! initialize
    integer, intent(in) :: bdynum

    ! internal
    include 'netcdf.inc'
    integer stat

    stat = nf_put_var_real(bdy(bdynum)%ofid,bdy(bdynum)%vals_nod_var1_id,bdy(bdynum)%vals_nod_var1)
    call ncderrcheck( __LINE__ ,stat )
    stat = nf_put_var_real(bdy(bdynum)%ofid,bdy(bdynum)%vals_nod_var2_id,bdy(bdynum)%vals_nod_var2)
    call ncderrcheck( __LINE__ ,stat )
    stat = nf_put_var_real(bdy(bdynum)%ofid,bdy(bdynum)%vals_nod_var3_id,bdy(bdynum)%vals_nod_var3)
    call ncderrcheck( __LINE__ ,stat )
    stat = nf_put_var_real(bdy(bdynum)%ofid,bdy(bdynum)%vals_nod_var4_id,bdy(bdynum)%vals_nod_var4)
    call ncderrcheck( __LINE__ ,stat )
    stat = nf_put_var_real(bdy(bdynum)%ofid,bdy(bdynum)%vals_nod_var5_id,bdy(bdynum)%vals_nod_var5)
    call ncderrcheck( __LINE__ ,stat )
    stat = nf_put_var_real(bdy(bdynum)%ofid,bdy(bdynum)%vals_nod_var6_id,bdy(bdynum)%vals_nod_var6)
    call ncderrcheck( __LINE__ ,stat )
    stat = nf_put_var_real(bdy(bdynum)%ofid,bdy(bdynum)%vals_nod_var7_id,bdy(bdynum)%vals_nod_var7)
    call ncderrcheck( __LINE__ ,stat )
    
  end subroutine write_vars_exodus
  
  
  !--------------------------------------------------------------------------------
  !
  !> @brief Calculate the rotation angle of the WRF w.r.t. true
  !> @param xlat      latitude
  !> @param dx        mesh spacing
  !> @param ids       start i index
  !> @param ide       end i index   
  !> @param jds       start j index 
  !> @param jde       end j index   
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

END MODULE module_exodus
