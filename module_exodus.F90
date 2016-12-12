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
     logical using_hfx
     integer num_nodes
     integer num_nod_var_dimid
     integer ofid, time_whole_id, coordx_id, coordy_id, coordz_id, &
          vals_nod_var1_id, vals_nod_var2_id, vals_nod_var3_id, &
          vals_nod_var4_id, vals_nod_var5_id, vals_nod_var6_id, &
          vals_nod_var7_id, vals_nod_var8_id
     integer,      dimension(:), allocatable :: exo_wrf_i, exo_wrf_j
     real,         dimension(:), allocatable :: coordx, coordy, coordz, lat, lon
     real,         dimension(:), allocatable ::vals_nod_var1, vals_nod_var2, vals_nod_var3, &
          vals_nod_var4, vals_nod_var5, vals_nod_var6, vals_nod_var7, vals_nod_var8

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
  !> @param ibdy         body index
  !> @param ofname       Exodus BC file data
  !> @param using_hfx    Whether to have a heat flux
  !
  !--------------------------------------------------------------------------------
  subroutine prep_exodus( ibdy, ofname, using_hfx)

    ! initialize
    integer,       intent(in) :: ibdy
    character*(*), intent(in) :: ofname
    logical,       intent(in) :: using_hfx

    ! internal
    include 'netcdf.inc'
    integer stat

    ! Store time information
    integer,dimension(8) :: date_values
    character(12)  :: date_str
    character(2) :: dd
    character(2) :: mm
    character(4) :: yyyy

    ! number of variables
    integer num_vars

    ! Variable ids and values
    integer num_nodes_dimid, time_step_dimid, len_name_dimid, &
         len_line_dimid, num_info_dimid
    integer info_records_id, eb_names_id,  &
         vals_nod_var1_dims(2), vals_nod_var2_dims(2), &
         vals_nod_var3_dims(2), vals_nod_var4_dims(2), &
         vals_nod_var5_dims(2), vals_nod_var6_dims(2), &
         vals_nod_var7_dims(2), vals_nod_var8_dims(2), &
         name_nod_var_id, name_nod_var_dims(2)

    character(len=255) :: info_records
    integer :: info_records_dims(2), info_records_start(2), info_records_count(2)

    character(len=255) :: eb_names
    integer :: eb_names_start(2), eb_names_count(2)

    character(len=33), dimension(:), allocatable :: name_nod_var
    integer :: name_nod_var_start(2), name_nod_var_count(2)

    ! heat flux?
    bdy(ibdy)%using_hfx = using_hfx
    if (bdy(ibdy)%using_hfx) then
       num_vars = 8
    else
       num_vars = 7
    endif

    WRITE(0,*)'opening output file ',TRIM(ofname)
    stat = nf_open(ofname, nf_write, bdy(ibdy)%ofid )
    call ncderrcheck( __FILE__, __LINE__ ,stat )

    !================================================================================
    ! Get system information
    call date_and_time(VALUES=date_values)
    write(  dd,'(i2)') date_values(3)
    write(  mm,'(i2)') date_values(2)
    write(yyyy,'(i4)') date_values(1)
    write(date_str,*),dd,'/',mm,'/',yyyy

    !================================================================================
    ! Get mesh information
    stat = nf_inq_dimid(bdy(ibdy)%ofid,"num_nodes", num_nodes_dimid)
    call ncderrcheck( __FILE__, __LINE__,stat)
    stat = nf_inq_dimlen(bdy(ibdy)%ofid, num_nodes_dimid, bdy(ibdy)%num_nodes)
    call ncderrcheck( __FILE__, __LINE__,stat)

    stat = nf_inq_dimid(bdy(ibdy)%ofid,"time_step", time_step_dimid)
    call ncderrcheck( __FILE__, __LINE__,stat)

    stat = nf_inq_dimid(bdy(ibdy)%ofid,"len_line", len_line_dimid)
    call ncderrcheck( __FILE__, __LINE__,stat)

    stat = nf_inq_dimid(bdy(ibdy)%ofid,"len_name", len_name_dimid)
    call ncderrcheck( __FILE__, __LINE__,stat)

    stat = nf_inq_varid(bdy(ibdy)%ofid,"info_records", info_records_id)
    if (stat .ne. 0) then ! if it doesn't exist, define it
       stat = nf_redef(bdy(ibdy)%ofid)
       call ncderrcheck( __FILE__, __LINE__,stat)

       stat = nf_def_dim(bdy(ibdy)%ofid, "num_info", 1 , num_info_dimid)
       call ncderrcheck( __FILE__, __LINE__,stat)

       info_records_dims(1) = len_line_dimid
       info_records_dims(2) = num_info_dimid
       stat = nf_def_var(bdy(ibdy)%ofid, "info_records", nf_char, 2, info_records_dims , info_records_id)
       call ncderrcheck( __FILE__, __LINE__,stat)

       stat = nf_enddef(bdy(ibdy)%ofid)
       call ncderrcheck( __FILE__, __LINE__,stat)
    endif
    call ncderrcheck( __FILE__, __LINE__,stat)

    stat = nf_inq_varid(bdy(ibdy)%ofid,"eb_names", eb_names_id)
    call ncderrcheck( __FILE__, __LINE__,stat)

    stat = nf_inq_varid(bdy(ibdy)%ofid,"coordx", bdy(ibdy)%coordx_id)
    call ncderrcheck( __FILE__, __LINE__,stat)

    stat = nf_inq_varid(bdy(ibdy)%ofid,"coordy", bdy(ibdy)%coordy_id)
    call ncderrcheck( __FILE__, __LINE__,stat)

    stat = nf_inq_varid(bdy(ibdy)%ofid,"coordz", bdy(ibdy)%coordz_id)
    call ncderrcheck( __FILE__, __LINE__,stat)

    !================================================================================
    ! Define new dimensions and variables

    ! put in define mode
    stat = nf_redef(bdy(ibdy)%ofid)
    call ncderrcheck( __FILE__, __LINE__,stat)

    ! Write out num_nod_var dimension
    stat = nf_def_dim(bdy(ibdy)%ofid, "num_nod_var", num_vars , bdy(ibdy)%num_nod_var_dimid)
    call ncderrcheck( __FILE__, __LINE__,stat)

    ! vals_nod_var1 variable
    vals_nod_var1_dims(1) = num_nodes_dimid
    vals_nod_var1_dims(2) = time_step_dimid
    stat = nf_def_var(bdy(ibdy)%ofid, "vals_nod_var1", nf_double, 2, vals_nod_var1_dims , bdy(ibdy)%vals_nod_var1_id)
    call ncderrcheck( __FILE__, __LINE__,stat)

    ! vals_nod_var2 variable
    vals_nod_var2_dims(1) = num_nodes_dimid
    vals_nod_var2_dims(2) = time_step_dimid
    stat = nf_def_var(bdy(ibdy)%ofid, "vals_nod_var2", nf_double, 2, vals_nod_var2_dims , bdy(ibdy)%vals_nod_var2_id)
    call ncderrcheck( __FILE__, __LINE__,stat)

    ! vals_nod_var3 variable
    vals_nod_var3_dims(1) = num_nodes_dimid
    vals_nod_var3_dims(2) = time_step_dimid
    stat = nf_def_var(bdy(ibdy)%ofid, "vals_nod_var3", nf_double, 2, vals_nod_var3_dims , bdy(ibdy)%vals_nod_var3_id)
    call ncderrcheck( __FILE__, __LINE__,stat)

    ! vals_nod_var4 variable
    vals_nod_var4_dims(1) = num_nodes_dimid
    vals_nod_var4_dims(2) = time_step_dimid
    stat = nf_def_var(bdy(ibdy)%ofid, "vals_nod_var4", nf_double, 2, vals_nod_var4_dims , bdy(ibdy)%vals_nod_var4_id)
    call ncderrcheck( __FILE__, __LINE__,stat)

    ! vals_nod_var5 variable
    vals_nod_var5_dims(1) = num_nodes_dimid
    vals_nod_var5_dims(2) = time_step_dimid
    stat = nf_def_var(bdy(ibdy)%ofid, "vals_nod_var5", nf_double, 2, vals_nod_var5_dims , bdy(ibdy)%vals_nod_var5_id)
    call ncderrcheck( __FILE__, __LINE__,stat)

    ! vals_nod_var6 variable
    vals_nod_var6_dims(1) = num_nodes_dimid
    vals_nod_var6_dims(2) = time_step_dimid
    stat = nf_def_var(bdy(ibdy)%ofid, "vals_nod_var6", nf_double, 2, vals_nod_var6_dims , bdy(ibdy)%vals_nod_var6_id)
    call ncderrcheck( __FILE__, __LINE__,stat)

    ! vals_nod_var7 variable
    vals_nod_var7_dims(1) = num_nodes_dimid
    vals_nod_var7_dims(2) = time_step_dimid
    stat = nf_def_var(bdy(ibdy)%ofid, "vals_nod_var7", nf_double, 2, vals_nod_var7_dims , bdy(ibdy)%vals_nod_var7_id)
    call ncderrcheck( __FILE__, __LINE__,stat)

    ! vals_nod_var8 variable
    if (bdy(ibdy)%using_hfx) then
       vals_nod_var8_dims(1) = num_nodes_dimid
       vals_nod_var8_dims(2) = time_step_dimid
       stat = nf_def_var(bdy(ibdy)%ofid, "vals_nod_var8", nf_double, 2, vals_nod_var8_dims , bdy(ibdy)%vals_nod_var8_id)
       call ncderrcheck( __FILE__, __LINE__,stat)
    endif
    
    ! name_nod_var variable
    name_nod_var_dims(1) = len_name_dimid
    name_nod_var_dims(2) = bdy(ibdy)%num_nod_var_dimid
    stat = nf_def_var(bdy(ibdy)%ofid, "name_nod_var", nf_char, 2, name_nod_var_dims , name_nod_var_id)
    call ncderrcheck( __FILE__, __LINE__,stat)

    ! End define mode
    stat = nf_enddef(bdy(ibdy)%ofid)
    call ncderrcheck( __FILE__, __LINE__,stat)

    !================================================================================
    ! Fill in some of the variables
   
    ! info_records variable
    info_records = "Made with WRFTONALU on"//trim(date_str)
    info_records_start(1) = 1                              ! start at beginning of variable
    info_records_start(2) = 1                              ! record number to write
    info_records_count(1) = len(trim(info_records))        ! number of chars to write
    info_records_count(2) = 1                              ! only write one record
    stat = nf_put_vara_text(bdy(ibdy)%ofid, info_records_id, info_records_start, info_records_count, info_records)
    call ncderrcheck( __FILE__, __LINE__ ,stat )

    ! eb_names
    if (ibdy .ne. interior) then
       eb_names =  "block_101"
       eb_names_start(1) = 1                              ! start at beginning of variable
       eb_names_start(2) = 1                              ! record number to write
       eb_names_count(1) = LEN(trim(eb_names))            ! number of chars to write
       eb_names_count(2) = 1                              ! only write one record
       stat = nf_put_vara_text(bdy(ibdy)%ofid, eb_names_id, eb_names_start, eb_names_count, eb_names)
       call ncderrcheck( __FILE__, __LINE__ ,stat )
    endif

    ! name_nod_var
    allocate(name_nod_var(num_vars))
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
    stat = nf_put_vara_text(bdy(ibdy)%ofid, name_nod_var_id, name_nod_var_start, name_nod_var_count, name_nod_var(1))
    call ncderrcheck( __FILE__, __LINE__ ,stat )
    name_nod_var_start(1) = 1
    name_nod_var_start(2) = 2
    name_nod_var_count(1) = LEN(trim(name_nod_var(2)))
    name_nod_var_count(2) = 1
    stat = nf_put_vara_text(bdy(ibdy)%ofid, name_nod_var_id, name_nod_var_start, name_nod_var_count, name_nod_var(2))
    call ncderrcheck( __FILE__, __LINE__ ,stat )
    name_nod_var_start(1) = 1
    name_nod_var_start(2) = 3
    name_nod_var_count(1) = LEN(trim(name_nod_var(3)))
    name_nod_var_count(2) = 1
    stat = nf_put_vara_text(bdy(ibdy)%ofid, name_nod_var_id, name_nod_var_start, name_nod_var_count, name_nod_var(3))
    call ncderrcheck( __FILE__, __LINE__ ,stat )
    name_nod_var_start(1) = 1
    name_nod_var_start(2) = 4
    name_nod_var_count(1) = LEN(trim(name_nod_var(4)))
    name_nod_var_count(2) = 1
    stat = nf_put_vara_text(bdy(ibdy)%ofid, name_nod_var_id, name_nod_var_start, name_nod_var_count, name_nod_var(4))
    call ncderrcheck( __FILE__, __LINE__ ,stat )
    name_nod_var_start(1) = 1
    name_nod_var_start(2) = 5
    name_nod_var_count(1) = LEN(trim(name_nod_var(5)))
    name_nod_var_count(2) = 1
    stat = nf_put_vara_text(bdy(ibdy)%ofid, name_nod_var_id, name_nod_var_start, name_nod_var_count, name_nod_var(5))
    call ncderrcheck( __FILE__, __LINE__ ,stat )
    name_nod_var_start(1) = 1
    name_nod_var_start(2) = 6
    name_nod_var_count(1) = LEN(trim(name_nod_var(6)))
    name_nod_var_count(2) = 1
    stat = nf_put_vara_text(bdy(ibdy)%ofid, name_nod_var_id, name_nod_var_start, name_nod_var_count, name_nod_var(6))
    call ncderrcheck( __FILE__, __LINE__ ,stat )
    name_nod_var_start(1) = 1
    name_nod_var_start(2) = 7
    name_nod_var_count(1) = LEN(trim(name_nod_var(7)))
    name_nod_var_count(2) = 1
    stat = nf_put_vara_text(bdy(ibdy)%ofid, name_nod_var_id, name_nod_var_start, name_nod_var_count, name_nod_var(7))
    call ncderrcheck( __FILE__, __LINE__ ,stat )
    if (bdy(ibdy)%using_hfx) then
       name_nod_var(8) = 'heat_flux'
       name_nod_var_start(1) = 1
       name_nod_var_start(2) = 8
       name_nod_var_count(1) = LEN(trim(name_nod_var(8)))
       name_nod_var_count(2) = 1
       stat = nf_put_vara_text(bdy(ibdy)%ofid, name_nod_var_id, name_nod_var_start, name_nod_var_count, name_nod_var(8))
    endif
    deallocate(name_nod_var)
    
    !================================================================================
    ! Allocate variables while we are at it
    allocate( bdy(ibdy)%coordx(bdy(ibdy)%num_nodes))
    allocate( bdy(ibdy)%coordy(bdy(ibdy)%num_nodes))
    allocate( bdy(ibdy)%coordz(bdy(ibdy)%num_nodes))
    allocate( bdy(ibdy)%utmzone(bdy(ibdy)%num_nodes))
    allocate( bdy(ibdy)%lat(bdy(ibdy)%num_nodes))
    allocate( bdy(ibdy)%lon(bdy(ibdy)%num_nodes))

    allocate( bdy(ibdy)%vals_nod_var1(bdy(ibdy)%num_nodes))
    allocate( bdy(ibdy)%vals_nod_var2(bdy(ibdy)%num_nodes))
    allocate( bdy(ibdy)%vals_nod_var3(bdy(ibdy)%num_nodes))
    allocate( bdy(ibdy)%vals_nod_var4(bdy(ibdy)%num_nodes))
    allocate( bdy(ibdy)%vals_nod_var5(bdy(ibdy)%num_nodes))
    allocate( bdy(ibdy)%vals_nod_var6(bdy(ibdy)%num_nodes))
    allocate( bdy(ibdy)%vals_nod_var7(bdy(ibdy)%num_nodes))
    if (bdy(ibdy)%using_hfx) then
       allocate( bdy(ibdy)%vals_nod_var8(bdy(ibdy)%num_nodes))
       bdy(ibdy)%vals_nod_var8 = 0.0
    endif
    
    bdy(ibdy)%coordx  = 0.0
    bdy(ibdy)%coordy  = 0.0
    bdy(ibdy)%coordz  = 0.0
    bdy(ibdy)%utmzone = ''
    bdy(ibdy)%lat     = 0.0
    bdy(ibdy)%lon     = 0.0

    bdy(ibdy)%vals_nod_var1 = 0.0
    bdy(ibdy)%vals_nod_var2 = 0.0
    bdy(ibdy)%vals_nod_var3 = 0.0
    bdy(ibdy)%vals_nod_var4 = 0.0
    bdy(ibdy)%vals_nod_var5 = 0.0
    bdy(ibdy)%vals_nod_var6 = 0.0
    bdy(ibdy)%vals_nod_var7 = 0.0
    
  end subroutine  prep_exodus


  !--------------------------------------------------------------------------------
  !
  !> @brief Close netCDF file
  !> @param ibdy  body index 
  !
  !--------------------------------------------------------------------------------
  subroutine close_exodus( ibdy )

    ! initialize
    integer,       intent(in) :: ibdy

    ! internal
    include 'netcdf.inc'
    integer stat

    write(*,*)'closing output file'
    
    ! Close the file
    stat = nf_close(bdy(ibdy)%ofid);
    call ncderrcheck( __FILE__, __LINE__ ,stat )

    ! Deallocate
    deallocate(bdy(ibdy)%exo_wrf_i)
    deallocate(bdy(ibdy)%exo_wrf_j)
    deallocate(bdy(ibdy)%coordx)
    deallocate(bdy(ibdy)%coordy)
    deallocate(bdy(ibdy)%coordz)
    deallocate(bdy(ibdy)%lat)
    deallocate(bdy(ibdy)%lon)
    deallocate(bdy(ibdy)%vals_nod_var1)
    deallocate(bdy(ibdy)%vals_nod_var2)
    deallocate(bdy(ibdy)%vals_nod_var3)
    deallocate(bdy(ibdy)%vals_nod_var4)
    deallocate(bdy(ibdy)%vals_nod_var5)
    deallocate(bdy(ibdy)%vals_nod_var6)
    deallocate(bdy(ibdy)%vals_nod_var7)
    if(bdy(ibdy)%using_hfx) then
       deallocate(bdy(ibdy)%vals_nod_var8)
    endif
    deallocate(bdy(ibdy)%utmzone)

  end subroutine close_exodus


  !--------------------------------------------------------------------------------
  !
  !> @brief Read in Exodus coordinates from a mesh
  !> @param ibdy          body index
  !> @param lat_offset    latitude offset
  !> @param lon_offset    longitude offset
  !
  !--------------------------------------------------------------------------------
  subroutine read_exodus_bdy_coords( ibdy, lat_offset, lon_offset)

    ! initialize
    integer, intent(in) :: ibdy
    real,    intent(in) :: lat_offset(1), lon_offset(1)

    ! internal
    include 'netcdf.inc'
    integer stat
    integer :: utmx_offset(1), utmy_offset(1)
    character(4) :: utmzone_offset(1)

    !================================================================================
    ! Read mesh coordinates
    stat = nf_get_var_real(bdy(ibdy)%ofid,bdy(ibdy)%coordx_id,bdy(ibdy)%coordx)
    call ncderrcheck( __FILE__, __LINE__ ,stat )
    stat = nf_get_var_real(bdy(ibdy)%ofid,bdy(ibdy)%coordy_id,bdy(ibdy)%coordy)
    call ncderrcheck( __FILE__, __LINE__ ,stat )
    stat = nf_get_var_real(bdy(ibdy)%ofid,bdy(ibdy)%coordz_id,bdy(ibdy)%coordz)
    call ncderrcheck( __FILE__, __LINE__ ,stat )

    !================================================================================
    ! Transform these to lat and lon with the offset
    call deg2utm(1, &
         lat_offset, lon_offset, &
         utmx_offset, utmy_offset, utmzone_offset)

    bdy(ibdy)%coordx(:)  = bdy(ibdy)%coordx(:) + utmx_offset(1)
    bdy(ibdy)%coordy(:)  = bdy(ibdy)%coordy(:) + utmy_offset(1)
    bdy(ibdy)%utmzone(:) = utmzone_offset(1)

    call utm2deg(bdy(ibdy)%num_nodes, &
         nint(bdy(ibdy)%coordx), nint(bdy(ibdy)%coordy), bdy(ibdy)%utmzone, &
         bdy(ibdy)%lat, bdy(ibdy)%lon)

  end subroutine read_exodus_bdy_coords


  !--------------------------------------------------------------------------------
  !
  !> @brief Find the closest matching point in WRF
  !> @param ibdy      body index
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
  subroutine relate_exodus_wrf( ibdy,xlat,xlong,ids,ide,jds,jde,ips,ipe,jps,jpe,ims,ime,jms,jme)

    ! initialize
    integer, intent(in) :: ibdy,ids,ide,jds,jde,ips,ipe,jps,jpe,ims,ime,jms,jme
    real, dimension(ids:ide-1,jds:jde-1), intent(in) :: xlat, xlong

    ! internal
    integer i,j,ipoint
    REAL exo_lat, exo_lon
    REAL dsw,dse,dnw,dne,lim,dmin

    allocate( bdy(ibdy)%exo_wrf_i(bdy(ibdy)%num_nodes))
    allocate( bdy(ibdy)%exo_wrf_j(bdy(ibdy)%num_nodes))

    ! Loop on all the nodes
    do ipoint = 1, bdy(ibdy)%num_nodes

       exo_lat = bdy(ibdy)%lat(ipoint)
       exo_lon = bdy(ibdy)%lon(ipoint)
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
                bdy(ibdy)%exo_wrf_i(ipoint) = i
                bdy(ibdy)%exo_wrf_j(ipoint) = j
                dmin = dsw
             endif
          enddo
       enddo
    enddo

  end subroutine relate_exodus_wrf


  !--------------------------------------------------------------------------------
  !
  !> @brief Write out the necessary variables to the new Exodus file
  !> @param ibdy          body index
  !> @param itime         Time step number
  !> @param sec           Time (seconds)

  !
  !--------------------------------------------------------------------------------
  subroutine write_vars_exodus( ibdy, itime, sec)

    ! initialize
    integer, intent(in) :: ibdy, itime, sec

    ! internal
    include 'netcdf.inc'
    integer stat, start(2), count(2)

    ! Initialize
    count(1) = 1
    count(2) = 1
    start(1) = 1
    start(2) = 1

    ! time_whole variable
    stat = nf_inq_varid(bdy(ibdy)%ofid,"time_whole", bdy(ibdy)%time_whole_id)
    call ncderrcheck( __FILE__, __LINE__,stat)
    count(1) = 1      ! number of doubles to write
    start(1) = itime  ! record to write
    stat = nf_put_vara_double(bdy(ibdy)%ofid, bdy(ibdy)%time_whole_id, start, count, dble(sec))
    call ncderrcheck( __FILE__, __LINE__ ,stat )
    
    ! All other variables
    count(1) = bdy(ibdy)%num_nodes ! number of doubles to write
    count(2) = 1                   ! only write one record
    start(1) = 1                   ! start at beginning of variable
    start(2) = itime               ! record number to write                        

    stat = nf_put_vara_real(bdy(ibdy)%ofid,bdy(ibdy)%vals_nod_var1_id, start, count,bdy(ibdy)%vals_nod_var1)
    call ncderrcheck( __FILE__, __LINE__ ,stat )
    stat = nf_put_vara_real(bdy(ibdy)%ofid,bdy(ibdy)%vals_nod_var2_id, start, count,bdy(ibdy)%vals_nod_var2)
    call ncderrcheck( __FILE__, __LINE__ ,stat )
    stat = nf_put_vara_real(bdy(ibdy)%ofid,bdy(ibdy)%vals_nod_var3_id, start, count,bdy(ibdy)%vals_nod_var3)
    call ncderrcheck( __FILE__, __LINE__ ,stat )
    stat = nf_put_vara_real(bdy(ibdy)%ofid,bdy(ibdy)%vals_nod_var4_id, start, count,bdy(ibdy)%vals_nod_var4)
    call ncderrcheck( __FILE__, __LINE__ ,stat )
    stat = nf_put_vara_real(bdy(ibdy)%ofid,bdy(ibdy)%vals_nod_var5_id, start, count,bdy(ibdy)%vals_nod_var5)
    call ncderrcheck( __FILE__, __LINE__ ,stat )
    stat = nf_put_vara_real(bdy(ibdy)%ofid,bdy(ibdy)%vals_nod_var6_id, start, count,bdy(ibdy)%vals_nod_var6)
    call ncderrcheck( __FILE__, __LINE__ ,stat )
    stat = nf_put_vara_real(bdy(ibdy)%ofid,bdy(ibdy)%vals_nod_var7_id, start, count,bdy(ibdy)%vals_nod_var7)
    call ncderrcheck( __FILE__, __LINE__ ,stat )

    if (bdy(ibdy)%using_hfx) then
       stat = nf_put_vara_real(bdy(ibdy)%ofid,bdy(ibdy)%vals_nod_var8_id, start, count,bdy(ibdy)%vals_nod_var8)
       call ncderrcheck( __FILE__, __LINE__ ,stat )
    endif

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
