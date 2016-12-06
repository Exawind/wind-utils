!------------------------------------------------------------------------------
!
! PROGRAM: wrftonalu
!
!> @author
!> Marc T. Henry de Frahan, National Renewable Energy Laboratory
!
!> @brief
!> Reads in WRF output files and outputs Nalu Exodus II boundary conditions
!
!> @date 01/12/2016 J. Michalakes and M. Churchfield
!> - Initial version from WRFTOOOF
! 
!------------------------------------------------------------------------------

PROGRAM wrftonalu
  USE module_dm
  USE module_exodus_bc
  USE module_utmdeg_converter
  IMPLICIT NONE

  !================================================================================
  !
  ! Initialize variables
  !
  !================================================================================

  ! domain specification variables
  INTEGER :: ids , ide , jds , jde , kds , kde
  INTEGER :: ims , ime , jms , jme , kms , kme ! some of this is unnecessary for serial
  INTEGER :: ips , ipe , jps , jpe , kps , kpe ! ditto
  ! variables for using netcdf
  INCLUDE 'netcdf.inc'
  INTEGER, PARAMETER :: MAXFILES = 20
  CHARACTER(LEN=255) :: flnm(MAXFILES),arg,outnameT,outnamePd,outnameU,outnameHFX,vname,comstr,dirpath
  CHARACTER(LEN=19) :: Times(100),tmpstr,secstr
  LOGICAL ic, ctrl, have_hfx, use_hfx ! whether or not to do an IC file too
  INTEGER it,ncid(MAXFILES),stat,iarg,narg,varid,strt(4),cnt(4),xtype,storeddim,dimids(4),natts
  REAL, EXTERNAL :: finterp
  INTEGER, EXTERNAL :: sec_of_day
  LOGICAL, EXTERNAL :: valid_date
  INTEGER sec, sec_start, sec_offset, nfiles
  REAL , PARAMETER :: g = 9.81 ! acceleration due to gravity (m {s}^-2)
  DOUBLE PRECISION of_lat, of_lon, of_lz ! lat and lon of desired openfoam point
  CHARACTER*32, DIMENSION(nbdys) :: bdynames

  REAL, DIMENSION(:,:,:), ALLOCATABLE :: zz & ! height in meters
       ,w & ! w at cell center
       ,ph & ! geop pert in wrf
       ,phb & ! geop base in wrf
       ,pres ! pressure in millibars
  ! half level WRF variables
  REAL, DIMENSION(:,:,:), ALLOCATABLE :: z & ! height in meters on cell ctrs
       ,p & ! pres pert in wrf
       ,pb & ! pres base in wrf
       ,t & ! temp in K
       ,u_ & ! staggered u_
       ,v_ & ! staggered v_
       ,u & ! u at cell center
       ,v ! v at cell center
  ! two-d WRF variables
  REAL, DIMENSION(:,:), ALLOCATABLE :: xlat, xlong , hfx
  ! temporaries
  REAL, DIMENSION(:,:,:), ALLOCATABLE :: zzcol
  REAL, DIMENSION(:,:,:), ALLOCATABLE :: zcol
  INTEGER kz(0:1,0:1), kzz(0:1,0:1), ibdy, ipoint
  REAL :: hfx_new, u_new, v_new, t_new, t_ground, pres_new, pd, w_new, theta, costheta, sintheta, dx
  REAL :: dx_check
  INTEGER :: ids_check , ide_check , jds_check , jde_check , kds_check , kde_check
  INTEGER :: i , j , k
  INTEGER :: ii, jj, kk

  ! borrowed from NCL for computing T from Theta
  DOUBLE PRECISION PI
  DOUBLE PRECISION P1000MB,R_D,CP, RHO0
  PARAMETER (P1000MB=100000.D0,R_D=266.9D0,CP=7.D0*R_D/2.D0)


  ! mine

  integer, dimension(6) :: test_x = (/458731,  407653,  239027,  230253,  343898,  362850/)
  integer, dimension(6) :: test_y = (/4462881, 5126290, 4163083, 3171843, 4302285, 2772478/)
  character(len=4), dimension(6) :: test_utmzone = (/'30 T', '32 T', '11 S', '28 R', '15 S', '51 R'/)
  real, dimension(6) :: test_lat, test_lon

  
  integer,dimension(8) :: date_values
  character(12)  :: date_str
  CHARACTER(2) :: dd
  CHARACTER(2) :: mm
  CHARACTER(4) :: yyyy

  ! Exodus mesh (lat,lon) offset
  double precision :: exo_lat_offset, exo_lon_offset
  
  character(len=255) :: meshname
  character(len=255) :: ofname
  integer ncmeshid, ofid
  integer num_side_sets_id, num_side_sets
  integer len_string_dimid, len_line_dimid, four_dimid, num_info_dimid &
       , num_qa_rec_dimid, len_name_dimid, num_dim_dimid, time_step_dimid &
       , num_nodes_dimid, num_elem_dimid, num_el_blk_dimid, num_el_in_blk1_dimid &
       , num_nod_per_el1_dimid, num_nod_var_dimid


  integer num_nodes
  
  integer info_records_id
  integer info_records_dims(2)
  character(len=255) :: info_records

  integer qa_records_id
  integer qa_records_dims(3)

  integer time_whole_id
  integer time_whole_dims
  double precision :: time_whole(6)

  integer node_num_map_id
  integer node_num_map_dims

  integer elem_num_map_id
  integer elem_num_map_dims

  integer eb_status_id
  integer eb_status_dims

  integer eb_prop1_id
  integer eb_prop1_dims

  integer eb_names_id
  integer eb_names_dims(2)
  character(len=255) :: eb_names
  
  integer coordx_id
  integer coordx_dims(2)

  integer coordy_id
  integer coordy_dims(2)

  integer coordz_id
  integer coordz_dims(2)

  integer coor_names_id
  integer coor_names_dims(2)
  character(len=1), dimension(3) :: coor_names = (/'x', 'y', 'z'/)
  integer :: coor_names_start(2) = (/1,1/)
  integer :: coor_names_count(2) = (/1,3/)
  
  integer connect1_id
  integer connect1_dims(2)

  integer vals_nod_var1_id
  integer vals_nod_var1_dims(2)

  integer vals_nod_var2_id
  integer vals_nod_var2_dims(2)
  
  integer vals_nod_var3_id
  integer vals_nod_var3_dims(2)

  integer vals_nod_var4_id
  integer vals_nod_var4_dims(2)

  integer vals_nod_var5_id
  integer vals_nod_var5_dims(2)

  integer vals_nod_var6_id
  integer vals_nod_var6_dims(2)

  integer vals_nod_var7_id
  integer vals_nod_var7_dims(2)

  integer name_nod_var_id
  integer name_nod_var_dims(2)
  character(len=33), dimension(7) :: name_nod_var
  integer :: name_nod_var_start(2), name_nod_var_count(2)


  real, dimension(:), allocatable :: exo_coordx, exo_coordy, exo_lat, exo_lon
  
  !================================================================================
  !
  ! Get system information
  !
  !================================================================================
  call date_and_time(VALUES=date_values)
  write(  dd,'(i2)') date_values(3)
  write(  mm,'(i2)') date_values(2)
  write(yyyy,'(i4)') date_values(1)
  write(date_str,*),dd,'/',mm,'/',yyyy

  !================================================================================
  !
  ! Parse input arguments
  !
  !================================================================================

  ic = .FALSE.
  use_hfx = .FALSE.
  it = 1
  sec_offset = 0
  sec_start = 0
  narg = iargc()

  IF ( narg .EQ. 0 ) THEN
     CALL help
     STOP 99
  ENDIF
  iarg = 1
  nfiles = 0
  DO WHILE ( .TRUE. )
     CALL getarg(iarg,arg)
     IF ( arg(1:1) .EQ. '-' ) THEN
        IF ( TRIM(arg) .EQ. '-startdate' ) THEN
           iarg = iarg + 1
           CALL getarg(iarg,arg)
           IF ( .NOT. valid_date( arg ) ) THEN
              WRITE(0,*)'Invalid data string in third argument to command: ',TRIM(arg)
              STOP 99
           ENDIF
           sec_start = sec_of_day(arg)
        ELSE IF ( TRIM(arg) .EQ. '-offset' ) THEN
           iarg = iarg + 1
           CALL getarg(iarg,arg)
           READ(arg,*)sec_offset
        ELSE IF ( TRIM(arg) .EQ. '-ic' ) THEN
           ic = .TRUE.
        ELSE IF ( TRIM(arg) .EQ. '-qwall' ) THEN
           use_hfx = .TRUE.
        ENDIF
     ELSE
        nfiles = nfiles + 1
        IF ( nfiles .GT. MAXFILES ) THEN
           write(0,*)'Too many input files'
           STOP
        ENDIF
        flnm(nfiles) = arg
     ENDIF
     iarg = iarg + 1
     IF ( iarg .GT. narg ) exit
  ENDDO

  bdynames(BDY_XS) = "west"
  bdynames(BDY_XE) = "east"
  bdynames(BDY_YS) = "south"
  bdynames(BDY_YE) = "north"
  bdynames(BDY_ZS) = "lower"
  bdynames(BDY_ZE) = "upper"
  bdynames(INTERIOR) = "interior"


  !================================================================================
  !
  ! Read the WRF data files using netCDF Fortran functions
  !
  !================================================================================

  ! Open files, get the mesh spacing and indices, and perform sanity checks
  DO i = 1, nfiles
     WRITE(0,*)'opening : flnm(i) ',i,TRIM(flnm(i))
     stat = NF_OPEN(flnm(i), NF_NOWRITE, ncid(i))
     CALL ncderrcheck( __LINE__ ,stat )
     stat=NF_GET_ATT_REAL(ncid(i),NF_GLOBAL,'DX',dx) ;
     IF( i .EQ. 1 ) dx_check = dx
     CALL ncderrcheck( __LINE__,stat )
     stat = NF_GET_ATT_INT (ncid(i),NF_GLOBAL,'WEST-EAST_PATCH_END_STAG',ide) ; ids = 1 ;
     IF( i .EQ. 1 ) ide_check = ide
     CALL ncderrcheck( __LINE__,stat )
     stat = NF_GET_ATT_INT (ncid(i),NF_GLOBAL,'SOUTH-NORTH_PATCH_END_STAG',jde) ; jds = 1 ;
     IF( i .EQ. 1 ) jde_check = jde
     CALL ncderrcheck( __LINE__,stat )
     stat = NF_GET_ATT_INT (ncid(i),NF_GLOBAL,'BOTTOM-TOP_PATCH_END_STAG',kde) ; kds = 1 ;
     IF( i .EQ. 1 ) kde_check = kde
     CALL ncderrcheck( __LINE__,stat )
     IF ( i .GT. 1 ) THEN
        stat = 0
        IF ( dx .NE. dx_check ) THEN
           stat = 1 ; write(0,*)'DX ',TRIM(flnm(i)),' does not match ',TRIM(flnm(i))
        ENDIF
        IF ( ide .NE. ide_check ) THEN
           stat = 1 ; write(0,*)'WEST-EAST_PATCH_END_STAG in ',TRIM(flnm(i)),' does not match ',TRIM(flnm(i))
        ENDIF
        IF ( jde .NE. jde_check ) THEN
           stat = 1 ; write(0,*)'SOUTH-NORTH_PATCH_END_STAG in ',TRIM(flnm(i)),' does not match ',TRIM(flnm(i))
        ENDIF
        IF ( kde .NE. kde_check ) THEN
           stat = 1 ; write(0,*)'BOTTOM-TOP_PATCH_END_STAG in ',TRIM(flnm(i)),' does not match ',TRIM(flnm(i))
        ENDIF
        IF ( stat .NE. 0 ) STOP
     ENDIF
  ENDDO

  strt = 1
  cnt(1) = 19

  ! Get info from the Times variable
  stat = NF_INQ_VARID(ncid,'Times',varid) ! get ID of variable Times
  CALL ncderrcheck( __LINE__,stat)
  stat = NF_INQ_VAR(ncid,varid,vname,xtype,storeddim,dimids,natts) ! get all information about Times
  write(*,*)'varid',varid
  write(*,*)'vname',vname
  write(*,*)'xtype',xtype
  write(*,*)'storeddim',storeddim
  write(*,*)'dimids',dimids
  write(*,*)'natts',natts
  CALL ncderrcheck( __LINE__,stat)
  stat = NF_INQ_DIMLEN(ncid,dimids(1),cnt(1))
  CALL ncderrcheck( __LINE__,stat)
  stat = NF_INQ_DIMLEN(ncid,dimids(2),cnt(2))
  CALL ncderrcheck( __LINE__,stat)
  stat = NF_GET_VARA_TEXT(ncid,varid,strt,cnt,Times) ! read in the Times data in the Times text
  CALL ncderrcheck( __LINE__,stat )
  DO WHILE (.TRUE.) ! just replace ':' char by '_' in the Times variable
     tmpstr = Times(it)
     i = INDEX(Times(it),':')
     IF ( i .EQ. 0 ) EXIT
     Times(it)(i:i) = '_'
  ENDDO
  write(*,*)'Times 1',trim(Times(1))
  write(*,*)'Times 2',trim(Times(2))
  write(*,*)'Times 3',trim(Times(3))
  write(*,*)'Times 4',trim(Times(4))
  write(*,*)'Times 5',trim(Times(5))
  write(*,*)'Times 6',trim(Times(6))
  write(*,*)'Times 7',trim(Times(7))
  write(*,*)'Times 8',trim(Times(8))
  write(*,*)'Times 9',trim(Times(9))
  write(*,*)'Times 10',trim(Times(10))
  write(*,*)'Times 11',trim(Times(11))
  write(*,*)'Times 12',trim(Times(12))
  write(*,*)'Times 13',trim(Times(13))
  write(*,*)'cnt',cnt(2)
  
  
  ! Allocate a lot of variables
  ips = ids ; ipe = ide
  jps = jds ; jpe = jde
  kps = kds ; kpe = kde
  ims = ids ; ime = ide
  jms = jds ; jme = jde
  kms = kds ; kme = kde

  ALLOCATE( xlat(ips:ipe-1,jps:jpe-1))
  ALLOCATE( xlong(ips:ipe-1,jps:jpe-1))
  ctrl = .TRUE. ! true value going in says this field is required
  CALL getvar_real(ctrl,ncid,nfiles,'XLAT' ,xlat ,it,2,ips,ipe-1,jps,jpe-1,1,1)
  CALL getvar_real(ctrl,ncid,nfiles,'XLONG',xlong,it,2,ips,ipe-1,jps,jpe-1,1,1)

  theta = rotation_angle ( xlat,dx,ids,ide,jds,jde,ips,ipe,jps,jpe,ims,ime,jms,jme )
  ! Computed theta is counterclockwise rotation in radians of the vector from X axis, so negate and
  ! convert to degrees for reporting rotation with respect to compass points
  write(*,'("WRF grid is clockwise rotated approx.",f9.5," deg. from true lat/lon. Compensating.")'),-theta*57.2957795
  costheta = cos(theta)
  sintheta = sin(theta)

  ALLOCATE( p(ips:ipe-1,jps:jpe-1,kps:kpe-1))
  ALLOCATE( pb(ips:ipe-1,jps:jpe-1,kps:kpe-1))
  ALLOCATE( ph(ips:ipe-1,jps:jpe-1,kps:kpe ))
  ALLOCATE( phb(ips:ipe-1,jps:jpe-1,kps:kpe ))
  ALLOCATE( zz(ips:ipe-1,jps:jpe-1,kps:kpe ))
  ALLOCATE( w(ips:ipe-1,jps:jpe-1,kps:kpe ))
  ALLOCATE( pres(ips:ipe-1,jps:jpe-1,kps:kpe-1))
  ALLOCATE( z(ips:ipe-1,jps:jpe-1,kps:kpe-1))
  ALLOCATE( t(ips:ipe-1,jps:jpe-1,kps:kpe-1))
  ALLOCATE( u(ips:ipe-1,jps:jpe-1,kps:kpe-1))
  ALLOCATE( v(ips:ipe-1,jps:jpe-1,kps:kpe-1))
  ALLOCATE( u_(ips:ipe ,jps:jpe-1,kps:kpe-1))
  ALLOCATE( v_(ips:ipe-1,jps:jpe ,kps:kpe-1))
  ALLOCATE(zzcol(0:1,0:1,kps:kpe ))
  ALLOCATE( zcol(0:1,0:1,kps:kpe-1))
  ALLOCATE( hfx(ips:ipe-1,jps:jpe-1))

  p = 0.
  pb = 0.
  ph = 0.
  phb = 0.
  zz = 0.
  w = 0.
  pres = 0.
  z = 0.
  t = 0.
  u = 0.
  v = 0.
  u_ = 0.
  v_ = 0.
  zzcol = 0.
  zcol = 0.
  hfx = 0.

  CALL getvar_real(ctrl,ncid,nfiles,'PH' ,ph ,it,3,ips,ipe-1,jps,jpe-1,kps,kpe )
  CALL getvar_real(ctrl,ncid,nfiles,'PHB',phb,it,3,ips,ipe-1,jps,jpe-1,kps,kpe )
  CALL getvar_real(ctrl,ncid,nfiles,'W' ,w ,it,3,ips,ipe-1,jps,jpe-1,kps,kpe )
  CALL getvar_real(ctrl,ncid,nfiles,'T' ,t ,it,3,ips,ipe-1,jps,jpe-1,kps,kpe-1)
  CALL getvar_real(ctrl,ncid,nfiles,'U' ,u_ ,it,3,ips,ipe ,jps,jpe-1,kps,kpe-1)
  CALL getvar_real(ctrl,ncid,nfiles,'V' ,v_ ,it,3,ips,ipe-1,jps,jpe ,kps,kpe-1)
  CALL getvar_real(ctrl,ncid,nfiles,'P' ,p ,it,3,ips,ipe-1,jps,jpe-1,kps,kpe-1)
  CALL getvar_real(ctrl,ncid,nfiles,'PB' ,pb ,it,3,ips,ipe-1,jps,jpe-1,kps,kpe-1)

  have_hfx = .FALSE. ! false value going in says this field is required
  CALL getvar_real(have_hfx,ncid,nfiles,'HFX',hfx,it,3,ips,ipe-1,jps,jpe-1,1,1)

  zz = (ph + phb )/g
  z = (zz(:,:,kps:kpe-1) + zz(:,:,kps+1:kpe))*0.5
  u = (u_(ips:ipe-1,:,:)+u_(ips+1:ipe,:,:))*0.5
  v = (v_(:,jps:jpe-1,:)+v_(:,jps+1:jpe,:))*0.5

  pres = p + pb

  t = t+300.

  DEALLOCATE(ph)
  DEALLOCATE(phb)
  DEALLOCATE(u_)
  DEALLOCATE(v_)

  !================================================================================
  !
  ! Read the Exodus mesh and interpolate WRF data there
  !
  !================================================================================
  write(*,*)"Reading Exodus mesh now"

  ! Copy the mesh file to an output file
  meshname = "front.g"
  ofname = "front.nc"
  comstr = "cp " // TRIM(meshname)//" "//TRIM(ofname)
  CALL system (TRIM(comstr))

  ! Open output file
  WRITE(0,*)'opening output file ',TRIM(ofname)
  stat = NF_OPEN(ofname, NF_WRITE, ofid)
  CALL ncderrcheck( __LINE__ ,stat )

  !================================================================================
  ! Get mesh information
  stat = NF_INQ_DIMID(ofid,"num_nodes", num_nodes_dimid)
  CALL ncderrcheck( __LINE__,stat)
  stat = NF_INQ_DIMLEN(ofid, num_nodes_dimid, num_nodes)
  CALL ncderrcheck( __LINE__,stat)

  stat = NF_INQ_DIMID(ofid,"time_step", time_step_dimid)
  CALL ncderrcheck( __LINE__,stat)

  stat = NF_INQ_DIMID(ofid,"len_name", len_name_dimid)
  CALL ncderrcheck( __LINE__,stat)

  stat = NF_INQ_VARID(ofid,"info_records", info_records_id)
  CALL ncderrcheck( __LINE__,stat)

  stat = NF_INQ_VARID(ofid,"eb_names", eb_names_id)
  CALL ncderrcheck( __LINE__,stat)

  stat = NF_INQ_VARID(ofid,"time_whole", time_whole_id)
  CALL ncderrcheck( __LINE__,stat)

  stat = NF_INQ_VARID(ofid,"coordx", coordx_id)
  CALL ncderrcheck( __LINE__,stat)

  stat = NF_INQ_VARID(ofid,"coordy", coordy_id)
  CALL ncderrcheck( __LINE__,stat)

  !================================================================================
  ! Define new dimensions and variables

  ! put in define mode
  stat = NF_REDEF(ofid)
  CALL ncderrcheck( __LINE__,stat)
  
  ! Write out num_nod_var dimension
  stat = nf_def_dim(ofid, "num_nod_var", 7 , num_nod_var_dimid) ! CHANGE
  CALL ncderrcheck( __LINE__,stat)

  ! vals_nod_var1 variable
  vals_nod_var1_dims(1) = num_nodes_dimid
  vals_nod_var1_dims(2) = time_step_dimid
  stat = nf_def_var(ofid, "vals_nod_var1", nf_double, 2, vals_nod_var1_dims , vals_nod_var1_id)
  CALL ncderrcheck( __LINE__,stat)

  ! vals_nod_var2 variable
  vals_nod_var2_dims(1) = num_nodes_dimid
  vals_nod_var2_dims(2) = time_step_dimid
  stat = nf_def_var(ofid, "vals_nod_var2", nf_double, 2, vals_nod_var2_dims , vals_nod_var2_id)
  CALL ncderrcheck( __LINE__,stat)

  ! vals_nod_var3 variable
  vals_nod_var3_dims(1) = num_nodes_dimid
  vals_nod_var3_dims(2) = time_step_dimid
  stat = nf_def_var(ofid, "vals_nod_var3", nf_double, 2, vals_nod_var3_dims , vals_nod_var3_id)
  CALL ncderrcheck( __LINE__,stat)

  ! vals_nod_var4 variable
  vals_nod_var4_dims(1) = num_nodes_dimid
  vals_nod_var4_dims(2) = time_step_dimid
  stat = nf_def_var(ofid, "vals_nod_var4", nf_double, 2, vals_nod_var4_dims , vals_nod_var4_id)
  CALL ncderrcheck( __LINE__,stat)

  ! vals_nod_var5 variable
  vals_nod_var5_dims(1) = num_nodes_dimid
  vals_nod_var5_dims(2) = time_step_dimid
  stat = nf_def_var(ofid, "vals_nod_var5", nf_double, 2, vals_nod_var5_dims , vals_nod_var5_id)
  CALL ncderrcheck( __LINE__,stat)

  ! vals_nod_var6 variable
  vals_nod_var6_dims(1) = num_nodes_dimid
  vals_nod_var6_dims(2) = time_step_dimid
  stat = nf_def_var(ofid, "vals_nod_var6", nf_double, 2, vals_nod_var6_dims , vals_nod_var6_id)
  CALL ncderrcheck( __LINE__,stat)

  ! vals_nod_var7 variable
  vals_nod_var7_dims(1) = num_nodes_dimid
  vals_nod_var7_dims(2) = time_step_dimid
  stat = nf_def_var(ofid, "vals_nod_var7", nf_double, 2, vals_nod_var7_dims , vals_nod_var7_id)
  CALL ncderrcheck( __LINE__,stat)

  ! name_nod_var variable
  name_nod_var_dims(1) = len_name_dimid
  name_nod_var_dims(2) = num_nod_var_dimid
  stat = nf_def_var(ofid, "name_nod_var", nf_char, 2, name_nod_var_dims , name_nod_var_id)
  CALL ncderrcheck( __LINE__,stat)

  !================================================================================
  ! Fill in some of the variables

  ! put in data mode
  stat = NF_ENDDEF(ofid)
  CALL ncderrcheck( __LINE__,stat)

  ! info_records variable
  info_records = "Made with WRFTONALU on"//date_str
  stat = nf_put_var_text(ofid, info_records_id, trim(info_records))
  CALL ncderrcheck( __LINE__ ,stat )

  ! eb_names
  eb_names =  "block_101"
  stat = nf_put_var_text(ofid, eb_names_id, trim(eb_names))
  CALL ncderrcheck( __LINE__ ,stat )

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
  stat = nf_put_vara_text(ofid, name_nod_var_id, name_nod_var_start, name_nod_var_count, name_nod_var(1))
  CALL ncderrcheck( __LINE__ ,stat )
  name_nod_var_start(1) = 1
  name_nod_var_start(2) = 2
  name_nod_var_count(1) = LEN(trim(name_nod_var(2)))
  name_nod_var_count(2) = 1
  stat = nf_put_vara_text(ofid, name_nod_var_id, name_nod_var_start, name_nod_var_count, name_nod_var(2))
  CALL ncderrcheck( __LINE__ ,stat )
  name_nod_var_start(1) = 1
  name_nod_var_start(2) = 3
  name_nod_var_count(1) = LEN(trim(name_nod_var(3)))
  name_nod_var_count(2) = 1
  stat = nf_put_vara_text(ofid, name_nod_var_id, name_nod_var_start, name_nod_var_count, name_nod_var(3))
  CALL ncderrcheck( __LINE__ ,stat )
  name_nod_var_start(1) = 1
  name_nod_var_start(2) = 4
  name_nod_var_count(1) = LEN(trim(name_nod_var(4)))
  name_nod_var_count(2) = 1
  stat = nf_put_vara_text(ofid, name_nod_var_id, name_nod_var_start, name_nod_var_count, name_nod_var(4))
  CALL ncderrcheck( __LINE__ ,stat )
  name_nod_var_start(1) = 1
  name_nod_var_start(2) = 5
  name_nod_var_count(1) = LEN(trim(name_nod_var(5)))
  name_nod_var_count(2) = 1
  stat = nf_put_vara_text(ofid, name_nod_var_id, name_nod_var_start, name_nod_var_count, name_nod_var(5))
  CALL ncderrcheck( __LINE__ ,stat )
  name_nod_var_start(1) = 1
  name_nod_var_start(2) = 6
  name_nod_var_count(1) = LEN(trim(name_nod_var(6)))
  name_nod_var_count(2) = 1
  stat = nf_put_vara_text(ofid, name_nod_var_id, name_nod_var_start, name_nod_var_count, name_nod_var(6))
  CALL ncderrcheck( __LINE__ ,stat )
  name_nod_var_start(1) = 1
  name_nod_var_start(2) = 7
  name_nod_var_count(1) = LEN(trim(name_nod_var(7)))
  name_nod_var_count(2) = 1
  stat = nf_put_vara_text(ofid, name_nod_var_id, name_nod_var_start, name_nod_var_count, name_nod_var(7))
  CALL ncderrcheck( __LINE__ ,stat )

  !================================================================================
  ! Read the mesh coordinates
  stat = NF_INQ_VARID(ofid,"coordx",coordx_id)
  CALL ncderrcheck( __LINE__ ,stat )

  ! Define an offset (lat,long) for the mesh 
  exo_lat_offset = 33
  exo_lon_offset = -101

  ! Check to make sure it is within the WRF data set
  if ( (exo_lat_offset .le. minval(xlat)) .or. &
       (exo_lat_offset .ge. maxval(xlat)) .or. &
       (exo_lon_offset .le. minval(xlong)) .or. &
       (exo_lon_offset .ge. maxval(xlong)) ) then

     write(0,*)"Offset (lat,lon) are not contained in the WRF data set"
     write(0,*)"Offset (lat,lon)=", exo_lat_offset, exo_lon_offset
     write(0,*)"WRF data bounds min (lat,lon)=", minval(xlat), minval(xlong)
     write(0,*)"                max (lat,lon)=", maxval(xlat), maxval(xlong)
     stop 99

  endif
  
  
  ! Read in mesh coordx and coordy
  allocate( exo_coordx(num_nodes))
  allocate( exo_coordy(num_nodes))
  allocate( exo_lat(num_nodes))
  allocate( exo_lon(num_nodes))

  stat = nf_get_var_real(ofid,coordx_id,exo_coordx)
  CALL ncderrcheck( __LINE__ ,stat )
  stat = nf_get_var_real(ofid,coordy_id,exo_coordy)
  CALL ncderrcheck( __LINE__ ,stat )


  
  ! Transform these to xlat and xlong with the offset
  call utm2deg(6,test_x,test_y,test_utmzone,test_lat,test_lon)
  write(*,*)test_lat
  write(*,*)test_lon

  

  
  ! 3. for each (xlat,xlong) pair in the mesh, figure out the closest point in the WRF data set and save the (i,j) index of the WRF data set

  
  

  
  ! Close the file
  stat = NF_CLOSE(ofid);
  CALL ncderrcheck( __LINE__ ,stat )

  
  ! !================================================================================
  ! !
  ! ! Create the ouput Nalu files for each BC
  ! !
  ! !================================================================================
  ! stat = NF_CREATE("front.nc", NF_CLOBBER, ofid);
  ! CALL ncderrcheck( __LINE__,stat)

  ! ! Write out some standard dimensions
  ! stat = nf_def_dim(ofid, "len_string", 33 , len_string_dimid)
  ! CALL ncderrcheck( __LINE__,stat)
  ! stat = nf_def_dim(ofid, "len_line", 81 , len_line_dimid)
  ! CALL ncderrcheck( __LINE__,stat)
  ! stat = nf_def_dim(ofid, "four", 4 , four_dimid)
  ! CALL ncderrcheck( __LINE__,stat)
  ! stat = nf_def_dim(ofid, "num_info", 3 , num_info_dimid)
  ! CALL ncderrcheck( __LINE__,stat)
  ! stat = nf_def_dim(ofid, "num_qa_rec", 2 , num_qa_rec_dimid)
  ! CALL ncderrcheck( __LINE__,stat)
  ! stat = nf_def_dim(ofid, "len_name", 33 , len_name_dimid)
  ! CALL ncderrcheck( __LINE__,stat)
  ! stat = nf_def_dim(ofid, "num_dim", 3 , num_dim_dimid)
  ! CALL ncderrcheck( __LINE__,stat)
  ! stat = nf_def_dim(ofid, "time_step", NF_UNLIMITED , time_step_dimid) ! CHANGE
  ! CALL ncderrcheck( __LINE__,stat)
  ! stat = nf_def_dim(ofid, "num_nodes", 1681 , num_nodes_dimid) ! CHANGE
  ! CALL ncderrcheck( __LINE__,stat)
  ! stat = nf_def_dim(ofid, "num_elem", 1600 , num_elem_dimid) ! CHANGE
  ! CALL ncderrcheck( __LINE__,stat)
  ! stat = nf_def_dim(ofid, "num_el_blk", 1 , num_el_blk_dimid) ! CHANGE
  ! CALL ncderrcheck( __LINE__,stat)
  ! stat = nf_def_dim(ofid, "num_el_in_blk1", 1600 , num_el_in_blk1_dimid) ! CHANGE
  ! CALL ncderrcheck( __LINE__,stat)
  ! stat = nf_def_dim(ofid, "num_nod_per_el1", 4 , num_nod_per_el1_dimid) ! CHANGE
  ! CALL ncderrcheck( __LINE__,stat)
  ! stat = nf_def_dim(ofid, "num_nod_var", 7 , num_nod_var_dimid) ! CHANGE
  ! CALL ncderrcheck( __LINE__,stat)
  
  ! ! write some global attributes
  ! stat = nf_put_att_text(ofid, nf_global, "api_version", 5, "6.36f") !CHANGE
  ! CALL ncderrcheck( __LINE__,stat)
  ! stat = nf_put_att_text(ofid, nf_global, "version", 5, "6.36f") !CHANGE
  ! CALL ncderrcheck( __LINE__,stat)
  ! stat = nf_put_att_int(ofid, nf_global, "floating_point_word_size", nf_int, 1, 8)
  ! CALL ncderrcheck( __LINE__,stat)
  ! stat = nf_put_att_int(ofid, nf_global, "file_size", nf_int, 1, 1)
  ! CALL ncderrcheck( __LINE__,stat)
  ! stat = nf_put_att_int(ofid, nf_global, "maximum_name_length", nf_int, 1, 18)
  ! CALL ncderrcheck( __LINE__,stat)
  ! stat = nf_put_att_int(ofid, nf_global, "int64_status", nf_int, 1, 0)
  ! CALL ncderrcheck( __LINE__,stat)
  ! stat = nf_put_att_text(ofid, nf_global, "title", 1, "")
  ! CALL ncderrcheck( __LINE__,stat)
  
  ! ! info_records variable
  ! info_records_dims(1) = len_line_dimid
  ! info_records_dims(2) = num_info_dimid
  ! stat = nf_def_var(ofid, "info_records", nf_char, 2, info_records_dims , info_records_id)
  ! CALL ncderrcheck( __LINE__,stat)
  
  ! ! qa_records variable
  ! qa_records_dims(1) = len_string_dimid
  ! qa_records_dims(2) = four_dimid
  ! qa_records_dims(3) = num_qa_rec_dimid
  ! stat = nf_def_var(ofid, "qa_records", nf_char, 3, qa_records_dims , qa_records_id)
  ! CALL ncderrcheck( __LINE__,stat)
  
  ! ! time_whole variable
  ! time_whole_dims = time_step_dimid
  ! stat = nf_def_var(ofid, "time_whole", nf_double, 1, time_whole_dims , time_whole_id)
  ! CALL ncderrcheck( __LINE__,stat)
  
  ! ! node_num_map variable
  ! node_num_map_dims = num_nodes_dimid
  ! stat = nf_def_var(ofid, "node_num_map", nf_int, 1, node_num_map_dims , node_num_map_id)
  ! CALL ncderrcheck( __LINE__,stat)
  
  ! ! elem_num_map variable
  ! elem_num_map_dims = num_elem_dimid
  ! stat = nf_def_var(ofid, "elem_num_map", nf_int, 1, elem_num_map_dims , elem_num_map_id)
  ! CALL ncderrcheck( __LINE__,stat)
  
  ! ! eb_status variable
  ! eb_status_dims = num_el_blk_dimid
  ! stat = nf_def_var(ofid, "eb_status", nf_int, 1, eb_status_dims , eb_status_id)
  ! CALL ncderrcheck( __LINE__,stat)
  
  ! ! eb_prop1 variable
  ! eb_prop1_dims = num_el_blk_dimid
  ! stat = nf_def_var(ofid, "eb_prop1", nf_int, 1, eb_prop1_dims , eb_prop1_id)
  ! CALL ncderrcheck( __LINE__,stat)
  ! stat = nf_put_att_text(ofid, eb_prop1_id, "name", 2, "ID")
  ! CALL ncderrcheck( __LINE__,stat)
  
  ! ! eb_names variable
  ! eb_names_dims(1) = len_name_dimid
  ! eb_names_dims(2) = num_el_blk_dimid
  ! stat = nf_def_var(ofid, "eb_names", nf_char, 2, eb_names_dims , eb_names_id)
  ! CALL ncderrcheck( __LINE__,stat)
  
  ! ! coordx variable
  ! coordx_dims = num_nodes_dimid
  ! stat = nf_def_var(ofid, "coordx", nf_double, 1, coordx_dims , coordx_id)
  ! CALL ncderrcheck( __LINE__,stat)
  
  ! ! coordy variable
  ! coordy_dims = num_nodes_dimid
  ! stat = nf_def_var(ofid, "coordy", nf_double, 1, coordy_dims , coordy_id)
  ! CALL ncderrcheck( __LINE__,stat)
  
  ! ! coordz variable
  ! coordz_dims = num_nodes_dimid
  ! stat = nf_def_var(ofid, "coordz", nf_double, 1, coordz_dims , coordz_id)
  ! CALL ncderrcheck( __LINE__,stat)
  
  ! ! coor_names variable
  ! coor_names_dims(1) = len_name_dimid
  ! coor_names_dims(2) = num_dim_dimid
  ! stat = nf_def_var(ofid, "coor_names", nf_char, 2, coor_names_dims , coor_names_id)
  ! CALL ncderrcheck( __LINE__,stat)
  
  ! ! connect1 variable
  ! connect1_dims(1) = num_nod_per_el1_dimid
  ! connect1_dims(2) = num_el_in_blk1_dimid
  ! stat = nf_def_var(ofid, "connect1", nf_int, 2, connect1_dims , connect1_id)
  ! CALL ncderrcheck( __LINE__,stat)
  ! stat = nf_put_att_text(ofid, connect1_id, "elem_type", 6, "SHELL4")
  ! CALL ncderrcheck( __LINE__,stat)
  
  ! ! vals_nod_var1 variable
  ! vals_nod_var1_dims(1) = num_nodes_dimid
  ! vals_nod_var1_dims(2) = time_step_dimid
  ! stat = nf_def_var(ofid, "vals_nod_var1", nf_double, 2, vals_nod_var1_dims , vals_nod_var1_id)
  ! CALL ncderrcheck( __LINE__,stat)
  
  ! ! vals_nod_var2 variable
  ! vals_nod_var2_dims(1) = num_nodes_dimid
  ! vals_nod_var2_dims(2) = time_step_dimid
  ! stat = nf_def_var(ofid, "vals_nod_var2", nf_double, 2, vals_nod_var2_dims , vals_nod_var2_id)
  ! CALL ncderrcheck( __LINE__,stat)
  
  ! ! vals_nod_var3 variable
  ! vals_nod_var3_dims(1) = num_nodes_dimid
  ! vals_nod_var3_dims(2) = time_step_dimid
  ! stat = nf_def_var(ofid, "vals_nod_var3", nf_double, 2, vals_nod_var3_dims , vals_nod_var3_id)
  ! CALL ncderrcheck( __LINE__,stat)
  
  ! ! vals_nod_var4 variable
  ! vals_nod_var4_dims(1) = num_nodes_dimid
  ! vals_nod_var4_dims(2) = time_step_dimid
  ! stat = nf_def_var(ofid, "vals_nod_var4", nf_double, 2, vals_nod_var4_dims , vals_nod_var4_id)
  ! CALL ncderrcheck( __LINE__,stat)
  
  ! ! vals_nod_var5 variable
  ! vals_nod_var5_dims(1) = num_nodes_dimid
  ! vals_nod_var5_dims(2) = time_step_dimid
  ! stat = nf_def_var(ofid, "vals_nod_var5", nf_double, 2, vals_nod_var5_dims , vals_nod_var5_id)
  ! CALL ncderrcheck( __LINE__,stat)
  
  ! ! vals_nod_var6 variable
  ! vals_nod_var6_dims(1) = num_nodes_dimid
  ! vals_nod_var6_dims(2) = time_step_dimid
  ! stat = nf_def_var(ofid, "vals_nod_var6", nf_double, 2, vals_nod_var6_dims , vals_nod_var6_id)
  ! CALL ncderrcheck( __LINE__,stat)
  
  ! ! vals_nod_var7 variable
  ! vals_nod_var7_dims(1) = num_nodes_dimid
  ! vals_nod_var7_dims(2) = time_step_dimid
  ! stat = nf_def_var(ofid, "vals_nod_var7", nf_double, 2, vals_nod_var7_dims , vals_nod_var7_id)
  ! CALL ncderrcheck( __LINE__,stat)
  
  ! ! name_nod_var variable
  ! name_nod_var_dims(1) = len_name_dimid
  ! name_nod_var_dims(2) = num_nod_var_dimid
  ! stat = nf_def_var(ofid, "name_nod_var", nf_char, 2, name_nod_var_dims , name_nod_var_id)
  ! CALL ncderrcheck( __LINE__,stat)
  
  ! ! Close the file
  ! stat = NF_CLOSE(ofid);
  ! CALL ncderrcheck( __LINE__ ,stat )


  ! ! Open file in data mode to write out the variables
  ! stat = NF_OPEN("front.nc", NF_WRITE, ofid);
  ! CALL ncderrcheck( __LINE__,stat)

  ! ! info_records variable
  ! info_records(1) = "Made with WRFTONALU on"//date_str
  ! info_records(2) = ""
  ! info_records(3) = ""
  ! info_records_start(1) = 1                              ! start at beginning of variable
  ! info_records_start(2) = 1                              ! record number to write
  ! info_records_count(1) = LEN(trim(info_records(1)))     ! number of chars to write
  ! info_records_count(2) = 1                              ! only write one record
  ! stat = nf_put_vara_text(ofid, info_records_id, info_records_start, info_records_count, info_records(1))
  ! CALL ncderrcheck( __LINE__ ,stat )
  ! info_records_start(1) = 1
  ! info_records_start(2) = 2
  ! info_records_count(1) = LEN(trim(info_records(2)))
  ! info_records_count(2) = 1
  ! stat = nf_put_vara_text(ofid, info_records_id, info_records_start, info_records_count, info_records(2))
  ! CALL ncderrcheck( __LINE__ ,stat )
  ! info_records_start(1) = 1
  ! info_records_start(2) = 3
  ! info_records_count(1) = LEN(trim(info_records(3)))
  ! info_records_count(2) = 1
  ! stat = nf_put_vara_text(ofid, info_records_id, info_records_start, info_records_count, info_records(3))
  ! CALL ncderrcheck( __LINE__ ,stat )

  ! ! eb_status 
  ! stat = nf_put_var_int(ofid, eb_status_id, 1) ! CHANGE?
  ! CALL ncderrcheck( __LINE__ ,stat )

  ! ! eb_prop1
  ! stat = nf_put_var_int(ofid, eb_prop1_id, 101) ! CHANGE?
  ! CALL ncderrcheck( __LINE__ ,stat )

  ! ! eb_names
  ! eb_names =  "block_101"
  ! stat = nf_put_var_text(ofid, eb_names_id, trim(eb_names))
  ! CALL ncderrcheck( __LINE__ ,stat )

  ! ! coor_names
  ! stat = nf_put_vara_text(ofid, coor_names_id, coor_names_start, coor_names_count, coor_names)
  ! CALL ncderrcheck( __LINE__ ,stat )

  ! ! name_nod_var
  ! name_nod_var(1) = 'cont_velocity_bc_x'
  ! name_nod_var(2) = 'cont_velocity_bc_y'
  ! name_nod_var(3) = 'cont_velocity_bc_z'
  ! name_nod_var(4) = 'temperature_bc'
  ! name_nod_var(5) = 'velocity_bc_x'
  ! name_nod_var(6) = 'velocity_bc_y'
  ! name_nod_var(7) = 'velocity_bc_z'  
  ! name_nod_var_start(1) = 1                              ! start at beginning of variable
  ! name_nod_var_start(2) = 1                              ! record number to write
  ! name_nod_var_count(1) = LEN(trim(name_nod_var(1)))     ! number of chars to write
  ! name_nod_var_count(2) = 1                              ! only write one record
  ! stat = nf_put_vara_text(ofid, name_nod_var_id, name_nod_var_start, name_nod_var_count, name_nod_var(1))
  ! CALL ncderrcheck( __LINE__ ,stat )
  ! name_nod_var_start(1) = 1
  ! name_nod_var_start(2) = 2
  ! name_nod_var_count(1) = LEN(trim(name_nod_var(2)))
  ! name_nod_var_count(2) = 1
  ! stat = nf_put_vara_text(ofid, name_nod_var_id, name_nod_var_start, name_nod_var_count, name_nod_var(2))
  ! CALL ncderrcheck( __LINE__ ,stat )
  ! name_nod_var_start(1) = 1
  ! name_nod_var_start(2) = 3
  ! name_nod_var_count(1) = LEN(trim(name_nod_var(3)))
  ! name_nod_var_count(2) = 1
  ! stat = nf_put_vara_text(ofid, name_nod_var_id, name_nod_var_start, name_nod_var_count, name_nod_var(3))
  ! CALL ncderrcheck( __LINE__ ,stat )
  ! name_nod_var_start(1) = 1
  ! name_nod_var_start(2) = 4
  ! name_nod_var_count(1) = LEN(trim(name_nod_var(4)))
  ! name_nod_var_count(2) = 1
  ! stat = nf_put_vara_text(ofid, name_nod_var_id, name_nod_var_start, name_nod_var_count, name_nod_var(4))
  ! CALL ncderrcheck( __LINE__ ,stat )
  ! name_nod_var_start(1) = 1
  ! name_nod_var_start(2) = 5
  ! name_nod_var_count(1) = LEN(trim(name_nod_var(5)))
  ! name_nod_var_count(2) = 1
  ! stat = nf_put_vara_text(ofid, name_nod_var_id, name_nod_var_start, name_nod_var_count, name_nod_var(5))
  ! CALL ncderrcheck( __LINE__ ,stat )
  ! name_nod_var_start(1) = 1
  ! name_nod_var_start(2) = 6
  ! name_nod_var_count(1) = LEN(trim(name_nod_var(6)))
  ! name_nod_var_count(2) = 1
  ! stat = nf_put_vara_text(ofid, name_nod_var_id, name_nod_var_start, name_nod_var_count, name_nod_var(6))
  ! CALL ncderrcheck( __LINE__ ,stat )
  ! name_nod_var_start(1) = 1
  ! name_nod_var_start(2) = 7
  ! name_nod_var_count(1) = LEN(trim(name_nod_var(7)))
  ! name_nod_var_count(2) = 1
  ! stat = nf_put_vara_text(ofid, name_nod_var_id, name_nod_var_start, name_nod_var_count, name_nod_var(7))
  ! CALL ncderrcheck( __LINE__ ,stat )

  ! ! Close the file
  ! stat = NF_CLOSE(ofid);
  ! CALL ncderrcheck( __LINE__ ,stat )
  
  write(0,*)"Done reading Exodus"

  ! !================================================================================
  ! !
  ! ! Read the OpenFoam BC files and interpolate WRF data there
  ! !
  ! !================================================================================

  ! DO ibdy = 1,nbdys
  !    IF ( ibdy .EQ. INTERIOR .AND. .NOT. ic ) cycle ! short circuit if we do not want to generate interior file
  !    IF ( ibdy .NE. INTERIOR ) THEN
  !       CALL read_openfoam_bdy_coords(ibdy,TRIM(bdynames(ibdy))//'_bc.dat')
  !    ELSE
  !       CALL read_openfoam_bdy_coords(ibdy,TRIM(bdynames(ibdy))//'.dat')
  !    ENDIF
  !    IF ( ALLOCATED(bdy(ibdy)%point) ) THEN
  !       CALL precompute_openfoam_points(ibdy,xlat,xlong,ids,ide,jds,jde,ips,ipe,jps,jpe,ims,ime,jms,jme )
  !       sec = sec_of_day(TRIM(Times(it)))
  !       sec = sec - sec_start + sec_offset
  !       ! create the time directory if it doesn't already exist
  !       IF ( ibdy .NE. INTERIOR ) THEN
  !          IF ( sec > 999999 ) THEN
  !             WRITE(0,*)sec,' is too many seconds from start.'
  !             WRITE(0,*)'Use -offset argument to make this a six digit number.'
  !             CALL help
  !             STOP 99
  !          ENDIF
  !          WRITE(secstr,'(I6.1)')sec
  !          dirpath = TRIM(bdynames(ibdy))//"/"
  !       ELSE
  !          WRITE(secstr,'(I6.1)')sec
  !          dirpath = ""
  !       ENDIF

  !       comstr = "mkdir -p " // TRIM(dirpath)//TRIM(ADJUSTL(secstr))
  !       CALL system (TRIM(comstr))

  !       ! Write the header of the OpenFoam files
  !       outnameT = TRIM(dirpath)//TRIM(ADJUSTL(secstr))//"/T"
  !       OPEN( 76 , file=TRIM(outnameT), form="formatted" )
  !       CALL write_header_scalar( 76, bdy(ibdy)%npoints , ibdy.NE.INTERIOR , sec , "T" )

  !       outnamePd = TRIM(dirpath)//TRIM(ADJUSTL(secstr))//"/p_rgh"
  !       OPEN( 77 , file=TRIM(outnamePd), form="formatted" )
  !       CALL write_header_scalar( 77, bdy(ibdy)%npoints , ibdy.NE.INTERIOR , sec , "p_rgh" )

  !       outnameU = TRIM(dirpath)//TRIM(ADJUSTL(secstr))//"/U"
  !       OPEN( 78 , file=TRIM(outnameU), form="formatted" )
  !       CALL write_header_vector( 78, bdy(ibdy)%npoints , ibdy.NE.INTERIOR , sec , "U" )

  !       IF ( ibdy .EQ. BDY_ZS .AND. have_hfx .AND. use_hfx ) THEN
  !          outnameHFX = TRIM(dirpath)//TRIM(ADJUSTL(secstr))//"/qwall"
  !          OPEN( 79 , file=TRIM(outnameHFX), form="formatted" )
  !          CALL write_header_vector( 79, bdy(ibdy)%npoints , ibdy.NE.INTERIOR , sec , "qwall" )
  !       ENDIF

  !       ! Linking WRF and OpenFoam data
  !       DO ipoint = 1,bdy(ibdy)%npoints
  !          of_lat = bdy(ibdy)%point(ipoint)%lat
  !          of_lon = bdy(ibdy)%point(ipoint)%lon
  !          of_lz = bdy(ibdy)%point(ipoint)%lz
  !          j = bdy(ibdy)%point(ipoint)%j ! precomputed jcoord of cell center corresponding to lat
  !          i = bdy(ibdy)%point(ipoint)%i ! precomputed icoord of cell center corresponding to lon

  !          DO kk = 1,size(zz,3)
  !             DO jj = 0,1
  !                DO ii = 0,1
  !                   zzcol(ii,jj,kk)=zz(i+ii,j+jj,kk) - zz(i+ii,j+jj,1) ! zz is full height at cell centers
  !                   IF ( kk .LE. kpe-1 ) THEN
  !                      zcol (ii,jj,kk)= z(i+ii,j+jj,kk) - zz(i+ii,j+jj,1) ! z is half height at cell centers
  !                   ENDIF
  !                ENDDO
  !             ENDDO
  !          ENDDO

  !          ! find the level index of the openfoam point in WRF, both in the full-level
  !          ! and half-level ranges. Lowest index is closest to surface. Also store the
  !          ! indices for the 3 neighbors to the north, east, and northeast, since these
  !          ! are needed for horizontally interpolating in the finterp function
  !          DO jj = 0,1
  !             DO ii = 0,1
  !                IF (zzcol(ii,jj,1).LE.of_lz.AND.of_lz.LT.zcol(ii,jj,1))THEN ! special case, of_lz is below first half-level
  !                   kzz(ii,jj) = 1 ! ignore other special case since open foam wont go that high
  !                   kz(ii,jj) = 0
  !                ELSE
  !                   DO k = kps+1,kpe
  !                      IF (zzcol(ii,jj,k-1).LE.of_lz.AND.of_lz.LT.zzcol(ii,jj,k)) kzz(ii,jj) = k-1 ! full level
  !                      IF (k.LT.kpe) THEN
  !                         IF (zcol(ii,jj,k-1).LE.of_lz.AND.of_lz.LT.zcol(ii,jj,k)) kz(ii,jj) = k-1 ! half level
  !                      ENDIF
  !                   ENDDO
  !                ENDIF
  !             ENDDO
  !          ENDDO

  !          !variables on half-levels open foam coords dims of field dims of lat lon arrays
  !          u_new = finterp(u ,zcol ,xlat,xlong,kz ,of_lat,of_lon,of_lz,i,j,ips,ipe-1,jps,jpe-1,kps,kpe-1, ips,ipe-1,jps,jpe-1)
  !          v_new = finterp(v ,zcol ,xlat,xlong,kz ,of_lat,of_lon,of_lz,i,j,ips,ipe-1,jps,jpe-1,kps,kpe-1, ips,ipe-1,jps,jpe-1)
  !          t_new = finterp(t ,zcol ,xlat,xlong,kz ,of_lat,of_lon,of_lz,i,j,ips,ipe-1,jps,jpe-1,kps,kpe-1, ips,ipe-1,jps,jpe-1)
  !          pres_new = finterp(pres,zzcol,xlat,xlong,kzz,of_lat,of_lon,of_lz,i,j,ips,ipe-1,jps,jpe-1,kps,kpe-1, ips,ipe-1,jps,jpe-1)

  !          !variables on full-levels open foam coords dims of field dims of lat lon arrays
  !          w_new = finterp(w ,zzcol,xlat,xlong,kzz,of_lat,of_lon,of_lz,i,j,ips,ipe-1,jps,jpe-1,kps,kpe , ips,ipe-1,jps,jpe-1)
  !          t_ground = t(i,j,1)
  !          ! compute "pd" which is defined as pressure divided by density at surface minus geopotential
  !          ! that is, pd = p / rho - g*z . Note, however, that we don.t have density so compute density at
  !          ! surface as rho0 = p0 / (R*T0), where R is 286.9 and T0 is surface temp. Substituting for rho
  !          ! into the above, this becomes:
  !          pd = (pres_new*R_D*t_ground)/pres(i,j,1) - g * of_lz
  !          IF ( ibdy .EQ. BDY_ZS .AND. have_hfx .AND. use_hfx ) THEN
  !             kzz = 0 ! turn off vertical interpolation in call to finterp
  !             hfx_new = finterp(hfx,zzcol,xlat,xlong,kzz,of_lat,of_lon,of_lz,i,j,ips,ipe-1,jps,jpe-1,1,1, ips,ipe-1,jps,jpe-1)
  !             rho0 = pres(i,j,1) / ( R_D * t_ground )
  !             hfx_new = -( hfx_new / ( rho0 * CP ) )
  !          ENDIF
  !          WRITE(76,*) t_new
  !          WRITE(77,*) pd
  !          ! note that positive theta angle (computed above) implies WRF grid is rotated counterclockwise w.r.t. true
  !          ! see http://en.wikipedia.org/wiki/Rotation_matrix
  !          WRITE(78,'("(",f12.7," ",f12.7," ",f12.7,")")')u_new*costheta-v_new*sintheta,u_new*sintheta+v_new*costheta,w_new
  !          IF ( ibdy .EQ. BDY_ZS .AND. have_hfx .AND. use_hfx ) THEN
  !             WRITE(79,'("(",f12.7," ",f12.7," ",f12.7,")")')0., 0., hfx_new
  !          ENDIF
  !       ENDDO
  !       ! close the unit and then run a system command to strip off the first space that Fortran insists on writing
  !       CALL write_trailer_T(76,ibdy.EQ.INTERIOR) ; CLOSE( 76 ) ;
  !       comstr = "sed 's/^ //' "//TRIM(outnameT) //"> foo_ ; /bin/mv -f foo_ "// TRIM(outnameT)
  !       CALL system(TRIM(comstr))
  !       CALL write_trailer_pd(77,ibdy.EQ.INTERIOR) ; CLOSE( 77 ) ;
  !       comstr = "sed 's/^ //' "//TRIM(outnamePd) //"> foo_ ; /bin/mv -f foo_ "// TRIM(outnamePd)
  !       CALL system(TRIM(comstr))
  !       CALL write_trailer_U(78,ibdy.EQ.INTERIOR) ; CLOSE( 78 ) ;
  !       comstr = "sed 's/^ //' "//TRIM(outnameU) //"> foo_ ; /bin/mv -f foo_ "// TRIM(outnameU)
  !       CALL system(TRIM(comstr))
  !       IF ( ibdy .EQ. BDY_ZS .AND. have_hfx .AND. use_hfx ) THEN
  !          CALL write_trailer_HFX(79,.FALSE.) ; CLOSE( 79 ) ;
  !          comstr = "sed 's/^ //' "//TRIM(outnameU) //"> foo_ ; /bin/mv -f foo_ "// TRIM(outnameU)
  !          CALL system(TRIM(comstr))
  !       ENDIF
  !    ENDIF
  ! ENDDO

  ! Nice final message of congratulations
  write(*,*)'End program. Congratulations!'
  
END PROGRAM wrftonalu

!--------------------------------------------------------------------------------
!
!> @brief finterp interpolates WRF data at OpenFoam nodes
!> @param f        variable to interpolate
!> @param zcol     
!> @param lat      latitude
!> @param lon      longitude
!> @param kz       
!> @param of_lat   latitude of desired openfoam point
!> @param of_lon   longitude of desired openfoam point
!> @param of_lz    height of desired openfoam point
!> @param i        precomputed icoord of cell center corresponding to lon
!> @param j        precomputed jcoord of cell center corresponding to lat
!> @param is       start i index
!> @param ie       end i index   
!> @param js       start j index 
!> @param je       end j index   
!> @param ks       start k index 
!> @param ke       end k index   
!> @param ims
!> @param ime
!> @param jms
!> @param jme
!
!--------------------------------------------------------------------------------
REAL FUNCTION finterp( f, zcol, lat, lon, kz, of_lat, of_lon, of_lz, i, j, is,ie,js,je,ks,ke,ims,ime,jms,jme )
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: i,j,is,ie,js,je,ks,ke,ims,ime,jms,jme
  INTEGER, INTENT(IN) :: kz(0:1,0:1)
  REAL, INTENT(IN) :: f(is:ie,js:je,ks:ke), zcol(0:1,0:1,ks:ke)
  DOUBLE PRECISION, INTENT(IN) :: of_lat, of_lon, of_lz
  REAL, INTENT(IN) :: lat(ims:ime,jms:jme),lon(ims:ime,jms:jme)
  ! local
  INTEGER k
  REAL f00,f10,f01,f11,rm
  k = kz(0,0)
  IF ( k .GE. 1 ) THEN
     f00 = f(i ,j ,k) + (of_lz-zcol(0 ,0 ,k))*(f(i ,j ,k+1)-f(i ,j ,k))/(zcol(0 ,0 ,k+1)-zcol(0 ,0 ,k))
  ELSE
     f00 = f(i ,j ,1)
  ENDIF
  k = kz(1,0)
  IF ( k .GE. 1 ) THEN
     f10 = f(i+1,j ,k) + (of_lz-zcol(0+1,0 ,k))*(f(i+1,j ,k+1)-f(i+1,j ,k))/(zcol(0+1,0 ,k+1)-zcol(0+1,0 ,k))
  ELSE
     f10 = f(i+1,j ,1)
  ENDIF
  k = kz(0,1)
  IF ( k .GE. 1 ) THEN
     f01 = f(i ,j+1,k) + (of_lz-zcol(0 ,0+1,k))*(f(i ,j+1,k+1)-f(i ,j+1,k))/(zcol(0 ,0+1,k+1)-zcol(0 ,0+1,k))
  ELSE
     f01 = f(i ,j+1,1)
  ENDIF
  k = kz(1,1)
  IF ( k .GE. 1 ) THEN
     f11 = f(i+1,j+1,k) + (of_lz-zcol(0+1,0+1,k))*(f(i+1,j+1,k+1)-f(i+1,j+1,k))/(zcol(0+1,0+1,k+1)-zcol(0+1,0+1,k))
  ELSE
     f11 = f(i+1,j+1,1)
  ENDIF
  !
  rm = 1.0/((lon(i+1,j)-lon(i,j))*(lat(i,j+1)-lat(i,j)))
  finterp = f00*rm*(lon(i+1,j)-of_lon )*(lat(i,j+1)-of_lat ) + &
       f10*rm*(of_lon -lon(i,j))*(lat(i,j+1)-of_lat ) + &
       f01*rm*(lon(i+1,j)-of_lon )*(of_lat -lat(i,j)) + &
       f11*rm*(of_lon -lon(i,j))*(of_lat -lat(i,j))
  RETURN
END FUNCTION finterp


!--------------------------------------------------------------------------------
!
!> @brief Help function for command line usage
!
!--------------------------------------------------------------------------------
SUBROUTINE help
  IMPLICIT NONE
  CHARACTER(LEN=120) :: cmd
  CALL getarg(0, cmd)
  WRITE(*,'(/,"Usage: ", A, " ncdfile [ncdfiles*] [-startdate startdate [-offset offset]] [-ic]")') trim(cmd)
  WRITE(*,'("       startdate   date string of form yyyy-mm-dd_hh_mm_ss or yyyy-mm-dd_hh:mm:ss")')
  WRITE(*,'("       offset      number of seconds to start OpenFOAM directory naming (default 0)")')
  WRITE(*,'("       -ic         program should generate init conditions too")')
  WRITE(*,'("       -qwall      program should generate temp flux in lower bc file ",/)')
  STOP
END SUBROUTINE help


!--------------------------------------------------------------------------------
!
!> @brief Check for netCDF function error
!> @param lineno line number of error
!> @param stat error status
!
!--------------------------------------------------------------------------------
SUBROUTINE ncderrcheck( lineno, stat )
  IMPLICIT NONE
  INCLUDE 'netcdf.inc'
  INTEGER, INTENT(IN) :: lineno,stat
  IF ( stat .NE. NF_NOERR ) THEN
     WRITE(0,*)'Line ',lineno,NF_STRERROR(stat)
     STOP 99
  ENDIF
END SUBROUTINE ncderrcheck


!--------------------------------------------------------------------------------
!
!> @brief Read data in netCDF file
!> @param ctrl      is this field required?
!> @param ncids     netCDF file identifier
!> @param numfiles  number of files to read
!> @param vname     variable name to be read
!> @param buf       where to put the variable
!> @param itime     time step to read the variable
!> @param ndim      number of dimensions
!> @param ids       start i index
!> @param ide       end i index   
!> @param jds       start j index 
!> @param jde       end j index   
!> @param kds       start k index 
!> @param kde       end k index   
!
!--------------------------------------------------------------------------------
SUBROUTINE getvar_real(ctrl,ncids,numfiles,vname,buf,itime,ndim,ids,ide,jds,jde,kds,kde)
  IMPLICIT NONE
  INCLUDE 'netcdf.inc'
  INTEGER, INTENT(IN), DIMENSION(*) :: ncids
  INTEGER, INTENT(IN) :: numfiles
  LOGICAL, INTENT(INOUT):: ctrl
  REAL, INTENT(INOUT) :: buf(*)
  CHARACTER*(*), INTENT(IN) :: vname
  INTEGER, INTENT(IN) :: itime,ndim,ids,ide,jds,jde,kds,kde
  INTEGER strt(4),cnt(4)
  INTEGER stat,varid, i,ncid
  LOGICAL found
  !
  found = .FALSE.
  DO i = 1,numfiles
     ncid = ncids(i)
     IF ( ncid .GT. 0 .AND. .NOT. found ) THEN
        stat = NF_INQ_VARID(ncid,vname,varid)
        IF ( stat .EQ. 0 ) THEN
           strt = 1
           IF ( ndim .EQ. 3 ) THEN
              cnt(1) = ide-ids+1
              cnt(2) = jde-jds+1
              cnt(3) = kde-kds+1
              cnt(4) = itime
           ELSE
              cnt(1) = ide-ids+1
              cnt(2) = jde-jds+1
              cnt(3) = itime
           ENDIF
           stat = NF_GET_VARA_REAL(ncid,varid,strt,cnt,buf)
           IF ( stat .EQ. 0 ) found = .TRUE.
        ENDIF
     ENDIF
  ENDDO
  IF ( .NOT. found .AND. ctrl ) THEN
     WRITE(0,*)'getvar_real: did not find ',TRIM(vname),' in any input file'
     STOP 99
  ENDIF
  ctrl = found
  RETURN
END SUBROUTINE getvar_real


!--------------------------------------------------------------------------------
!
!> @brief Write OpenFoam header for a scalar field
!> @param unit        file output unit
!> @param n           number of nodes
!> @param bc          is this a boundary file or interior
!> @param location    location of variable
!> @param objectname  name of variable
!
!--------------------------------------------------------------------------------
SUBROUTINE write_header_scalar( unit, n, bc , location, objectname )
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: unit, n
  LOGICAL, INTENT(IN) :: bc ! is this a boundary file or interior
  INTEGER, INTENT(IN) :: location
  CHARACTER*(*), INTENT(IN) :: objectname
  WRITE(unit,*)'/*--------------------------------*- C++ -*----------------------------------*\'
  WRITE(unit,*)'| =========                 |                                                 |'
  WRITE(unit,*)'| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |'
  WRITE(unit,*)'|  \\    /   O peration     | Version:  1.6                                   |'
  WRITE(unit,*)'|   \\  /    A nd           | Web:      http://www.OpenFOAM.org               |'
  WRITE(unit,*)'|    \\/     M anipulation  |                                                 |'
  WRITE(unit,*)'\*---------------------------------------------------------------------------*/'
  WRITE(unit,*)'FoamFile'
  WRITE(unit,*)'{'
  WRITE(unit,*)'    version     2.0;'
  WRITE(unit,*)'    format      ascii;'
  IF ( bc ) THEN
     WRITE(unit,*)'    class       scalarAverageField;'
     WRITE(unit,*)'    object      values;'
  ELSE
     WRITE(unit,*)'    class       volScalarField;'
     WRITE(unit,'("    location    """,i4.4,""";")')location
     WRITE(unit,*)'    object      ',TRIM(objectname),' ;'
  ENDIF
  WRITE(unit,*)'}'
  WRITE(unit,*)'// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //'
  WRITE(unit,*)' '
  IF ( bc ) THEN
     WRITE(unit,*)'// Average'
     IF ( TRIM(objectname).EQ.'U' ) THEN
        WRITE(unit,*)'(0 0 0)'
     ELSE
        WRITE(unit,*)'0'
     ENDIF
     WRITE(unit,*)' '
     WRITE(unit,*)' '
  ELSE
     IF ( TRIM(objectname).EQ.'p_rgh' ) THEN
        WRITE(unit,*)'dimensions      [0 2 -2 0 0 0 0];' ! And this means m^2/s^2 !
     ELSE IF ( TRIM(objectname).EQ.'T' ) THEN
        WRITE(unit,*)'dimensions      [0 0 0 1 0 0 0];' ! This means degrees Kelvin !
     ELSE
        WRITE(0,*)'write_header_scalar does not know about this variable: ',TRIM(objectname)
        STOP 122
     ENDIF
     WRITE(unit,*)' '
     WRITE(unit,*)'internalField   nonuniform List<scalar>'
  ENDIF
  WRITE(unit,'(I10)')n
  WRITE(unit,*)'('
END SUBROUTINE write_header_scalar


!--------------------------------------------------------------------------------
!
!> @brief Write OpenFoam header for a vector field
!> @param unit        file output unit
!> @param n           number of nodes
!> @param bc          is this a boundary file or interior
!> @param location    location of variable
!> @param objectname  name of variable
!
!--------------------------------------------------------------------------------
SUBROUTINE write_header_vector( unit, n, bc, location, objectname )
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: unit, n
  LOGICAL, INTENT(IN) :: bc ! boundary file or interior
  INTEGER, INTENT(IN) :: location
  CHARACTER*(*), INTENT(IN) :: objectname
  WRITE(unit,*)'/*--------------------------------*- C++ -*----------------------------------*\'
  WRITE(unit,*)'| =========                 |                                                 |'
  WRITE(unit,*)'| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |'
  WRITE(unit,*)'|  \\    /   O peration     | Version:  1.6                                   |'
  WRITE(unit,*)'|   \\  /    A nd           | Web:      http://www.OpenFOAM.org               |'
  WRITE(unit,*)'|    \\/     M anipulation  |                                                 |'
  WRITE(unit,*)'\*---------------------------------------------------------------------------*/'
  WRITE(unit,*)'FoamFile'
  WRITE(unit,*)'{'
  WRITE(unit,*)'    version     2.0;'
  WRITE(unit,*)'    format      ascii;'
  IF ( bc ) THEN
     WRITE(unit,*)'    class       vectorAverageField;'
     WRITE(unit,*)'    object      values;'
  ELSE
     WRITE(unit,*)'    class       volVectorField;'
     WRITE(unit,'("    location    """,i4.4,""";")')location
     WRITE(unit,*)'    object      ',TRIM(objectname),';'
  ENDIF
  WRITE(unit,*)'}'
  WRITE(unit,*)'// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //'
  WRITE(unit,*)' '
  IF ( bc ) THEN
     WRITE(unit,*)'// Average'
     WRITE(unit,*)'(0 0 0)'
     WRITE(unit,*)' '
     WRITE(unit,*)' '
  ELSE
     IF ( TRIM(objectname).EQ.'U' ) THEN
        WRITE(unit,*)'dimensions      [0 1 -1 0 0 0 0];' ! this means m/s
     ELSE IF ( TRIM(objectname).EQ.'qwall' ) THEN
        WRITE(unit,*)'dimensions      [0 1 -1 1 0 0 0];' ! this means ?
     ENDIF
     WRITE(unit,*)' '
     WRITE(unit,*)'internalField   nonuniform List<vector>'
  ENDIF
  WRITE(unit,'(I10)')n
  WRITE(unit,*)'('
END SUBROUTINE write_header_vector

!--------------------------------------------------------------------------------
!
!> @brief Write OpenFoam file trailer for variable U
!> @param unit        file output unit
!> @param interior    is this a boundary file or interior
!
!--------------------------------------------------------------------------------
SUBROUTINE write_trailer_U ( unit, interior )
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: unit
  LOGICAL, INTENT(IN) :: interior
  WRITE(unit,*)')'
  WRITE(unit,*)';'
  IF ( interior ) THEN
     WRITE(unit,*)
     WRITE(unit,*)'boundaryField'
     WRITE(unit,*)'{'
     WRITE(unit,*)'    lower'
     WRITE(unit,*)'    {'
     WRITE(unit,*)'        type            velocityABLWallFunction;'
     WRITE(unit,*)'        print           1;'
     WRITE(unit,*)'        U               U;'
     WRITE(unit,*)'        value           uniform (0 0 0);'
     WRITE(unit,*)'    }'
     WRITE(unit,*)'    upper'
     WRITE(unit,*)'    {'
     WRITE(unit,*)'        type            timeVaryingMappedFixedValue;'
     WRITE(unit,*)'        setAverage      0;'
     WRITE(unit,*)'        offset          (0 0 0);'
     WRITE(unit,*)'    }'
     WRITE(unit,*)'    west'
     WRITE(unit,*)'    {'
     WRITE(unit,*)'        type            timeVaryingMappedFixedValue;'
     WRITE(unit,*)'        setAverage      0;'
     WRITE(unit,*)'        offset          (0 0 0);'
     WRITE(unit,*)'    }'
     WRITE(unit,*)'    east'
     WRITE(unit,*)'    {'
     WRITE(unit,*)'        type            timeVaryingMappedFixedValue;'
     WRITE(unit,*)'        setAverage      0;'
     WRITE(unit,*)'        offset          (0 0 0);'
     WRITE(unit,*)'    }'
     WRITE(unit,*)'    south'
     WRITE(unit,*)'    {'
     WRITE(unit,*)'        type            timeVaryingMappedFixedValue;'
     WRITE(unit,*)'        setAverage      0;'
     WRITE(unit,*)'        offset          (0 0 0);'
     WRITE(unit,*)'    }'
     WRITE(unit,*)'    north'
     WRITE(unit,*)'    {'
     WRITE(unit,*)'        type            timeVaryingMappedFixedValue;'
     WRITE(unit,*)'        setAverage      0;'
     WRITE(unit,*)'        offset          (0 0 0);'
     WRITE(unit,*)'    }'
     WRITE(unit,*)'}'
     WRITE(unit,*)
     WRITE(unit,*)
  ENDIF
  WRITE(unit,*)
  WRITE(unit,*)'// ************************************************************************* //'
  RETURN
END SUBROUTINE write_trailer_U


!--------------------------------------------------------------------------------
!
!> @brief Write OpenFoam file trailer for variable T
!> @param unit        file output unit
!> @param interior    is this a boundary file or interior
!
!--------------------------------------------------------------------------------
SUBROUTINE write_trailer_T ( unit, interior )
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: unit
  LOGICAL, INTENT(IN) :: interior
  WRITE(unit,*)')'
  WRITE(unit,*)';'
  WRITE(unit,*)
  IF ( interior ) THEN
     WRITE(unit,*)'boundaryField'
     WRITE(unit,*)'{'
     WRITE(unit,*)'    lower'
     WRITE(unit,*)'    {'
     WRITE(unit,*)'        type            zeroGradient;'
     WRITE(unit,*)'    }'
     WRITE(unit,*)'    upper'
     WRITE(unit,*)'    {'
     WRITE(unit,*)'        type            timeVaryingMappedFixedValue;'
     WRITE(unit,*)'        setAverage      0;'
     WRITE(unit,*)'        offset          0;'
     WRITE(unit,*)'    }'
     WRITE(unit,*)'    west'
     WRITE(unit,*)'    {'
     WRITE(unit,*)'        type            timeVaryingMappedFixedValue;'
     WRITE(unit,*)'        setAverage      0;'
     WRITE(unit,*)'        offset          0;'
     WRITE(unit,*)'    }'
     WRITE(unit,*)'    east'
     WRITE(unit,*)'    {'
     WRITE(unit,*)'        type            timeVaryingMappedFixedValue;'
     WRITE(unit,*)'        setAverage      0;'
     WRITE(unit,*)'        offset          0;'
     WRITE(unit,*)'    }'
     WRITE(unit,*)'    south'
     WRITE(unit,*)'    {'
     WRITE(unit,*)'        type            timeVaryingMappedFixedValue;'
     WRITE(unit,*)'        setAverage      0;'
     WRITE(unit,*)'        offset          0;'
     WRITE(unit,*)'    }'
     WRITE(unit,*)'    north'
     WRITE(unit,*)'    {'
     WRITE(unit,*)'        type            timeVaryingMappedFixedValue;'
     WRITE(unit,*)'        setAverage      0;'
     WRITE(unit,*)'        offset          0;'
     WRITE(unit,*)'    }'
     WRITE(unit,*)'}'
     WRITE(unit,*)
  ENDIF
  WRITE(unit,*)
  WRITE(unit,*)'// ************************************************************************* //'
  RETURN
END SUBROUTINE write_trailer_T


!--------------------------------------------------------------------------------
!
!> @brief Write OpenFoam file trailer for variable pd
!> @param unit        file output unit
!> @param interior    is this a boundary file or interior
!
!--------------------------------------------------------------------------------
SUBROUTINE write_trailer_pd( unit, interior )
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: unit
  LOGICAL, INTENT(IN) :: interior
  WRITE(unit,*)')'
  WRITE(unit,*)';'
  WRITE(unit,*)
  IF ( interior ) THEN
     WRITE(unit,*)'boundaryField'
     WRITE(unit,*)'{'
     WRITE(unit,*)'    lower'
     WRITE(unit,*)'    {'
     WRITE(unit,*)'        type            buoyantPressureMod;'
     WRITE(unit,*)'        rho             rhok;'
     WRITE(unit,*)'        value           uniform 0;'
     WRITE(unit,*)'    }'
     WRITE(unit,*)'    upper'
     WRITE(unit,*)'    {'
     WRITE(unit,*)'        type            buoyantPressureMod;'
     WRITE(unit,*)'        rho             rhok;'
     WRITE(unit,*)'        value           uniform 0;'
     WRITE(unit,*)'    }'
     WRITE(unit,*)'    east'
     WRITE(unit,*)'    {'
     WRITE(unit,*)'        type            buoyantPressureMod;'
     WRITE(unit,*)'        rho             rhok;'
     WRITE(unit,*)'        value           uniform 0;'
     WRITE(unit,*)'    }'
     WRITE(unit,*)'    west'
     WRITE(unit,*)'    {'
     WRITE(unit,*)'        type            buoyantPressureMod;'
     WRITE(unit,*)'        rho             rhok;'
     WRITE(unit,*)'        value           uniform 0;'
     WRITE(unit,*)'    }'
     WRITE(unit,*)'    south'
     WRITE(unit,*)'    {'
     WRITE(unit,*)'        type            buoyantPressureMod;'
     WRITE(unit,*)'        rho             rhok;'
     WRITE(unit,*)'        value           uniform 0;'
     WRITE(unit,*)'    }'
     WRITE(unit,*)'    north'
     WRITE(unit,*)'    {'
     WRITE(unit,*)'        type            buoyantPressureMod;'
     WRITE(unit,*)'        rho             rhok;'
     WRITE(unit,*)'        value           uniform 0;'
     WRITE(unit,*)'    }'
     WRITE(unit,*)'}'
     WRITE(unit,*)
  ENDIF
  WRITE(unit,*)
  WRITE(unit,*)'// ************************************************************************* //'
  RETURN
END SUBROUTINE write_trailer_pd


!--------------------------------------------------------------------------------
!
!> @brief Write OpenFoam file trailer for variable HFX (heat flux)
!> @param unit        file output unit
!> @param interior    is this a boundary file or interior
!
!--------------------------------------------------------------------------------
SUBROUTINE write_trailer_HFX ( unit, interior )
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: unit
  LOGICAL, INTENT(IN) :: interior
  WRITE(unit,*)')'
  WRITE(unit,*)';'
  WRITE(unit,*)
  WRITE(unit,*)'// ************************************************************************* //'
  RETURN
END SUBROUTINE write_trailer_HFX

!--------------------------------------------------------------------------------
!
!> @brief Find the second of the day
!> @param s        seconds
!> @warning OVERLY SIMPLE -- assumes same day! and assumes 19 char WRF style date str
!
!--------------------------------------------------------------------------------
INTEGER FUNCTION sec_of_day ( s )
  IMPLICIT NONE
  CHARACTER*(*), INTENT(IN) :: s
  INTEGER hh,mm,ss
  ! 0000000001111111111
  ! 1234567890123456789
  ! 2005-01-15_02_04_31
  READ(s(12:13),*)hh
  READ(s(15:16),*)mm
  READ(s(18:19),*)ss
  sec_of_day = hh*3600 + mm*60 + ss
END FUNCTION sec_of_day


!--------------------------------------------------------------------------------
!
!> @brief Check if this is a valid date
!> @param s    seconds
!
!--------------------------------------------------------------------------------
LOGICAL FUNCTION valid_date ( s )
  IMPLICIT NONE
  CHARACTER*(*), INTENT(IN) :: s
  LOGICAL, EXTERNAL :: isnum
  LOGICAL retval
  retval = .FALSE.
  IF ( LEN(TRIM(s)) .EQ. 19 ) THEN
     IF ( isnum(1,s) .AND. isnum(2,s) .AND. isnum(3,s) .AND. isnum(4,s) .AND. &
          s(5:5).EQ.'-' .AND. &
          isnum(6,s) .AND. isnum(7,s) .AND. &
          s(8:8).EQ.'-' .AND. &
          isnum(9,s) .AND. isnum(10,s) .AND. &
          s(11:11).EQ.'_' .AND. &
          isnum(12,s) .AND. isnum(13,s) .AND. &
          (s(14:14).EQ.'_' .OR. s(14:14).EQ.':') .AND. &
          isnum(15,s) .AND. isnum(16,s) .AND. &
          (s(17:17).EQ.'_' .OR. s(17:17).EQ.':') .AND. &
          isnum(18,s) .AND. isnum(19,s) ) THEN
        retval = .TRUE.
     ENDIF
  ENDIF
  valid_date = retval
END FUNCTION valid_date


!--------------------------------------------------------------------------------
!
!> @brief Not sure what this does
!
!--------------------------------------------------------------------------------
LOGICAL FUNCTION isnum ( i, str )
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: i
  CHARACTER*(*), INTENT(IN) :: str
  isnum = (ICHAR('0').LE. ICHAR(str(i:i)).AND.ICHAR(str(i:i)) .LE. ICHAR('9'))
END FUNCTION isnum
