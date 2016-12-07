!------------------------------------------------------------------------------
!
! MODULE: module_dm
!
!> @author
!> J. Michalakes and M. Churchfield, National Renewable Energy Laboratory
!
!> @brief
!> Single-process stubs that would be actual parallel operations
!> in an MPI implementation if you want real versions of these,
!> link them in from WRF
!
!> @date 01/12/2016 J. Michalakes and M. Churchfield
!> - Initial version from WRFTOOOF
! 
!------------------------------------------------------------------------------


MODULE module_dm

CONTAINS
  
  !> @brief wrf_dm_max_real
  !> @param inval real
  REAL FUNCTION wrf_dm_max_real ( inval )
    IMPLICIT NONE
    REAL inval
    wrf_dm_max_real = inval
  END FUNCTION wrf_dm_max_real
  
END MODULE module_dm

!> @brief wrf_dm_bcast_integer
!> @param buf buffer
!> @param n1 integer
SUBROUTINE wrf_dm_bcast_integer( buf, n1 )
  IMPLICIT NONE
  INTEGER n1
  INTEGER buf(*)
  RETURN
END SUBROUTINE wrf_dm_bcast_integer

!> @brief wrf_dm_bcast_real
!> @param buf buffer
!> @param n1 integer
SUBROUTINE wrf_dm_bcast_real( buf, n1 )
  IMPLICIT NONE
  INTEGER n1
  REAL buf(*)
  RETURN
END SUBROUTINE wrf_dm_bcast_real

!> @brief wrf_dm_bcast_double
!> @param buf buffer
!> @param n1 integer
SUBROUTINE wrf_dm_bcast_double( buf, n1 )
  IMPLICIT NONE
  INTEGER n1
  REAL buf(*)
  RETURN
END SUBROUTINE wrf_dm_bcast_double

!> @brief wrf_dm_on_monitor
LOGICAL FUNCTION wrf_dm_on_monitor()
  wrf_dm_on_monitor = .true.
END FUNCTION wrf_dm_on_monitor
