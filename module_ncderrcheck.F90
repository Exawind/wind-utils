!------------------------------------------------------------------------------
!
! MODULE: module_ncderrcheck
!
!> @author
!> J. Michalakes and M. Churchfield, National Renewable Energy Laboratory
!
!> @brief
!> A netCDF error checker
! 
!------------------------------------------------------------------------------


module module_ncderrcheck
contains
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

end module module_ncderrcheck
