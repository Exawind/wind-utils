!------------------------------------------------------------------------------
!
! MODULE: module_str2int
!
!> @author
!> Alexander Vogt
!
!> @brief
!> Convert a string to an integer variable
!
!> @details
!> This is from
!> href="http://stackoverflow.com/questions/24071722/converting-a-string-to-an-integer-in-fortran-90">Stack
!> Overflow</a>
!
! 
!------------------------------------------------------------------------------


module module_str2int
contains

  !--------------------------------------------------------------------------------
  !
  !> @brief Convert string to integer
  !
  !> @param[in]  str     input string
  !> @param[out] int     output integer
  !> @param[out] stat    output status
  !
  !--------------------------------------------------------------------------------
  elemental subroutine str2int(str,int,stat)
    implicit none
    ! Arguments
    character(len=*),intent(in) :: str
    integer,intent(out)         :: int
    integer,intent(out)         :: stat

    read(str,*,iostat=stat)  int
  end subroutine str2int

end module module_str2int
