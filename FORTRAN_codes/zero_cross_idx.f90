


! find zero crossings
subroutine zero_cross_idx(y_arr,zero_idx,num_cross)
  ! Finds all zeros of y. Based on "Numerical Recipes in Fortran 90" (1996),
  ! chapter B9, p.1184.
  implicit none
  real(kind=8), dimension(:), intent(in)  :: y_arr
  real :: lam
  real(kind=8), dimension(size(y_arr))    :: roll_arr
  logical, dimension(size(y_arr))         :: mask
  integer, intent(out) :: num_cross
  real, dimension (:), allocatable, intent(out) :: zero_idx
  integer :: ix

  ! Make shifted array (assumes periodic boundary conditions)
  roll_arr  = CSHIFT(y_arr,1)

  ! check zero crossings, returns TRUE left of the crossing. Also, count number
  ! of crossings
  mask      = y_arr*roll_arr <= 0
  num_cross = COUNT(mask)

  ! Allocate arrays
  allocate ( zero_idx(num_cross) )

  ! Return indices where mask is true, so left of crossings
  zero_idx = PACK([(ix,ix=1,size(mask))],mask)
end subroutine zero_cross_idx