


! find zero crossings
subroutine zero_cross_idx(y_arr,num_cross,N,zero_idx)
  ! Finds all zeros of y. Based on "Numerical Recipes in Fortran 90" (1996),
  ! chapter B9, p.1184.
  implicit none
  integer, intent(in)                     :: num_cross, N
  real(kind=8), dimension(N), intent(in)  :: y_arr
  real(kind=8), dimension(N)              :: roll_arr
  logical, dimension(N)                   :: mask
  integer, dimension(num_cross), intent(out) :: zero_idx
  integer :: ix


  ! print *,y_arr
  ! Make shifted array (assumes periodic boundary conditions)
  roll_arr  = CSHIFT(y_arr,1)
  
  ! check zero crossings, returns TRUE left of the crossing. Also, count number
  ! of crossings
  mask      = y_arr*roll_arr .LE. 0
  ! Return indices where mask is true, so left of crossings
  zero_idx = PACK([(ix,ix=1,SIZE(mask))],mask)
end subroutine zero_cross_idx
