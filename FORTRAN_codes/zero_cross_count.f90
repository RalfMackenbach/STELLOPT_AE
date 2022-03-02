


! find zero crossings
subroutine zero_cross_count(y_arr,N,num_cross)
  ! Finds all zeros of y. Based on "Numerical Recipes in Fortran 90" (1996),
  ! chapter B9, p.1184.
  implicit none
  integer, intent(in)                     :: N
  real(kind=8), dimension(N), intent(in)  :: y_arr
  real(kind=8), dimension(size(y_arr))    :: roll_arr
  logical, dimension(size(y_arr))         :: mask
  integer, intent(out) :: num_cross
  integer :: ix


  ! print *,y_arr
  ! Make shifted array (assumes periodic boundary conditions)
  roll_arr  = CSHIFT(y_arr,1)


  ! check zero crossings, returns TRUE left of the crossing. Also, count number
  ! of crossings
  mask      = y_arr*roll_arr .LE. 0
  num_cross = COUNT(mask)
end subroutine zero_cross_count
