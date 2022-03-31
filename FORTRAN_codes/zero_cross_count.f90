


! find zero crossings
subroutine zero_cross_count(y_arr,N,num_cross)
  ! Finds all zeros of y. Based on "Numerical Recipes in Fortran 90" (1996),
  ! chapter B9, p.1184.
  implicit none
  integer, intent(in)                     :: N
  real(kind=8), dimension(N), intent(in)  :: y_arr
  real(kind=8), dimension(N-1)            :: l_arr, r_arr
  logical, dimension(N-1)                 :: mask
  integer, intent(out)                    :: num_cross
  integer :: ix


  ! print *,y_arr
  ! Make shifted array (assumes periodic boundary conditions)
  l_arr  = y_arr(1:(N-1))
  r_arr  = y_arr(2:N)

  ! check zero crossings, returns TRUE left of the crossing. Also, count number
  ! of crossings
  mask      = l_arr*r_arr .LE. 0
  num_cross = COUNT(mask)
end subroutine zero_cross_count
