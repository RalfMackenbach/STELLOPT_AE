


! find zero crossings
subroutine zero_cross_idx(y_arr,num_cross,N,zero_idx)
  ! Finds all zeros of y. Based on "Numerical Recipes in Fortran 90" (1996),
  ! chapter B9, p.1184.
  implicit none
  integer, intent(in)                     :: num_cross, N
  real(kind=8), dimension(N), intent(in)  :: y_arr
  real(kind=8), dimension(N-1)            :: l_arr, r_arr
  logical, dimension(N-1)                 :: mask
  integer, dimension(num_cross), intent(out) :: zero_idx
  integer :: ix


  ! Make shifted array (assumes periodic boundary conditions)
  l_arr  = y_arr(1:(N-1))
  r_arr  = y_arr(2:N)

  ! check zero crossings, returns TRUE left of the crossing. Also, count number
  ! of crossings
  mask      = l_arr*r_arr .LE. 0
  ! Return indices where mask is true, so left of crossings
  zero_idx = PACK([(ix,ix=1,SIZE(mask))],mask)
end subroutine zero_cross_idx
