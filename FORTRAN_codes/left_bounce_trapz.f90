
! trapz-like estimate of bounce-ave function, left edge
subroutine left_bounce_trapz(h_l,h_0,f_l,f_0,theta_l,theta_0,y)
  ! Estimation of integral of edge integral h/sqrt(f) dtheta, from theta_l where
  ! f = 0, to the first theta node to its right theta_0
  implicit none
  real(kind=8), intent(in) :: h_l
  real(kind=8), intent(in) :: h_0
  real(kind=8), intent(in) :: f_l
  real(kind=8), intent(in) :: f_0
  real(kind=8), intent(in) :: theta_l
  real(kind=8), intent(in) :: theta_0
  real(kind=8), intent(out) :: y
  y = (2*(h_0 - 4*h_l)*(-theta_0 + theta_l))/(3.*Sqrt(f_0))
end subroutine left_bounce_trapz
