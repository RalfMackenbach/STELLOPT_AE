
! trapz-like estimate of bounce-ave function, right edge
subroutine right_bounce_trapz(h_n,h_r,f_n,f_r,theta_n,theta_r,y)
  ! Estimation of integral of edge integral h/sqrt(f) dtheta, from theta_n to theta_r
  ! where f = 0
  implicit none
  real(kind=8), intent(in) :: h_n
  real(kind=8), intent(in) :: h_r
  real(kind=8), intent(in) :: f_n
  real(kind=8), intent(in) :: f_r
  real(kind=8), intent(in) :: theta_n
  real(kind=8), intent(in) :: theta_r
  real(kind=8), intent(out):: y
  y = (2*(h_n + 2*h_r)*(-theta_n + theta_r))/(3.*Sqrt(f_n))
end subroutine right_bounce_trapz
