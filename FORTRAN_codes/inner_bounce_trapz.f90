



! trapz-like estimate of bounce-ave function
subroutine inner_bounce_trapz(h_i,h_j,f_i,f_j,theta_i,theta_j,y)
  ! Estimation of integral of h/sqrt(f) dtheta, where we assume h and f to be
  ! well approximated by their linear interpolation between theta_i and theta_j.
  ! Only estimate the INNER part (i.e. not accounting for the edges where f
  ! vanshises). These are handled seperately.
  implicit none
  real(kind=8), dimension(:), intent(in) :: h_i
  real(kind=8), dimension(:), intent(in) :: h_j
  real(kind=8), dimension(:), intent(in) :: f_i
  real(kind=8), dimension(:), intent(in) :: f_j
  real(kind=8), dimension(:), intent(in) :: theta_i
  real(kind=8), dimension(:), intent(in) :: theta_j
  real(kind=8), intent(out) :: y
  y = SUM((-2*(Sqrt(f_j)*(2*h_i + h_j) + Sqrt(f_i)*(h_i + 2*h_j))* &
          (theta_i - theta_j))/(3.*(f_i + f_j + 2*Sqrt(f_i*f_j))))
end subroutine inner_bounce_trapz