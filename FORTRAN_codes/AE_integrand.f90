! integrand of AE as in paper ( https://arxiv.org/pdf/2109.01042.pdf )
subroutine AE_integrand(num_wells,dlnTdx,dlnndx,Delta_x,z,walpha,wpsi,G,AE)
  ! Given array of w_alpha, w_psi, and G (decided by lambda), and w_dia and z
  ! return the AE integrand at that (lambda,z).
  implicit none
  integer, intent(in)                                :: num_wells
  real(kind=8), dimension(num_wells), intent(in)     :: walpha, wpsi, G
  real(kind=8), intent(in)                           :: dlnTdx, dlnndx, Delta_x, z
  real(kind=8)                                       :: wdia
  real(kind=8), intent(out)                          :: AE
  ! problem with variable size again, unsure how this is handled.
  wdia = Delta_x * ( dlnndx/z + dlnTdx * ( 1.0 - 3.0 / (2.0 * z) ) )
  AE = SUM((G*(walpha*(-walpha + wdia) - wpsi**2 + &
       Sqrt(walpha**2 + wpsi**2)*Sqrt((walpha - wdia)**2 + wpsi**2))*z**2.5)/exp(z))
end subroutine AE_integrand
