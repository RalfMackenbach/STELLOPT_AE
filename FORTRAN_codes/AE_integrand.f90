! integrand of AE
subroutine AE_integrand(num_wells,dlnTdx,dlnndx,Delta_x,z,w_alpha,w_psi,G,AE)
  ! Given array of w_alpha, w_psi, and G (decided by lambda), and w_dia and z
  ! return the AE integrand at that (lambda,z).
  implicit none
  integer, intent(in)                                :: num_wells
  real(kind=8), dimension(num_wells), intent(in)     :: w_alpha, w_psi, G
  real(kind=8), intent(in)                           :: dlnTdx, dlnndx, Delta_x, z
  real(kind=8)                                       :: w_dia
  real(kind=8), intent(out)                          :: AE
  ! problem with variable size again, unsure how this is handled.
  w_dia = Delta_x * ( dlnndx/z + dlnTdx * ( 1.0 - 3.0 / (2.0 * z) ) )
  AE = SUM((G*(w_psi*(-1 + Sqrt((-w_alpha + w_dia)**2 + w_psi**2)/ &
       Sqrt(w_alpha**2 + w_psi**2)) + w_alpha**2*(-1 + w_dia/w_alpha + &
       Sqrt((-w_alpha + w_dia)**2 + w_psi**2)/ &
       Sqrt(w_alpha**2 + w_psi**2)))*z**2.5)/exp(z))
end subroutine AE_integrand
