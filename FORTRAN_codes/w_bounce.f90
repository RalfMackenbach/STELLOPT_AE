! returns drifts and bounce time per well
subroutine w_bounce(b_arr,dbdx_arr,dbdy_arr,sqrtg_arr,theta_arr,L_tot,lam, &
                    Delta_x,Delta_y,num_wells,w_psi_arr,w_alpha_arr,G_arr)
  ! Returns w_alpha, w_psi, and G.
  implicit none
  real(kind=8), dimension(:), intent(in)  :: b_arr,dbdx_arr,dbdy_arr,sqrtg_arr,theta_arr
  real(kind=8), intent(in)                :: L_tot,lam,Delta_x,Delta_y
  integer,      intent(in)                :: num_wells
  real(kind=8), dimension(size(b_arr))    :: h_arr
  real(kind=8), dimension(:), allocatable :: denom_arr, numer_arr_alpha, numer_arr_psi
  real(kind=8), dimension(:), allocatable, intent(out) :: w_psi_arr,w_alpha_arr,G_arr


  ! Allocate
  allocate(w_psi_arr(num_wells),w_alpha_arr(num_wells),G_arr(num_wells), &
          denom_arr(num_wells),numer_arr_alpha(num_wells),numer_arr_psi(num_wells))

  ! Let us do the denominator first
  ! Make array with ones
  h_arr = b_arr * sqrtg_arr
  CALL bounce_average(h_arr,b_arr,theta_arr,lam,denom_arr)

  ! Now construct omega_alpha and omega_psi
  h_arr = lam * Delta_x * dbdx_arr * b_arr * sqrtg_arr
  CALL bounce_average(h_arr,b_arr,theta_arr,lam,numer_arr_alpha)
  h_arr = -1.0 * lam * Delta_y * dbdy_arr * b_arr * sqrtg_arr
  CALL bounce_average(h_arr,b_arr,theta_arr,lam,numer_arr_psi)

  ! Make arrays for w_psi, w_alpha, and G
  w_psi_arr   = numer_arr_psi   / denom_arr
  w_alpha_arr = numer_arr_alpha / denom_arr
  G_arr       = denom_arr       / L_tot
end subroutine w_bounce