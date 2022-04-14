! returns drifts and bounce time per well
subroutine w_bounce(q0,L_tot,b_arr,dbdx_arr,dbdy_arr,sqrtg_arr,theta_arr,lam, &
                    Delta_x,Delta_y,num_wells,w_psi_arr,w_alpha_arr,G_arr, &
                    N,bounce_idx,bounce_arr)
  ! Returns w_alpha, w_psi, and G.
  implicit none
  real(kind=8), dimension(N), intent(in)  :: b_arr,dbdx_arr,dbdy_arr,sqrtg_arr,theta_arr
  real(kind=8), intent(in)                :: q0, lam,Delta_x,Delta_y,L_tot
  integer,      intent(in)                :: num_wells, N
  real(kind=8), dimension(N)              :: h_arr
  real(kind=8), dimension(num_wells,2)    :: bounce_arr
  integer,      dimension(num_wells,2)    :: bounce_idx
  real(kind=8), dimension(num_wells)      :: denom_arr, numer_arr_alpha, numer_arr_psi
  real(kind=8), dimension(num_wells), intent(out) :: w_psi_arr,w_alpha_arr,G_arr


  ! Let us do the denominator first
  ! print *,'calculating denom'
  h_arr = q0 * b_arr * sqrtg_arr
  CALL bounce_average(h_arr,b_arr,theta_arr,lam,num_wells,bounce_idx, &
                      bounce_arr,denom_arr,N)
  ! Now construct omega_alpha and omega_psi
  ! print *,'calculating alpha num'
  h_arr = lam * Delta_x * dbdx_arr * q0 * b_arr * sqrtg_arr
  CALL bounce_average(h_arr,b_arr,theta_arr,lam,num_wells,bounce_idx, &
                      bounce_arr,numer_arr_alpha,N)

  ! print *,'calculating psi num'
  h_arr = -1.0 * lam * Delta_y * dbdy_arr * q0 * b_arr * sqrtg_arr
  CALL bounce_average(h_arr,b_arr,theta_arr,lam,num_wells,bounce_idx, &
                      bounce_arr,numer_arr_psi,N)


  ! Make arrays for w_psi, w_alpha, and G
  w_psi_arr   = numer_arr_psi   / denom_arr
  w_alpha_arr = numer_arr_alpha / denom_arr
  G_arr       = denom_arr       / L_tot
end subroutine w_bounce
