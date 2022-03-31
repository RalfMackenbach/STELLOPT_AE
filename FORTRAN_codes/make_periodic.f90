! Makes periodic array (just appends first value to the end)
subroutine make_per(N,Delta_theta,b_arr,dbdx_arr,dbdy_arr,sqrtg_arr,theta_arr, &
                    b_arr_p, dbdx_arr_p, dbdy_arr_p, sqrtg_arr_p, theta_arr_p)
  integer, intent(in)                       :: N
  real(kind=8), intent(in)                  :: Delta_theta
  real(kind=8), dimension(N), intent(in)    :: b_arr, dbdx_arr, dbdy_arr, sqrtg_arr, theta_arr
  real(kind=8), dimension(N+1), intent(out) :: b_arr_p, dbdx_arr_p, dbdy_arr_p, sqrtg_arr_p, theta_arr_p

  ! Set all but last index equal
  b_arr_p(1:N)    = b_arr
  dbdx_arr_p(1:N) = dbdx_arr
  dbdy_arr_p(1:N) = dbdy_arr
  sqrtg_arr_p(1:N)= sqrtg_arr
  theta_arr_p(1:N)= theta_arr

  ! Now set last value equal to first
  b_arr_p(N+1)    = b_arr(1)
  dbdx_arr_p(N+1) = dbdx_arr(1)
  dbdy_arr_p(N+1) = dbdy_arr(1)
  sqrtg_arr_p(N+1)= sqrtg_arr(1)
  ! Special case for theta_arr. Last val is simply increased by Delta_theta
  ! True results is limit Delta theta -> 0^(+)
  theta_arr_p(N+1)= theta_arr(N) + Delta_theta
end subroutine make_per
