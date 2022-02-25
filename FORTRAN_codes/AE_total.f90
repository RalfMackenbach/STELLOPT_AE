! Calculate total AE
subroutine AE_total(dlnTdx,dlnndx,Delta_x,Delta_y,b_arr,dbdx_arr,dbdy_arr, &
                    sqrtg_arr,theta_arr,lam_res,z_res,z_min,z_max,AE_tot)
  implicit none
  real(kind=8), intent(in)              :: dlnTdx, dlnndx, Delta_x, Delta_y, z_min, z_max
  real(kind=8), dimension(:), intent(in):: b_arr, dbdx_arr, dbdy_arr, sqrtg_arr, theta_arr
  integer, intent(in)                   :: lam_res, z_res
  real(kind=8), dimension(:),allocatable:: z_arr, lam_arr, w_psi_arr, w_alpha_arr, G_arr, AE_per_lam, AE_over_z
  integer                               :: lam_idx, z_idx, num_cross, num_wells
  real(kind=8)                          :: z_val, lam_val
  real(kind=8), intent(out)             :: AE_tot



  allocate(z_arr(z_res),AE_over_z(z_res),lam_arr(lam_res),AE_per_lam(lam_res))
  z_arr = (/(( z_min + (z_max-z_min)*( i )/z_res ), i=0,z_res-1)/)
  lam_arr =  (/(( lam_min + (lam_max-lam_min)*( i )/(lam_res+1) ), i=1,lam_res)/)

  do lam_idx = 1, lam_res, 1

    lam_val = lam_arr(lam_idx)
    call zero_cross_idx(1.0 - lam_val * b_arr,zero_idx,num_cross)
    num_wells = num_cross / 2
    allocate(w_psi_arr(num_wells),w_alpha_arr(num_wells),G_arr(num_wells))
    call w_bounce(b_arr,dbdx_arr,dbdy_arr,sqrtg_arr,theta_arr,L_tot,lam,Delta_x,Delta_y,num_wells,w_psi_arr,w_alpha_arr,G_arr)

    do z_idx = 1, z_res, 1
      z_val   = z_arr(z_idx)
      call AE_integrand(dlnTdx,dlnndx,Delta_x,z_val,w_alpha_arr,w_psi_arr,G_arr,AE)
      AE_over_z(z_idx) = AE
    end do

    deallocate(w_psi_arr,w_alpha_arr,G_arr)
    AE_per_lam(lam_idx) = trapz(z_arr, AE_over_z)

  end do

  AE_tot = trapz(lam_arr, AE_per_lam)
end subroutine AE_total