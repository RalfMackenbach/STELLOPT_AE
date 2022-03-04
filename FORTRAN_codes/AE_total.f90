! Calculate total AE
subroutine AE_total(q0,dlnTdx,dlnndx,Delta_x,Delta_y,b_arr,dbdx_arr,dbdy_arr, &
                    sqrtg_arr,theta_arr,lam_res,z_res,z_min,z_max,N,L_tot,AE_tot)
  implicit none
  integer, intent(in)                   :: N
  real(kind=8), intent(in)              :: q0, dlnTdx, dlnndx, Delta_x, Delta_y, z_min, z_max, L_tot
  real(kind=8), dimension(N), intent(in):: b_arr, dbdx_arr, dbdy_arr, sqrtg_arr, theta_arr
  integer, intent(in)                   :: lam_res, z_res
  real(kind=8), dimension(:),allocatable:: z_arr, lam_arr, w_psi_arr, w_alpha_arr, G_arr, AE_per_lam, AE_over_z
  integer                               :: lam_idx, z_idx, num_cross, num_wells, i
  real(kind=8), dimension(:,:), allocatable :: bounce_arr
  integer,      dimension(:,:), allocatable :: bounce_idx
  integer, dimension(:), allocatable    :: zero_idx
  real(kind=8)                          :: z_val, lam_val, AE, lam, lam_max, lam_min, AE_int
  real(kind=8), intent(out)             :: AE_tot


  ! Allocate arrays
  allocate(z_arr(z_res),AE_over_z(z_res),lam_arr(lam_res),AE_per_lam(lam_res))

  ! Find min and max values of lambda
  lam_min = 1/(MAXVAL(b_arr))
  lam_max = 1/(MINVAL(b_arr))

  ! Make arrays for lambda and normalized energy
  z_arr =    (/(( z_min + (z_max-z_min)*( i )/z_res ), i=0,z_res-1)/)
  lam_arr =  (/(( lam_min + (lam_max-lam_min)*( i )/(lam_res+1) ), i=1,lam_res)/)

  ! Loop over lambda indices
  do lam_idx = 1, lam_res, 1

    ! If bounce array is already allocated from previous loop,
    ! deallocate
    IF (ALLOCATED(bounce_arr)) THEN
      DEALLOCATE(bounce_arr,bounce_idx)
    END IF

    ! Assign lambda value
    lam_val = lam_arr(lam_idx)

    ! Count the number of zero crossings
    call zero_cross_count(1.0 - lam_val*b_arr,N,num_cross)
    ! Allocate the bounce arrays
    allocate(bounce_arr(num_cross/2,2),bounce_idx(num_cross/2,2))
    ! Find the bounce well indices and crossings
    call bounce_wells(b_arr,theta_arr,lam_val,N,num_cross,bounce_arr,bounce_idx)


    ! Number of wells should be number of crossings divides by 2
    num_wells = num_cross / 2
    ! Allocate the drifts arrays
    allocate(w_psi_arr(num_wells),w_alpha_arr(num_wells),G_arr(num_wells))
    ! Make drift arrays
    call w_bounce(q0,L_tot,b_arr,dbdx_arr,dbdy_arr,sqrtg_arr,theta_arr,lam_val, &
                  Delta_x,Delta_y,num_wells,w_psi_arr,w_alpha_arr,G_arr, &
                  N,bounce_idx,bounce_arr)
    ! print *,w_psi_arr,w_alpha_arr,G_arr
    ! print *,'lam = ',lam_val,'w_alp = ',w_alpha_arr,'G = ',G_arr
    print *,lam_val,w_psi_arr,w_alpha_arr,G_arr

    ! Loop over indices
    do z_idx = 1, z_res, 1
      z_val   = z_arr(z_idx)
      call AE_integrand(num_wells,dlnTdx,dlnndx,Delta_x,z_val,w_alpha_arr,w_psi_arr,G_arr,AE)
      AE_over_z(z_idx) = AE
    end do


    deallocate(w_psi_arr,w_alpha_arr,G_arr)
    call trapz(z_arr, AE_over_z,size(z_arr),AE_int)
    AE_per_lam(lam_idx) = AE_int


  end do

  call trapz(lam_arr, AE_per_lam, size(lam_arr),AE_tot)
end subroutine AE_total
