! Author: r.j.j.mackenbach@tue.nl
! This code computes bounce averaged precession frequencies, and
! the available energy (AE) of trapped particles. It assumes all relevant arrays
! have been supplied (i.e. reading GIST should be done seperately)


! just regular trapz
pure function trapz(x, y) result(r)
   !! Found this little routine online
   !! Calculates the integral of an array y with respect to x using the trapezoid
   !! approximation. Note that the mesh spacing of x does not have to be uniform.
   real(kind=8), intent(in)  :: x(:)         !! Variable x
   real(kind=8), intent(in)  :: y(size(x))   !! Function y(x)
   real(kind=8)              :: r            !! Integral ∫y(x)·dx

   ! Integrate using the trapezoidal rule
   associate(n => size(x))
     r = SUM((y(1+1:n-0) + y(1+0:n-1))*(x(1+1:n-0) - x(1+0:n-1)))/2
   end associate
end function


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


! trapz-like estimate of bounce-ave function, left edge
subroutine left_bounce_trapz(h_l,h_0,f_l,f_0,theta_l,theta_0,y)
  ! Estimation of integral of edge integral h/sqrt(f) dtheta, from theta_l where
  ! f = 0, to the first theta node to its right theta_0
  implicit none
  real(kind=8), dimension(:), intent(in) :: h_l
  real(kind=8), dimension(:), intent(in) :: h_0
  real(kind=8), dimension(:), intent(in) :: f_l
  real(kind=8), dimension(:), intent(in) :: f_0
  real(kind=8), dimension(:), intent(in) :: theta_l
  real(kind=8), dimension(:), intent(in) :: theta_0
  real(kind=8), dimension(:), intent(out) :: y
  y = (2*(h_0 - 4*h_l)*(-theta_0 + theta_l))/(3.*Sqrt(f_0))
end subroutine left_bounce_trapz


! trapz-like estimate of bounce-ave function, right edge
subroutine right_bounce_trapz(h_n,h_r,f_n,f_r,theta_n,theta_r,y)
  ! Estimation of integral of edge integral h/sqrt(f) dtheta, from theta_n to theta_r
  ! where f = 0
  implicit none
  real(kind=8), dimension(:), intent(in) :: h_n
  real(kind=8), dimension(:), intent(in) :: h_r
  real(kind=8), dimension(:), intent(in) :: f_n
  real(kind=8), dimension(:), intent(in) :: f_r
  real(kind=8), dimension(:), intent(in) :: theta_n
  real(kind=8), dimension(:), intent(in) :: theta_r
  real(kind=8), dimension(:), intent(out):: y
  y = (2*(h_n + 2*h_r)*(-theta_n + theta_r))/(3.*Sqrt(f_n))
end subroutine right_bounce_trapz


! find zero crossings
subroutine zero_cross_idx(y_arr,zero_idx,num_cross)
  ! Finds all zeros of y. Based on "Numerical Recipes in Fortran 90" (1996),
  ! chapter B9, p.1184.
  implicit none
  real(kind=8), dimension(:), intent(in)  :: y_arr
  real :: lam
  real(kind=8), dimension(size(y_arr))    :: roll_arr
  logical, dimension(size(y_arr))         :: mask
  integer, intent(out) :: num_cross
  real, dimension (:), allocatable, intent(out) :: zero_idx
  integer :: ix

  ! Make shifted array (assumes periodic boundary conditions)
  roll_arr  = CSHIFT(y_arr,1)

  ! check zero crossings, returns TRUE left of the crossing. Also, count number
  ! of crossings
  mask      = y_arr*roll_arr <= 0
  num_cross = COUNT(mask)

  ! Allocate arrays
  allocate ( zero_idx(num_cross) )

  ! Return indices where mask is true, so left of crossings
  zero_idx = PACK([(ix,ix=1,size(mask))],mask)
end subroutine zero_cross_idx


! find bounce wells
subroutine bounce_wells(b_arr,theta_arr,lam,bounce_arr,bounce_idx,num_wells)
  ! Constructs an array of bounce wells.  Structure is as follows
  ! bounce_arr [x_idx_l, x_idx_r, x_val_l, x_val_r]. Each row corresponds
  ! to a unique bounce well.
  real(kind=8), dimension(:), intent(in)                  :: b_arr, theta_arr
  real(kind=8), dimension(size(b_arr))                    :: zero_arr, roll_arr
  logical, dimension(size(b_arr))                         :: mask
  real(kind=8), dimension(:), allocatable                 :: zero_idx
  real(kind=8), dimension(:,:), allocatable, intent(out)  :: bounce_arr
  integer,      dimension(:,:), allocatable, intent(out)  :: bounce_idx
  integer, intent(out)                                    :: num_wells
  integer :: first_well_end, well_disc, num_cross, do_idx, l_idx, r_idx

  ! define function of which we need zero crossings, and retrieve crossings
  zero_arr = 1.0 - lam * b_arr
  roll_arr  = CSHIFT(zero_arr,1)
  mask      = zero_arr*roll_arr <= 0
  num_cross = COUNT(mask)
  allocate(zero_idx(num_cross))
  CALL zero_cross_idx(zero_arr,zero_idx,num_cross)

  ! Check if even number of wells
  IF (MOD(num_cross,2).NE.0) THEN
    ! print('Odd number of wells crossings, brace for crash...')
  END IF
  num_wells = num_cross/2

  ! Allocate arrays
  allocate ( bounce_arr(num_wells,2) )
  allocate ( bounce_idx(num_wells,2) )

  ! Check if the first crossing is end of well
  IF  ( zero_arr(zero_idx(1)+1) - zero_arr(zero_idx(1)) < 0 ) THEN
    first_well_end = 1
  ELSE
    first_well_end = 0
  END IF

  ! If it is an end, there are two options for the type of well it is:
  ! 1) the well crosses the periodicity boundary.
  ! 2) the well ends discontinously at the periodicity boundary
  ! Let's if a zero crossing is at the discontinous periodicity boundary
  IF ( zero_idx(SIZE(zero_idx)) == SIZE(zero_arr) ) THEN
    well_disc = 1
  ELSE
    well_disc = 0
  END IF

  ! If the first crossing is end AND the discontinuity is a zero crossing, then
  ! the first well is defined by the range between this discontinuity and the
  ! first crossing
  IF (well_disc == 1 .AND. first_well_end == 1) THEN
    r_idx            = zero_idx(1)
    bounce_idx(1, 1) = 0 ! 0 because indexing starts at 1 so left of it = 0
    bounce_idx(1, 2) = r_idx
    bounce_arr(1, 1) = theta_arr(1)
    ! linear interpolation to find theta_val
    bounce_arr(1, 2) = (-(zero_arr(r_idx+1)*theta_arr(r_idx)) + &
                        zero_arr(r_idx)*theta_arr(r_idx+1))/(zero_arr(r_idx) - &
                        zero_arr(r_idx+1))

  ! If the first crossing is end AND the discontinuity is NOT a zero crossing,
  ! then the first well is defined by the last zero crossing to the first.
  ELSE IF (well_disc == 0 .AND. first_well_end == 1) THEN
    l_idx            = zero_idx(SIZE(zero_idx))
    r_idx            = zero_idx(1)
    bounce_idx(1, 1) = l_idx
    bounce_idx(1, 2) = r_idx
    ! linear interpolation to find theta_val
    bounce_arr(1, 1) = (-(zero_arr(l_idx+1)*theta_arr(l_idx)) + &
                        zero_arr(l_idx)*theta_arr(l_idx+1))/(zero_arr(l_idx) - &
                        zero_arr(l_idx+1))
    bounce_arr(1, 2) = (-(zero_arr(r_idx+1)*theta_arr(r_idx)) + &
                        zero_arr(r_idx)*theta_arr(r_idx+1))/(zero_arr(r_idx) - &
                        zero_arr(r_idx+1))
  END IF

  ! Fill up bounce_well array
  DO  do_idx = 1 + first_well_end, num_wells, 1
    l_idx                 = zero_idx(2*(do_idx) -1 + first_well_end)
    r_idx                 = zero_idx(2*(do_idx)    + first_well_end)
    bounce_idx(do_idx, 1) = l_idx
    bounce_idx(do_idx, 2) = r_idx
    ! linear interpolation to find theta_val
    bounce_arr(do_idx, 1) = (-(zero_arr(l_idx+1)*theta_arr(l_idx)) + &
                            zero_arr(l_idx)*theta_arr(l_idx+1))/ &
                            (zero_arr(l_idx) - zero_arr(l_idx+1))
    bounce_arr(do_idx, 2) = (-(zero_arr(r_idx+1)*theta_arr(r_idx)) + &
                            zero_arr(r_idx)*theta_arr(r_idx+1))/ &
                            (zero_arr(r_idx) - zero_arr(r_idx+1))
  END DO

  ! Finally, if there is a discont. crossing, but first crossing is start of a
  ! well, then the last crossing point is at last theta val
  IF (well_disc == 1 .AND. first_well_end == 0) THEN
    bounce_arr(num_wells, 2) = theta_arr(SIZE(theta_arr))
  END IF
end subroutine bounce_wells


! calculate bounce average
subroutine bounce_average(h_arr,b_arr,theta_arr,lam,num_wells,bounce_idx, &
                          bounce_arr,bounce_ave)
  ! Calculates the bounce averaged precession frequencies for all bounce wells
  ! of a given lambda.
  implicit none
  real(kind=8), dimension(:), intent(in)              :: h_arr, b_arr, theta_arr
  integer,      dimension(:), intent(in)              :: bounce_idx
  real(kind=8), dimension(:), intent(in)              :: bounce_arr
  real(kind=8), dimension(size(b_arr))                :: f_arr
  real(kind=8), intent(in)                            :: lam
  integer,      intent(in)                            :: num_wells
  integer                                             :: do_idx, first_well_end
  integer,      dimension(:), allocatable             :: l_idx, r_idx
  real(kind=8), dimension(:), allocatable             :: l_cross, r_cross
  real(kind=8), dimension(:), allocatable, intent(out):: bounce_ave
  real(kind=8)  :: h_i, h_j, f_i, f_j, theta_i, theta_j, inner, left, right, &
                   h_l, h_r, y_l, y_r
  ! Allocate
  allocate ( bounce_ave(num_wells), l_idx(num_wells), r_idx(num_wells), &
            l_cross(num_wells), r_cross(num_wells) )

  l_idx     = bounce_idx(:,1)
  r_idx     = bounce_idx(:,2)
  l_cross   = bounce_arr(:,1)
  r_cross   = bounce_arr(:,2)

  ! Allocate
  allocate ( bounce_ave(num_wells) )

  ! make f_arr
  f_arr = 1.0 - lam*b_arr

  num_wells = SIZE(l_idx)
  DO  do_idx = 1 + first_well_end, num_wells, 1
    ! Check if inner int crosses periodicity boundary
    IF (l_idx > r_idx) THEN
      ! Split up inner int into two parts
      ! first left-to-end
      h_i     = h_arr((l_idx(do_idx) + 1):(SIZE(h_arr)-1))
      h_j     = h_arr((l_idx(do_idx) + 2):(SIZE(h_arr)  ))
      f_i     = f_arr((l_idx(do_idx) + 1):(SIZE(f_arr)-1))
      f_j     = h_arr((l_idx(do_idx) + 2):(SIZE(f_arr)  ))
      theta_i = theta_arr((l_idx(do_idx) + 1):(SIZE(theta_arr)-1))
      theta_j = theta_arr((l_idx(do_idx) + 2):(SIZE(theta_arr)  ))
      CALL inner_bounce_trapz(h_i,h_j,f_i,f_j,theta_i,theta_j,y_l)
      ! then start-to-right
      h_i     = h_arr((1):(r_idx(do_idx)-1))
      h_j     = h_arr((2):(r_idx(do_idx)))
      f_i     = f_arr((1):(r_idx(do_idx)-1))
      f_j     = h_arr((2):(r_idx(do_idx)))
      theta_i = theta_arr((1):(r_idx(do_idx)-1))
      theta_j = theta_arr((2):(r_idx(do_idx)))
      CALL inner_bounce_trapz(h_i,h_j,f_i,f_j,theta_i,theta_j,y_r)
      inner = y_l + y_r
    ! otherwise business as usual
    ELSE
      h_i     = h_arr((l_idx(do_idx) + 1):(r_idx(do_idx)-1))
      h_j     = h_arr((l_idx(do_idx) + 2):(r_idx(do_idx)  ))
      f_i     = f_arr((l_idx(do_idx) + 1):(r_idx(do_idx)-1))
      f_j     = f_arr((l_idx(do_idx) + 2):(r_idx(do_idx)  ))
      theta_i = theta_arr((l_idx(do_idx) + 1):(r_idx(do_idx)-1))
      theta_j = theta_arr((l_idx(do_idx) + 2):(r_idx(do_idx)  ))
      CALL inner_bounce_trapz(h_i,h_j,f_i,f_j,theta_i,theta_j,inner)
    END IF

    ! Now do the edge integrals

    ! check if left point is at the discont. boundary
    IF (l_idx(do_idx) == 0) THEN
      left = 0.0
    ! otherwise do the left int as usual
    ELSE
      h_l = h_arr(l_idx(do_idx)) + (l_cross(do_idx) - &
            theta_arr(l_idx(do_idx)))/(theta_arr(l_idx(do_idx)+1) - &
            theta_arr(l_idx(do_idx) )) * ( h_arr(l_idx(do_idx)+1) - &
            h_arr(l_idx(do_idx)) )
      CALL left_bounce_trapz(h_l,h_arr(l_idx(do_idx)+1),0.0, &
                             f_arr(l_idx(do_idx)+1),l_cross(do_idx), &
                             theta_arr(l_idx(do_idx)+1),left)
    END IF

    ! check if right point is at the discont. boundary
    IF (r_idx(do_idx) == SIZE(theta_arr)) THEN
      right = 0.0
    ! otherwise do the right int as usual
    ELSE
      h_r = h_arr(r_idx(do_idx)) + (r_cross(do_idx) - &
            theta_arr(r_idx(do_idx)))/(theta_arr(r_idx(do_idx)+1) - &
            theta_arr(r_idx(do_idx) )) * ( h_arr(r_idx(do_idx)+1) - &
            h_arr(r_idx(do_idx)) )
      CALL right_bounce_trapz(h_arr(r_idx(do_idx)),h_r,f_arr(r_idx(do_idx)), &
                              0.0,r_cross(do_idx),theta_arr(r_idx(do_idx)), &
                              right)
    END IF

    ! finally, fill in full integral!
    bounce_ave(do_idx)= left + inner + right
  END DO
end subroutine bounce_average


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


! integrand of AE
subroutine AE_integrand(dlnTdx,dlnndx,Delta_x,z,w_alpha,w_psi,G,AE)
  ! Given array of w_alpha, w_psi, and G (decided by lambda), and w_dia and z
  ! return the AE integrand at that (lambda,z).
  implicit none
  real(kind=8), dimension(:), intent(in)             :: w_alpha, w_psi, G
  real, intent(in)                                   :: dlnTdx, dlnndx, Delta_x, z
  real                                               :: w_dia
  real, intent(out)                                  :: AE
  ! problem with variable size again, unsure how this is handled.

  w_dia = Delta_x * ( dlnndx/z + dlnTdx * ( 1.0 - 3.0 / (2.0 * z) ) )
  AE = SUM((G*(w_psi*(-1 + Sqrt((-w_alpha + w_dia)**2 + w_psi**2)/ &
       Sqrt(w_alpha**2 + w_psi**2)) + w_alpha**2*(-1 + w_dia/w_alpha + &
       Sqrt((-w_alpha + w_dia)**2 + w_psi**2)/ &
       Sqrt(w_alpha**2 + w_psi**2)))*z**2.5)/E**z)
end subroutine AE_integrand


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





!! READ FILES AND COMPUTE AE !!

subroutine compute_AE(AE_tot)
  real(kind=8), dimension(:,:),allocatable:: g11, g12, g22, Bhat, abs_jac, L2, L1, dBdt, dummy
  integer                               :: gridpoints, n_pol, iunit, i, ialpha
  real(kind=8)                          :: my_dpdx, q0, shat
  real(kind=8), intent(out)             :: AE_tot
  namelist /PARAMETERS/ my_dpdx, q0, shat, gridpoints, n_pol


  iunit  = 1
  ialpha = 1
  OPEN(UNIT=iunit, FILE='gist_genet_wout_hsx.63e_pest_HSX_vac_s05_nz96.txt')
  !!!!! READING DATA !!!!!
  READ(iunit,NML=PARAMETERS)
  allocate(g11(1,gridpoints), g12(1,gridpoints), g22(1,gridpoints),  &
           Bhat(1,gridpoints), abs_jac(1,gridpoints), L2(1,gridpoints), &
           L1(1,gridpoints), dBdt(1,gridpoints), dummy(1,gridpoints))
  DO i = 1, gridpoints
     READ(iunit,"(9ES20.10)") g11(ialpha,i),g12(ialpha,i),g22(ialpha,i),&
                              Bhat(ialpha,i),abs_jac(ialpha,i),L2(ialpha,i),&
                              L1(ialpha,i),dBdt(ialpha,i), dummy(ialpha,i)
  END DO
  print(Bhat)
  CLOSE(iunit)


  AE_tot = 1.0
end subroutine compute_AE

program debug
  real(kind=8) :: AE_tot
  CALL compute_AE(AE_tot)
end program debug
