
! calculate bounce average
subroutine bounce_average(h_arr,b_arr,theta_arr,lam,num_wells,bounce_idx, &
                          bounce_arr,bounce_ave)
  ! Calculates the bounce averaged precession frequencies for all bounce wells
  ! of a given lambda.
  implicit none
  real(kind=8), dimension(:), intent(in)              :: h_arr, b_arr, theta_arr
  integer,      dimension(:,:), intent(in)              :: bounce_idx
  real(kind=8), dimension(:,:), intent(in)              :: bounce_arr
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