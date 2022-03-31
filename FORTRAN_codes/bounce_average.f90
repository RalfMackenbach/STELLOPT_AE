
! calculate bounce average
subroutine bounce_average(h_arr,b_arr,theta_arr,lam,num_wells,bounce_idx, &
                          bounce_arr,bounce_ave,N)
  ! Calculates the bounce averaged precession frequencies for all bounce wells
  ! of a given lambda.
  implicit none
  integer,      intent(in)                            :: N, num_wells
  real(kind=8), dimension(N), intent(in)              :: h_arr, b_arr, theta_arr
  integer,      dimension(num_wells,2), intent(in)    :: bounce_idx
  real(kind=8), dimension(num_wells,2), intent(in)    :: bounce_arr
  real(kind=8), dimension(N)                          :: f_arr
  real(kind=8), intent(in)                            :: lam
  integer                                             :: do_idx, N_alloc
  integer,      dimension(num_wells)                  :: l_idx, r_idx
  real(kind=8), dimension(num_wells)                  :: l_cross, r_cross
  real(kind=8), dimension(num_wells), intent(out)     :: bounce_ave
  real(kind=8), dimension(:), allocatable             :: h_i, h_j, f_i, f_j, theta_i, theta_j
  real(kind=8)  :: inner, left, right, h_l, h_r, y_l, y_r
  ! If allocate, deallocate


  l_idx     = reshape(bounce_idx(:,1), shape(l_idx))
  r_idx     = reshape(bounce_idx(:,2), shape(r_idx))
  l_cross   = reshape(bounce_arr(:,1), shape(l_cross))
  r_cross   = reshape(bounce_arr(:,2), shape(r_cross))

  ! make f_arr
  f_arr = 1.0 - lam*b_arr


  DO  do_idx = 1, num_wells, 1
    ! If allocated, deallocate
    IF (ALLOCATED(h_i)) THEN
      DEALLOCATE(h_i, h_j, f_i, f_j, theta_i, theta_j)
    END IF

    ! Check if inner int crosses periodicity boundary
    IF (l_idx(do_idx).GT.r_idx(do_idx)) THEN
      ! print *,'well crosses periodicity boundary'
      ! Split up inner int into two parts
      ! first left-to-end
      N_alloc = SIZE(h_arr) - l_idx(do_idx) -1
      allocate( h_i(N_alloc), h_j(N_alloc), f_i(N_alloc), f_j(N_alloc), &
                theta_i(N_alloc), theta_j(N_alloc) )
      h_i     = h_arr((l_idx(do_idx) + 1):(SIZE(h_arr)-1))
      h_j     = h_arr((l_idx(do_idx) + 2):(SIZE(h_arr)  ))
      f_i     = f_arr((l_idx(do_idx) + 1):(SIZE(f_arr)-1))
      f_j     = f_arr((l_idx(do_idx) + 2):(SIZE(f_arr)  ))
      theta_i = theta_arr((l_idx(do_idx) + 1):(SIZE(theta_arr)-1))
      theta_j = theta_arr((l_idx(do_idx) + 2):(SIZE(theta_arr)  ))
      CALL inner_bounce_trapz(N_alloc,h_i,h_j,f_i,f_j,theta_i,theta_j,y_l)
      ! then start-to-right
      N_alloc = size(h_arr((1):(r_idx(do_idx)-1)))
      deallocate(h_i, h_j, f_i, f_j, theta_i, theta_j)
      allocate( h_i(N_alloc), h_j(N_alloc), f_i(N_alloc), f_j(N_alloc), &
                theta_i(N_alloc), theta_j(N_alloc) )
      h_i     = h_arr((1):(r_idx(do_idx)-1))
      h_j     = h_arr((2):(r_idx(do_idx)  ))
      f_i     = f_arr((1):(r_idx(do_idx)-1))
      f_j     = f_arr((2):(r_idx(do_idx)  ))
      theta_i = theta_arr((1):(r_idx(do_idx)-1))
      theta_j = theta_arr((2):(r_idx(do_idx)  ))
      CALL inner_bounce_trapz(N_alloc,h_i,h_j,f_i,f_j,theta_i,theta_j,y_r)
      inner = y_l + y_r



    ! otherwise business as usual
    ELSE
      ! print *,'well does not cross periodicity boundary'
      N_alloc = r_idx(do_idx) - l_idx(do_idx) - 1
      allocate( h_i(N_alloc), h_j(N_alloc),f_i(N_alloc), f_j(N_alloc), &
                theta_i(N_alloc), theta_j(N_alloc) )
      h_i     = h_arr((l_idx(do_idx) + 1):(r_idx(do_idx)-1))
      h_j     = h_arr((l_idx(do_idx) + 2):(r_idx(do_idx)  ))
      f_i     = f_arr((l_idx(do_idx) + 1):(r_idx(do_idx)-1))
      f_j     = f_arr((l_idx(do_idx) + 2):(r_idx(do_idx)  ))
      theta_i = theta_arr((l_idx(do_idx) + 1):(r_idx(do_idx)-1))
      theta_j = theta_arr((l_idx(do_idx) + 2):(r_idx(do_idx)  ))


      CALL inner_bounce_trapz(N_alloc,h_i,h_j,f_i,f_j,theta_i,theta_j,inner)
    END IF

    ! Now do the edge integrals
    h_l = h_arr(l_idx(do_idx)) + (l_cross(do_idx) - &
          theta_arr(l_idx(do_idx)))/(theta_arr(l_idx(do_idx)+1) - &
          theta_arr(l_idx(do_idx) )) * ( h_arr(l_idx(do_idx)+1) - &
          h_arr(l_idx(do_idx)) )
    CALL left_bounce_trapz(h_l,h_arr(l_idx(do_idx)+1),0.0, &
                           f_arr(l_idx(do_idx)+1),l_cross(do_idx), &
                           theta_arr(l_idx(do_idx)+1),left)


    h_r = h_arr(r_idx(do_idx)) + (r_cross(do_idx) - &
          theta_arr(r_idx(do_idx)))/(theta_arr(r_idx(do_idx)+1) - &
          theta_arr(r_idx(do_idx) )) * ( h_arr(r_idx(do_idx)+1) - &
          h_arr(r_idx(do_idx)) )
    CALL right_bounce_trapz(h_arr(r_idx(do_idx)),h_r,f_arr(r_idx(do_idx)), &
                            0.0,theta_arr(r_idx(do_idx)),r_cross(do_idx), &
                            right)

    ! finally, fill in full integral!
    ! print *,'left, inner, right =',left,inner,right
    bounce_ave(do_idx)= left + inner + right

  END DO
end subroutine bounce_average
