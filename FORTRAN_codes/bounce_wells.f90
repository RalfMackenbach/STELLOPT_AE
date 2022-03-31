


! find bounce wells
subroutine bounce_wells(b_arr,theta_arr,lam_val,N,num_cross,bounce_arr,bounce_idx)
  ! Constructs an array of bounce wells.  Structure is as follows
  ! bounce_arr [x_idx_l, x_idx_r, x_val_l, x_val_r]. Each row corresponds
  ! to a unique bounce well.
  ! Works with PERIODIC arrays (which have been fed through make periodic
  ! function).
  integer, intent(in)                                     :: N, num_cross
  real(kind=8), intent(in)                                :: lam_val
  real(kind=8), dimension(N), intent(in)                  :: b_arr, theta_arr
  real(kind=8), dimension(N)                              :: zero_arr, roll_arr
  logical, dimension(N)                                   :: mask
  integer, dimension(num_cross)                           :: zero_idx
  real(kind=8), dimension(num_cross/2,2), intent(out)     :: bounce_arr
  integer,      dimension(num_cross/2,2), intent(out)     :: bounce_idx
  integer                                                 :: num_wells
  integer :: first_well_end, do_idx, l_idx, r_idx

  ! define function of which we need zero crossings, and retrieve crossings
  zero_arr = 1.0 - lam_val * b_arr



  CALL zero_cross_idx(zero_arr,num_cross,N,zero_idx)
  ! print *,''
  ! print *,'zero indices are: ',zero_idx
  ! print *,''


  
  ! Check if even number of wells
  IF (MOD(num_cross,2).NE.0) THEN
    print *,'ERROR: Odd number of well crossings, please adjust lambda resolution'
  END IF
  num_wells = num_cross/2


  ! Check if the first crossing is end of well
  IF  ( b_arr(zero_idx(1)+1) - b_arr(zero_idx(1)) .GT. 0 ) THEN
    first_well_end = 1
  ELSE
    first_well_end = 0
  END IF




  ! First let's fill up the bounce_idx array
  ! If well crosses periodicity we must shift the indices
  IF (first_well_end .EQ. 1) THEN
    zero_idx = CSHIFT(zero_idx,-1)
  END IF

  ! Fill up bounce array
  DO  do_idx = 1, num_wells,1
      l_idx                     = zero_idx(2*do_idx-1)
      r_idx                     = zero_idx(2*do_idx)
      bounce_idx(do_idx  ,1)    = l_idx
      bounce_idx(do_idx  ,2)    = r_idx
      bounce_arr(do_idx, 1)     = (-(zero_arr(l_idx+1)*theta_arr(l_idx)) + &
                                  zero_arr(l_idx)*theta_arr(l_idx+1))/(zero_arr(l_idx) - &
                                  zero_arr(l_idx+1))
      bounce_arr(do_idx, 2)     = (-(zero_arr(r_idx+1)*theta_arr(r_idx)) + &
                                  zero_arr(r_idx)*theta_arr(r_idx+1))/(zero_arr(r_idx) - &
                                  zero_arr(r_idx+1))
  END DO



end subroutine bounce_wells
