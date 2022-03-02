


! find bounce wells
subroutine bounce_wells(b_arr,theta_arr,lam_val,N,num_cross,bounce_arr,bounce_idx)
  ! Constructs an array of bounce wells.  Structure is as follows
  ! bounce_arr [x_idx_l, x_idx_r, x_val_l, x_val_r]. Each row corresponds
  ! to a unique bounce well.
  integer, intent(in)                                     :: N, num_cross
  real(kind=8), intent(in)                                :: lam_val
  real(kind=8), dimension(N), intent(in)                  :: b_arr, theta_arr
  real(kind=8), dimension(N)                              :: zero_arr, roll_arr
  logical, dimension(N)                                   :: mask
  integer, dimension(num_cross)                           :: zero_idx, zero_idx_shift
  real(kind=8), dimension(num_cross/2,2), intent(out)     :: bounce_arr
  integer,      dimension(num_cross/2,2), intent(out)     :: bounce_idx
  integer                                                 :: num_wells
  integer :: first_well_end, well_disc, do_idx, l_idx, r_idx

  ! define function of which we need zero crossings, and retrieve crossings
  zero_arr = 1.0 - lam_val * b_arr

  CALL zero_cross_idx(zero_arr,num_cross,N,zero_idx)
  ! print *,''
  ! print *,'zero indices are: ',zero_idx
  ! print *,''
  ! Check if even number of wells
  IF (MOD(num_cross,2).NE.0) THEN
    print *,'Odd number of well crossings, please adjust lambda resolution'
  END IF
  num_wells = num_cross/2


  ! Check if the first crossing is end of well
  IF  ( b_arr(zero_idx(1)+1) - b_arr(zero_idx(1)) .GT. 0 ) THEN
    first_well_end = 1
  ELSE
    first_well_end = 0
  END IF
  ! print *,''
  ! print *,'first_well_end=',first_well_end
  ! print *,''
  ! If it is an end, there are two options for the type of well it is:
  ! 1) the well crosses the periodicity boundary.
  ! 2) the well ends discontinously at the periodicity boundary
  ! Let's if a zero crossing is at the discontinous periodicity boundary
  IF ( zero_idx(SIZE(zero_idx)) .EQ. SIZE(zero_arr) ) THEN
    well_disc = 1
  ELSE
    well_disc = 0
  END IF
  ! print *,''
  ! print *,'well_disc=',well_disc
  ! print *,''




  ! First let's fill up the bounce_idx array
  ! If well crosses periodicity we must shift the indices
  IF (well_disc .EQ. 0 .AND. first_well_end .EQ. 1) THEN
    zero_idx_shift = CSHIFT(zero_idx,-1)
    DO  do_idx = 1, num_wells,1
      l_idx                     = zero_idx_shift(2*do_idx-1)
      r_idx                     = zero_idx_shift(2*do_idx)
      bounce_idx(do_idx  ,1)    = l_idx
      bounce_idx(do_idx  ,2)    = r_idx
      bounce_arr(do_idx, 1)     = (-(zero_arr(l_idx+1)*theta_arr(l_idx)) + &
                                  zero_arr(l_idx)*theta_arr(l_idx+1))/(zero_arr(l_idx) - &
                                  zero_arr(l_idx+1))
      bounce_arr(do_idx, 2)     = (-(zero_arr(r_idx+1)*theta_arr(r_idx)) + &
                                  zero_arr(r_idx)*theta_arr(r_idx+1))/(zero_arr(r_idx) - &
                                  zero_arr(r_idx+1))
    END DO
  END IF





  ! If first index is end of a well and the discontinuity is a well crossing,
  ! then the well is defined by theta=0 to first crossing.
  IF (well_disc .EQ. 1 .AND. first_well_end .EQ. 1) THEN
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
    bounce_idx(1  ,1) = 0 ! left of first index is zero
    bounce_arr(1, 1)  = theta_arr(1) ! first theta val is bounce point
  END IF







  ! If first index is not end of well and the discontinuity is a well crossing,
  ! then last index is end of well
  IF (well_disc .EQ. 1 .AND. first_well_end .EQ. 0) THEN
    DO  do_idx = 1, num_wells-1, 1
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
    l_idx                     = zero_idx(2*do_idx-1)
    bounce_idx(num_wells, 1) = l_idx
    bounce_arr(num_wells, 1)  = (-(zero_arr(l_idx+1)*theta_arr(l_idx)) + &
                                zero_arr(l_idx)*theta_arr(l_idx+1))/(zero_arr(l_idx) - &
                                zero_arr(l_idx+1))
    bounce_idx(num_wells, 2)  = size(theta_arr) ! left of first index is zero
    bounce_arr(num_wells, 2)  = theta_arr(size(theta_arr)) ! first theta val is bounce point
  END IF





  ! Otherewise, it's business as usual
  IF (well_disc .EQ. 0 .AND. first_well_end .EQ. 0) THEN
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
  END IF



end subroutine bounce_wells
