


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