! Author: r.j.j.mackenbach@tue.nl
! This code computes bounce averaged precession frequencies, and
! the available energy (AE) of trapped particles. It assumes all relevant arrays
! have been supplied (i.e. reading GIST should be done seperately)

!! READ FILES AND COMPUTE AE !!


subroutine compute_AE(AE_tot)
  real(kind=8), dimension(:,:),allocatable  :: g11, g12, g22, Bhat, abs_jac, L2, L1, dBdt, dummy
  real(kind=8), dimension(:),  allocatable  :: Bhat_array, abs_jac_array, L2_array, L1_array

  integer                                   :: gridpoints, n_pol, iunit, i, ialpha, N
  real(kind=8)                              :: my_dpdx, q0, shat
  real(kind=8), dimension(:), allocatable   :: y
  real(kind=8)                              :: z_min, z_max, dlnndx, dlnTdx
  real(kind=8)                              :: a, pi, B_ave, r1, r2
  real(kind=8), intent(out)                 :: AE_tot

  ! Import namelist
  namelist /PARAMETERS/ my_dpdx, q0, shat, gridpoints, n_pol

  ! Import arrays
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

  ! define scalars
  pi= 4.D0*DATAN(1.D0)
  a = n_pol*2*pi ! max of theta array
  N = size(Bhat)

  ! allocate data
  allocate(y(N),Bhat_array(N),L1_array(N),L2_array(N),abs_jac_array(N))

  ! Reshape arrays
  Bhat_array    = reshape(Bhat(1,:),    shape(Bhat_array))
  L1_array      = reshape(L1(1,:),      shape(Bhat_array))
  L2_array      = reshape(L2(1,:),      shape(Bhat_array))
  abs_jac_array = reshape(abs_jac(1,:), shape(Bhat_array))
  y = (/((i*a)/(N-1),i=0,N-1)/) ! assign to y an array for theta [0,2*pi]


  ! Calculate average B (arclength-averaged)
  call trapz(y,q0*Bhat_array*Bhat_array*abs_jac_array,N,r1)
  call trapz(y,q0*Bhat_array*abs_jac_array,N,r2)
  B_ave = r1/r2

  z_min = 0.001
  z_max = 10.0
  dlnndx = 3.0
  dlnTdx = 0.0
  CALL AE_total(q0,dlnTdx,dlnndx,q0,q0,Bhat_array/B_ave,L1_array/B_ave,L2_array/B_ave, &
                abs_jac_array,y,1000,1000,z_min,z_max,N,AE_tot)

  CLOSE(iunit)
end subroutine compute_AE

program debug
  real(kind=8) :: AE_tot
  CALL compute_AE(AE_tot)
  print *,AE_tot
end program debug
