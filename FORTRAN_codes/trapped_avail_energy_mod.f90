! Author: r.j.j.mackenbach@tue.nl
! This code computes bounce averaged precession frequencies, and
! the available energy (AE) of trapped particles. It assumes all relevant arrays
! have been supplied (i.e. reading GIST should be done seperately)


subroutine compute_AE(g11,g12,g22,Bhat,abs_jac,L1,L2,dBdt, &
                      gridpoints, n_pol, my_dpdx, q0, shat, omn, omt, &
                      z_min, z_max, z_res, lam_res, AE_tot)
  ! Some clarification on the inputs:
  ! g11,g12,g22,Bhat,abs_jac,L1,L2,dBdt are all arrays provided by GIST
  ! gridpoints, n_pol, my_dpdx, q0, shat are all scalars provided by GIST
  ! omn and omt are the density and temperature  gradient respectively, as defined in GENE
  ! z_min and z_max are minimal and maximal values of range of normalized energies
  ! z_res is the resolution of the integration range over normalized energies
  ! lam_res is the  resolution of the integration range over pitch angles.
  ! I typically use z_min = 1E-5, z_max = 40, z_res = lam_res = 10000.
  ! Please increase gridpoints for convergence checks.
  integer, intent(in)                       :: gridpoints, n_pol, z_res, lam_res
  real(kind=8), intent(in)                  :: my_dpdx, q0, shat, omn, omt
  real(kind=8), intent(in)                  :: z_min, z_max
  real(kind=8), dimension(N_arr), intent(in):: g11, g12, g22, Bhat, abs_jac, L2, L1, dBdt
  real(kind=8)                              :: Delta_x, Delta_y
  real(kind=8), dimension(N_arr)            :: y, sqrt_g
  real(kind=8)                              :: dlnndx, dlnTdx
  real(kind=8)                              :: a, pi, B_ave, r1, r2, L_tot
  integer                                   :: N
  real(kind=8), intent(out)                 :: AE_tot

  ! define scalars
  pi= 4.D0*DATAN(1.D0)
  a = n_pol*2*pi ! max of theta array
  N = size(Bhat)


  y = (/((i*a)/(N-1),i=0,N-1)/) ! assign to y an array for theta [0,a]
  sqrt_g = 1.0/abs_jac
  ! print *,y


  ! Calculate average B (arclength-averaged)
  call trapz(y,q0*Bhat_array*Bhat_array*sqrtg_arr,N,r1)
  call trapz(y,q0*Bhat_array*sqrtg_arr,N,r2)
  B_ave = r1/r2
  Delta_x = q0
  Delta_y = q0
  dlnndx = -omn
  dlnTdx = -omt

  call trapz(y,q0*Bhat_array*sqrtg_arr/B_ave,N,L_tot)


  CALL AE_total(q0,dlnTdx,dlnndx,Delta_x,Delta_y,Bhat_array/B_ave, &
                L1_array/B_ave,L2_array/B_ave,sqrtg_arr,y,lam_res,z_res, &
                z_min,z_max,N,L_tot,AE_tot)

end subroutine compute_AE
