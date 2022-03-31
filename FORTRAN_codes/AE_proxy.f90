! Author: r.j.j.mackenbach@tue.nl
! This code computes bounce averaged precession frequencies, and
! the available energy (AE) of trapped particles. It assumes all relevant arrays
! have been supplied (i.e. reading GIST should be done seperately)

!! READ FILES AND COMPUTE AE !!


! subroutine compute_AE(AE_tot)
!   real(kind=8), dimension(:,:),allocatable  :: g11, g12, g22, Bhat, abs_jac, L2, L1, dBdt, dummy
!   real(kind=8), dimension(:),  allocatable  :: Bhat_array, sqrtg_arr, L2_array, L1_array
!
!   integer                                   :: gridpoints, n_pol, iunit, i, ialpha, N
!   real(kind=8)                              :: my_dpdx, q0, shat, Delta_x, Delta_y
!   real(kind=8), dimension(:), allocatable   :: y
!   real(kind=8)                              :: z_min, z_max, dlnndx, dlnTdx
!   real(kind=8)                              :: a, pi, B_ave, r1, r2, L_tot, omn, omt
!   real(kind=8), intent(out)                 :: AE_tot
!
!   ! Import namelist
!   namelist /PARAMETERS/ my_dpdx, q0, shat, gridpoints, n_pol
!
!   ! Import arrays
!   iunit  = 1
!   ialpha = 1
!
!   OPEN(UNIT=iunit, FILE='gist_d3d_vac.txt')
!   !!!!! READING DATA !!!!!
!   READ(iunit,NML=PARAMETERS)
!   allocate(g11(1,gridpoints), g12(1,gridpoints), g22(1,gridpoints),  &
!            Bhat(1,gridpoints), abs_jac(1,gridpoints), L2(1,gridpoints), &
!            L1(1,gridpoints), dBdt(1,gridpoints), dummy(1,gridpoints))
!   DO i = 1, gridpoints
!      READ(iunit,"(9ES20.10)") g11(ialpha,i),g12(ialpha,i),g22(ialpha,i),&
!                               Bhat(ialpha,i),abs_jac(ialpha,i),L2(ialpha,i),&
!                               L1(ialpha,i),dBdt(ialpha,i), dummy(ialpha,i)
!   END DO
!
!   ! define scalars
!   pi= 4.D0*DATAN(1.D0)
!   a = n_pol*2*pi ! max of theta array
!   N = size(Bhat)
!
!   ! allocate data
!   allocate(y(N),Bhat_array(N),L1_array(N),L2_array(N),sqrtg_arr(N))
!   ! Reshape arrays
!   Bhat_array    = reshape(Bhat(1,:),    shape(Bhat_array))
!   L1_array      = reshape(L1(1,:),      shape(Bhat_array))
!   L2_array      = reshape(L2(1,:),      shape(Bhat_array))
!   sqrtg_arr = 1.0/reshape(abs_jac(1,:), shape(Bhat_array))
!   y = (/((i*a)/(N-1),i=0,N-1)/) ! assign to y an array for theta [0,a]
!   ! print *,y
!
!
!   ! Calculate average B (arclength-averaged)
!   call trapz(y,q0*Bhat_array*Bhat_array*sqrtg_arr,N,r1)
!   call trapz(y,q0*Bhat_array*sqrtg_arr,N,r2)
!   B_ave = r1/r2
!   Delta_x = q0
!   Delta_y = q0
!   z_min = 0.0001
!   z_max = 40.0
!   omn   = 1.0
!   omt   = 0.0
!   dlnndx = -omn
!   dlnTdx = -omt
!
!   call trapz(y,q0*Bhat_array*sqrtg_arr/B_ave,N,L_tot)
!
!
!   CALL AE_total(q0,dlnTdx,dlnndx,Delta_x,Delta_y,Bhat_array/B_ave, &
!                 L1_array/B_ave,L2_array/B_ave,sqrtg_arr,y,10000,10000, &
!                 z_min,z_max,N,L_tot,AE_tot)
!
!   CLOSE(iunit)
! end subroutine compute_AE


! subroutine AE_benchmark(N_theta,AE_tot)
!   integer , intent(in)                      :: N_theta
!   real(kind=8), dimension(N_theta)          :: theta, mod_b, dbdx, dbdy, sqrtg
!   real(kind=8)                              :: q0, z_min, z_max, dlnndx, dlnTdx, &
!                                                Delta_x, Delta_y, theta_b, &
!                                                theta_x, theta_y, one, L_tot
!   integer                                   :: i, shift_int
!   real(kind=8), intent(out)                 :: AE_tot
!
!
!
!   ! define scalars
!   q0     = 1D0
!   z_min  = 1D-5
!   z_max  = 1D1
!   dlnndx = -4D0
!   dlnTdx = 0D0
!   Delta_x= q0
!   Delta_y= q0
!   theta_b= 1D0
!   theta_x= 1D0
!   theta_y= 1D0
!   one    = 1D0
!   L_tot  = 1D0
!
!
!   shift_int = 0
!   ! Make theta_arr
!   theta = (2.0)*(/(i, i=0,N_theta-1, 1)/)/(N_theta-1) - one ! assign to y an array for theta [-1,1]
!   mod_b = CSHIFT(+one + (theta/theta_b)**2.0,                     shift_int)
!   dbdx  = CSHIFT(-one + (theta/theta_x)**2.0,                     shift_int)
!   dbdy  = CSHIFT(       (theta/theta_y),                          shift_int)
!   sqrtg = CSHIFT((/(one,i=0,N_theta-1)/),                         shift_int)
!
!   ! (q0, dlnTdx, dlnndx, Delta_x, Delta_y, b_arr, dbdx_arr, dbdy_arr, sqrtg_arr, theta_arr, lam_res, z_res, z_min, z_max, N, AE_tot)
!   call AE_total(q0,dlnTdx,dlnndx,Delta_x,Delta_y,mod_b,dbdx,dbdy, &
!                 sqrtg,theta,10,100,z_min,z_max,1D-10, &
!                 N_theta,L_tot,AE_tot)
! end subroutine AE_benchmark


! program debug
!  real(kind=8) :: AE_tot
!  CALL AE_benchmark(10000,AE_tot)
!  ! CALL compute_AE(AE_tot)
!  print *,AE_tot
! end program debug

program debug_SAL
  USE trapped_avail_energy_mod
  real(kind=8) :: AE_tot
  CHARACTER(256) :: filename
  filename='gist_files/gist_STELLOPT.txt'
  CALL compute_AE_GIST(filename,AE_tot)
  print *,AE_tot

  filename='gist_files/gist_HSX.txt'
  CALL compute_AE_GIST(filename,AE_tot)
  print *,AE_tot

  filename='gist_files/gist_NCSX.txt'
  CALL compute_AE_GIST(filename,AE_tot)
  print *,AE_tot

  filename='gist_files/gist_W7XHM.txt'
  CALL compute_AE_GIST(filename,AE_tot)
  print *,AE_tot

  filename='gist_files/gist_W7XSC.txt'
  CALL compute_AE_GIST(filename,AE_tot)
  print *,AE_tot
end program debug_SAL
