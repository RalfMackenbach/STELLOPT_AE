! Author: r.j.j.mackenbach@tue.nl
! This code computes bounce averaged precession frequencies, and
! the available energy (AE) of trapped particles. It assumes all relevant arrays
! have been supplied (i.e. reading GIST should be done seperately)

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
  print *,Bhat
  CLOSE(iunit)


  AE_tot = 1.0
end subroutine compute_AE

program debug
  real(kind=8) :: AE_tot
  CALL compute_AE(AE_tot)
end program debug
