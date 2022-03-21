!-----------------------------------------------------------------------
!     Module:        trapped_avail_energy_mod
!     Authors:       R. Mackenbach (r.j.j.mackenbach@tue.nl)
!     Date:          03/21/2022
!     Description:   This module is used to compute the available energy
!                    of trapped particles three dimensional equilibira.
!
!-----------------------------------------------------------------------
      MODULE trapped_avail_energy_mod
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     Module Variables
!     PUBLIC
!        Y             This is for example
!                      
!     PRIVATE
!        X             This is for example
!-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(kind=8), PRIVATE :: X
      REAL(kind=8)          :: Y

!-----------------------------------------------------------------------
!     Private Subroutines
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     Public Subroutines
!-----------------------------------------------------------------------
      CONTAINS

      subroutine compute_AE_GIST(gist_file,AE_tot)
         real(kind=8), dimension(:,:),allocatable  :: g11, g12, g22, Bhat, abs_jac, L2, L1, dBdt, dummy
         real(kind=8), dimension(:),  allocatable  :: Bhat_array, sqrtg_arr, L2_array, L1_array

         integer                                   :: gridpoints, n_pol, iunit, i, ialpha, N
         real(kind=8)                              :: my_dpdx, q0, shat, Delta_x, Delta_y
         real(kind=8), dimension(:), allocatable   :: y
         real(kind=8)                              :: z_min, z_max, dlnndx, dlnTdx
         real(kind=8)                              :: a, pi, B_ave, r1, r2, L_tot, omn, omt
         real(kind=8), intent(out)                 :: AE_tot
         CHARACTER(*), INTENT(inout)               :: gist_file

         ! Import namelist
         namelist /PARAMETERS/ my_dpdx, q0, shat, gridpoints, n_pol

         ! Import arrays
         iunit  = 1
         ialpha = 1

         OPEN(UNIT=iunit, FILE=TRIM(gist_file))
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
         CLOSE(iunit)

         ! define scalars
         pi= 4.D0*DATAN(1.D0)
         a = n_pol*2*pi ! max of theta array
         N = size(Bhat)

         ! allocate data
         allocate(y(N),Bhat_array(N),L1_array(N),L2_array(N),sqrtg_arr(N))
         ! Reshape arrays
         Bhat_array    = reshape(Bhat(1,:),    shape(Bhat_array))
         L1_array      = reshape(L1(1,:),      shape(Bhat_array))
         L2_array      = reshape(L2(1,:),      shape(Bhat_array))
         sqrtg_arr = 1.0/reshape(abs_jac(1,:), shape(Bhat_array))
         y = (/((i*a)/(N-1),i=0,N-1)/) ! assign to y an array for theta [0,a]
         ! print *,y


         ! Calculate average B (arclength-averaged)
         call trapz(y,q0*Bhat_array*Bhat_array*sqrtg_arr,N,r1)
         call trapz(y,q0*Bhat_array*sqrtg_arr,N,r2)
         B_ave = r1/r2
         Delta_x = q0
         Delta_y = q0
         z_min = 0.0001
         z_max = 40.0
         omn   = 1.0
         omt   = 0.0
         dlnndx = -omn
         dlnTdx = -omt

         call trapz(y,q0*Bhat_array*sqrtg_arr/B_ave,N,L_tot)


         CALL AE_total(q0,dlnTdx,dlnndx,Delta_x,Delta_y,Bhat_array/B_ave, &
                       L1_array/B_ave,L2_array/B_ave,sqrtg_arr,y,10000,10000, &
                       z_min,z_max,N,L_tot,AE_tot)

       end subroutine compute_AE_GIST

      END MODULE trapped_avail_energy_mod
