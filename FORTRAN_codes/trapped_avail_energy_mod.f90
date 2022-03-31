!-----------------------------------------------------------------------
!     Module:        trapped_avail_energy_mod
!     Authors:       R. Mackenbach (r.j.j.mackenbach@tue.nl)
!     Date:          03/21/2022
!     Description:   This module is used to compute the available energy
!                    of trapped particles three dimensional equilibira.
!                    ( https://arxiv.org/pdf/2109.01042.pdf )
!
!-----------------------------------------------------------------------
      MODULE trapped_avail_energy_mod
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     Module Variables
!-----------------------------------------------------------------------
      IMPLICIT NONE

!-----------------------------------------------------------------------
!     Private Subroutines
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     Public Subroutines
!-----------------------------------------------------------------------
      CONTAINS

      SUBROUTINE compute_AE_GIST(gist_file,AE_tot)
         IMPLICIT NONE
      !-----------------------------------------------------------------
      !        Subroutine Input Variables
      !           gist_file :  Filename of gist file.
      !           AE_tot :     Total Free Energy.
      !-----------------------------------------------------------------
         real(kind=8), intent(out)                 :: AE_tot
         CHARACTER(*), INTENT(INOUT)               :: gist_file
      !-----------------------------------------------------------------
      !        Subroutine Variables
      !-----------------------------------------------------------------
         real(kind=8), dimension(:,:),allocatable  :: g11, g12, g22, Bhat, abs_jac, L2, L1, dBdt, dummy
         real(kind=8), dimension(:),  allocatable  :: Bhat_array, sqrtg_arr, L2_array, L1_array

         integer                                   :: gridpoints, n_pol, iunit, i, ialpha, N
         real(kind=8)                              :: my_dpdx, q0, shat, Delta_x, Delta_y
         real(kind=8), dimension(:), allocatable   :: y
         real(kind=8)                              :: z_min, z_max, dlnndx, dlnTdx
         real(kind=8)                              :: a, pi, B_ave, r1, r2, L_tot, omn, omt

         ! GIST/GENE namelist
         namelist /PARAMETERS/ my_dpdx, q0, shat, gridpoints, n_pol
      !-----------------------------------------------------------------
      !        Begin Subroutine
      !-----------------------------------------------------------------

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


         CALL compute_AE(N,g11,g12,g22,Bhat_array,abs_jac,L1,L2,dBdt, &
                      n_pol, my_dpdx, q0, shat, 3.0D0, 0.0D0, &
                      1D-4, 1.0D+1, 10000, 10000, 1D-10, AE_tot)

         RETURN
      !-----------------------------------------------------------------
      !        End Subroutine
      !-----------------------------------------------------------------
       end subroutine compute_AE_GIST


       SUBROUTINE compute_AE(N_arr,g11,g12,g22,Bhat,abs_jac,L1,L2,dBdt, &
                      n_pol, my_dpdx, q0, shat, omn, omt, &
                      z_min, z_max, z_res, lam_res, Delta_t, AE_tot)
          IMPLICIT NONE
      !-----------------------------------------------------------------
      !        Subroutine Input Variables
      !           g11      g_s-s GIST Array
      !           g12      g_s-theta GIST Array
      !           g22      g_theta-theta GIST Array
      !           Bhat     Normalized magnetic field Array
      !           abs_jac  Normalized Jacobian
      !           L1
      !           L2
      !           dBdt     dB/dtheta Array
      !           n_pol    Poloidal Turns (usually 1)
      !           my_dpdx  Normalized pressure gradient
      !           q0       1/iota
      !           shat     Normalized radial coordinate
      !           omn      Density Gradient (GENE Deff)
      !           omt      Temperature Gradient (GENE Deff)
      !           z_min    Minimum normalized Energy (1E-5)
      !           z_max    Maximum normalized Energy (40)
      !           z_res    Resolution of Integration Range (10000)
      !           lam_res  Resolution of Pitch Angle Integration (10000)
      !           Delta_t  Padding for periodicity boundary condition (1E-10)
      !           AE_tot :     Total Free Energy.
      !-----------------------------------------------------------------
          integer, intent(in)                       :: N_arr, n_pol, z_res, lam_res
          real(kind=8), intent(in)                  :: my_dpdx, q0, shat, omn, omt
          real(kind=8), intent(in)                  :: z_min, z_max, Delta_t
          real(kind=8), dimension(N_arr), intent(in):: g11, g12, g22, Bhat, abs_jac, L2, L1, dBdt
          real(kind=8), intent(out)                 :: AE_tot
      !-----------------------------------------------------------------
      !        Subroutine Variables
      !-----------------------------------------------------------------
          real(kind=8)                              :: Delta_x, Delta_y
          real(kind=8), dimension(N_arr)            :: y, sqrt_g
          real(kind=8)                              :: dlnndx, dlnTdx
          real(kind=8)                              :: a, B_ave, r1, r2, L_tot
          integer                                   :: N, i

      !-----------------------------------------------------------------
      !        Begin Subroutine
      !-----------------------------------------------------------------
          ! define scalars
          a = n_pol*2*4.D0*DATAN(1.D0) ! max of theta array
          N = N_arr


          y = (/((i*a)/(N-1),i=0,N-1)/) ! assign to y an array for theta [0,a]
          sqrt_g = 1.0/abs_jac


          ! Calculate average B (arclength-averaged)
          call trapz(y,q0*Bhat*Bhat*sqrt_g,N,r1)
          call trapz(y,q0*Bhat*sqrt_g,N,r2)
          B_ave = r1/r2
          Delta_x = q0
          Delta_y = q0
          dlnndx = -omn
          dlnTdx = -omt

          call trapz(y,q0*Bhat*sqrt_g/B_ave,N,L_tot)


          CALL AE_total(q0,dlnTdx,dlnndx,Delta_x,Delta_y,Bhat/B_ave, &
            L1/B_ave,L2/B_ave,sqrt_g,y,lam_res,z_res, &
            z_min,z_max,Delta_t,N,L_tot,AE_tot)

         RETURN
      !-----------------------------------------------------------------
      !        End Subroutine
      !-----------------------------------------------------------------
        END SUBROUTINE compute_AE

      END MODULE trapped_avail_energy_mod
