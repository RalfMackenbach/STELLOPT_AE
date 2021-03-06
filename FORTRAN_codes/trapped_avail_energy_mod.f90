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
      INTEGER, PRIVATE      :: z_res, lam_res
      REAL(KIND=8), PRIVATE :: omn, omt, z_min, z_max, Delta_t


      NAMELIST /AVAIL_ENERGY_OPTIONS/ omn, omt, z_min, z_max, Delta_t, z_res, lam_res

      PRIVATE :: AVAIL_ENERGY_OPTIONS

!-----------------------------------------------------------------------
!     Private Subroutines
!-----------------------------------------------------------------------
      CONTAINS

!-----------------------------------------------------------------------
!     Public Subroutines
!-----------------------------------------------------------------------

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
         real(kind=8), dimension(:,:),allocatable  :: g11, g12, g22, Bhat, abs_jac, dBdx, dBdy, dBdt, dummy
         real(kind=8), dimension(:),  allocatable  :: Bhat_array, sqrtg_arr, dBdx_array, dBdy_array

         integer                                   :: gridpoints, n_pol, iunit, i, ialpha, N
         real(kind=8)                              :: my_dpdx, q0, shat, Delta_x, Delta_y
         real(kind=8), dimension(:), allocatable   :: y
         real(kind=8)                              :: dlnndx, dlnTdx
         real(kind=8)                              :: a, pi, B_ave, r1, r2, L_tot

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
         CALL read_avail_energy_nml(iunit,i)
         READ(iunit,NML=PARAMETERS)
         allocate(g11(1,gridpoints), g12(1,gridpoints), g22(1,gridpoints),  &
                  Bhat(1,gridpoints), abs_jac(1,gridpoints), dBdx(1,gridpoints), &
                  dBdy(1,gridpoints), dBdt(1,gridpoints), dummy(1,gridpoints))
         DO i = 1, gridpoints
            READ(iunit,"(9ES20.10)") g11(ialpha,i),g12(ialpha,i),g22(ialpha,i),&
                                     Bhat(ialpha,i),abs_jac(ialpha,i),dBdx(ialpha,i),&
                                     dBdy(ialpha,i),dBdt(ialpha,i), dummy(ialpha,i)
         END DO
         CLOSE(iunit)

         ! define scalars
         pi= 4.D0*DATAN(1.D0)
         a = n_pol*2*pi ! max of theta array
         N = size(Bhat)

         ! allocate data
         allocate(y(N),Bhat_array(N),dBdx_array(N),dBdy_array(N),sqrtg_arr(N))
         ! Reshape arrays
         Bhat_array    = reshape(Bhat(1,:),    shape(Bhat_array))
         dBdx_array    = reshape(dBdx(1,:),    shape(Bhat_array))
         dBdy_array    = reshape(dBdy(1,:),    shape(Bhat_array))
         sqrtg_arr = 1.0/reshape(abs_jac(1,:), shape(Bhat_array))
         y = (/((i*a)/(N-1),i=0,N-1)/) ! assign to y an array for theta [0,a]
         ! print *,y


         CALL compute_AE(N,g11,g12,g22,Bhat_array,abs_jac,dBdx_array,dBdy_array,dBdt, &
                      n_pol, my_dpdx, q0, shat, AE_tot)

         RETURN
      !-----------------------------------------------------------------
      !        End Subroutine
      !-----------------------------------------------------------------
       end subroutine compute_AE_GIST


       SUBROUTINE compute_AE(N_arr,g11,g12,g22,Bhat,abs_jac,dBdx,dBdy,dBdt, &
                      n_pol, my_dpdx, q0, shat, AE_tot)
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
          integer, intent(in)                       :: N_arr, n_pol
          real(kind=8), intent(in)                  :: my_dpdx, q0, shat
          real(kind=8), dimension(N_arr), intent(in):: g11, g12, g22, Bhat, abs_jac, dBdx, dBdy, dBdt
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
            dBdx/B_ave,dBdy/B_ave,sqrt_g,y,lam_res,z_res, &
            z_min,z_max,Delta_t,N,L_tot,AE_tot)

         RETURN
      !-----------------------------------------------------------------
      !        End Subroutine
      !-----------------------------------------------------------------
        END SUBROUTINE compute_AE

      SUBROUTINE read_avail_energy_nml(iunit,ier)
        IMPLICIT NONE
      !-----------------------------------------------------------------
      !        Subroutine Input Variables
      !           iunit :  File unit number.
      !           ier :    Error flag
      !-----------------------------------------------------------------
         INTEGER, intent(INOUT)  :: iunit
         INTEGER, intent(OUT)    :: ier
      !-----------------------------------------------------------------
      !        Begin Subroutine
      !-----------------------------------------------------------------
         ! Defaults
         omn = 4.0D0
         omt = 0.0D0
         z_min = 1.0D-4
         z_max = 4.0D+1
         Delta_t = 1D-10
         z_res   = 1000
         lam_res = 1000
         ier     = 0

         ! Read the Namelist
         READ(iunit, NML=AVAIL_ENERGY_OPTIONS, IOSTAT=ier)
         REWIND(iunit)

         RETURN
      END SUBROUTINE read_avail_energy_nml

      SUBROUTINE write_avail_energy_nml(iunit,ier)
        IMPLICIT NONE
      !-----------------------------------------------------------------
      !        Subroutine Input Variables
      !           iunit :  File unit number.
      !           ier :    Error flag
      !-----------------------------------------------------------------
         INTEGER, intent(INOUT)  :: iunit
         INTEGER, intent(OUT)    :: ier
      !-----------------------------------------------------------------
      !        Subroutine Variables
      !-----------------------------------------------------------------
      !-----------------------------------------------------------------
      !        Begin Subroutine
      !-----------------------------------------------------------------
         WRITE(iunit,'(A)') '&AVAIL_ENERGY_OPTIONS'
         WRITE(iunit,'(2X,A,D20.10)') 'OMN = ',omn
         WRITE(iunit,'(2X,A,D20.10)') 'OMT = ',omt
         WRITE(iunit,'(2X,A,D20.10)') 'Z_MIN = ',z_min
         WRITE(iunit,'(2X,A,D20.10)') 'Z_MAX = ',z_max
         WRITE(iunit,'(2X,A,I7)')     'Z_RES = ',z_res
         WRITE(iunit,'(2X,A,I7)')     'LAM_RES = ',lam_res
         WRITE(iunit,'(2X,A,D20.10)') 'DELTA_T = ',Delta_t
         WRITE(iunit,'(A)') '/'

         RETURN
      END SUBROUTINE write_avail_energy_nml

      END MODULE trapped_avail_energy_mod
