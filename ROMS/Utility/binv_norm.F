#include "cppdefs.h"
      MODULE binv_norm_mod

#if defined HESSIAN_SV && defined BINVNORM
!
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2013 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE

      PUBLIC :: binv_norm

      CONTAINS
!
!***********************************************************************
      SUBROUTINE binv_norm (ng, tile, Kinp, Kout, Ninp, Nout, StateNorm)
!***********************************************************************
!
      USE mod_param
# ifdef ADJUST_BOUNDARY
      USE mod_boundary
# endif
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      USE mod_forces
# endif
      USE mod_grid
      USE mod_ocean
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, Kinp, Kout, Ninp, Nout

      real(r8), intent(out) :: StateNorm
!
!  Local variable declarations.
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iTLM, 2)
# endif
      CALL binv_norm_tile (ng, tile,                                    &
     &                       LBi, UBi, LBj, UBj, LBij, UBij,            &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       Kinp, Kout, Ninp, Nout,                    &
     &                       StateNorm,                                 &
# ifdef MASKING
     &                       GRID(ng) % rmask,                          &
     &                       GRID(ng) % umask,                          &
     &                       GRID(ng) % vmask,                          &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                       BOUNDARY(ng) % t_obc,                      &
     &                       BOUNDARY(ng) % u_obc,                      &
     &                       BOUNDARY(ng) % v_obc,                      &
#  endif
     &                       BOUNDARY(ng) % ubar_obc,                   &
     &                       BOUNDARY(ng) % vbar_obc,                   &
     &                       BOUNDARY(ng) % zeta_obc,                   &
# endif
# ifdef ADJUST_WSTRESS
     &                       FORCES(ng) % ustr,                         &
     &                       FORCES(ng) % vstr,                         &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                       FORCES(ng) % tflux,                        &
# endif
# ifdef SOLVE3D
     &                       OCEAN(ng) % t,                             &
     &                       OCEAN(ng) % u,                             &
     &                       OCEAN(ng) % v,                             &
# else
     &                       OCEAN(ng) % ubar,                          &
     &                       OCEAN(ng) % vbar,                          &
# endif
     &                       OCEAN(ng) % zeta,                          &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                       BOUNDARY(ng) % ad_t_obc,                   &
     &                       BOUNDARY(ng) % ad_u_obc,                   &
     &                       BOUNDARY(ng) % ad_v_obc,                   &
#  endif
     &                       BOUNDARY(ng) % ad_ubar_obc,                &
     &                       BOUNDARY(ng) % ad_vbar_obc,                &
     &                       BOUNDARY(ng) % ad_zeta_obc,                &
# endif
# ifdef ADJUST_WSTRESS
     &                       FORCES(ng) % ad_ustr,                      &
     &                       FORCES(ng) % ad_vstr,                      &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                       FORCES(ng) % ad_tflux,                     &
# endif
# ifdef SOLVE3D
     &                       OCEAN(ng) % ad_t,                          &
     &                       OCEAN(ng) % ad_u,                          &
     &                       OCEAN(ng) % ad_v,                          &
# else
     &                       OCEAN(ng) % ad_ubar,                       &
     &                       OCEAN(ng) % ad_vbar,                       &
# endif
     &                       OCEAN(ng) % ad_zeta,                       &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                       BOUNDARY(ng) % tl_t_obc,                   &
     &                       BOUNDARY(ng) % tl_u_obc,                   &
     &                       BOUNDARY(ng) % tl_v_obc,                   &
#  endif
     &                       BOUNDARY(ng) % tl_ubar_obc,                &
     &                       BOUNDARY(ng) % tl_vbar_obc,                &
     &                       BOUNDARY(ng) % tl_zeta_obc,                &
# endif
# ifdef ADJUST_WSTRESS
     &                       FORCES(ng) % tl_ustr,                      &
     &                       FORCES(ng) % tl_vstr,                      &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                       FORCES(ng) % tl_tflux,                     &
# endif
# ifdef SOLVE3D
     &                       OCEAN(ng) % tl_t,                          &
     &                       OCEAN(ng) % tl_u,                          &
     &                       OCEAN(ng) % tl_v,                          &
# else
     &                       OCEAN(ng) % tl_ubar,                       &
     &                       OCEAN(ng) % tl_vbar,                       &
# endif
     &                       OCEAN(ng) % tl_zeta)
# ifdef PROFILE
      CALL wclock_off (ng, iTLM, 2)
# endif

      RETURN
      END SUBROUTINE binv_norm
!
!***********************************************************************
      SUBROUTINE binv_norm_tile (ng, tile,                              &
     &                             LBi, UBi, LBj, UBj, LBij, UBij,      &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             Kinp, Kout, Ninp, Nout,              &
     &                             StateNorm,                           &
# ifdef MASKING
     &                             rmask, umask, vmask,                 &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                             nl_t_obc, nl_u_obc, nl_v_obc,        &
#  endif
     &                             nl_ubar_obc, nl_vbar_obc,            &
     &                             nl_zeta_obc,                         &
# endif
# ifdef ADJUST_WSTRESS
     &                             nl_ustr, nl_vstr,                    &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                             nl_tflux,                            &
# endif
# ifdef SOLVE3D
     &                             nl_t, nl_u, nl_v,                    &
# else
     &                             nl_ubar, nl_vbar,                    &
# endif
     &                             nl_zeta,                             &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                             ad_t_obc, ad_u_obc, ad_v_obc,        &
#  endif
     &                             ad_ubar_obc, ad_vbar_obc,            &
     &                             ad_zeta_obc,                         &
# endif
# ifdef ADJUST_WSTRESS
     &                             ad_ustr, ad_vstr,                    &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                             ad_tflux,                            &
# endif
# ifdef SOLVE3D
     &                             ad_t, ad_u, ad_v,                    &
# else
     &                             ad_ubar, ad_vbar,                    &
# endif
     &                             ad_zeta,                             &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                             tl_t_obc, tl_u_obc, tl_v_obc,        &
#  endif
     &                             tl_ubar_obc, tl_vbar_obc,            &
     &                             tl_zeta_obc,                         &
# endif
# ifdef ADJUST_WSTRESS
     &                             tl_ustr, tl_vstr,                    &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                             tl_tflux,                            &
# endif
# ifdef SOLVE3D
     &                             tl_t, tl_u, tl_v,                    &
# else
     &                             tl_ubar, tl_vbar,                    &
# endif
     &                             tl_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
!
      USE state_addition_mod, ONLY : state_addition
      USE state_dotprod_mod, ONLY : state_dotprod
      USE state_initialize_mod, ONLY : state_initialize
      USE state_scale_mod, ONLY : state_scale
      USE inner2state_mod, ONLY : read_state
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBij, UBij
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Kinp, Kout, Ninp, Nout
      real(r8), intent(out) :: StateNorm
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: nl_t_obc(LBij:,:,:,:,:,:)
      real(r8), intent(inout) :: nl_u_obc(LBij:,:,:,:,:)
      real(r8), intent(inout) :: nl_v_obc(LBij:,:,:,:,:)
#   endif
      real(r8), intent(inout) :: nl_ubar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: nl_vbar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: nl_zeta_obc(LBij:,:,:,:)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: nl_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: nl_vstr(LBi:,LBj:,:,:)
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(inout) :: nl_tflux(LBi:,LBj:,:,:,:)
#  endif
#  ifdef SOLVE3D
      real(r8), intent(inout) :: nl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: nl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: nl_v(LBi:,LBj:,:,:)
#  else
      real(r8), intent(inout) :: nl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: nl_vbar(LBi:,LBj:,:)
#  endif
      real(r8), intent(inout) :: nl_zeta(LBi:,LBj:,:)
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t_obc(LBij:,:,:,:,:,:)
      real(r8), intent(inout) :: ad_u_obc(LBij:,:,:,:,:)
      real(r8), intent(inout) :: ad_v_obc(LBij:,:,:,:,:)
#   endif
      real(r8), intent(inout) :: ad_ubar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: ad_vbar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: ad_zeta_obc(LBij:,:,:,:)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ad_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_vstr(LBi:,LBj:,:,:)
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(inout) :: ad_tflux(LBi:,LBj:,:,:,:)
#  endif
#  ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_v(LBi:,LBj:,:,:)
#  else
      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
#  endif
      real(r8), intent(inout) :: ad_zeta(LBi:,LBj:,:)
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t_obc(LBij:,:,:,:,:,:)
      real(r8), intent(inout) :: tl_u_obc(LBij:,:,:,:,:)
      real(r8), intent(inout) :: tl_v_obc(LBij:,:,:,:,:)
#   endif
      real(r8), intent(inout) :: tl_ubar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: tl_vbar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: tl_zeta_obc(LBij:,:,:,:)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: tl_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_vstr(LBi:,LBj:,:,:)
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(inout) :: tl_tflux(LBi:,LBj:,:,:,:)
#  endif
#  ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_v(LBi:,LBj:,:,:)
#  else
      real(r8), intent(inout) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_vbar(LBi:,LBj:,:)
#  endif
      real(r8), intent(inout) :: tl_zeta(LBi:,LBj:,:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: nl_t_obc(LBij:UBij,N(ng),4,            &
     &                                    Nbrec(ng),2,NT(ng))
      real(r8), intent(inout) :: nl_u_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
      real(r8), intent(inout) :: nl_v_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
#   endif
      real(r8), intent(inout) :: nl_ubar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: nl_vbar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: nl_zeta_obc(LBij:UBij,4,Nbrec(ng),2)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: nl_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: nl_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(inout) :: nl_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
#  endif
#  ifdef SOLVE3D
      real(r8), intent(inout) :: nl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: nl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: nl_v(LBi:UBi,LBj:UBj,N(ng),2)
#  else
      real(r8), intent(inout) :: nl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: nl_vbar(LBi:UBi,LBj:UBj,3)
#  endif
      real(r8), intent(inout) :: nl_zeta(LBi:UBi,LBj:UBj,3)
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t_obc(LBij:UBij,N(ng),4,            &
     &                                    Nbrec(ng),2,NT(ng))
      real(r8), intent(inout) :: ad_u_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
      real(r8), intent(inout) :: ad_v_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
#   endif
      real(r8), intent(inout) :: ad_ubar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: ad_vbar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: ad_zeta_obc(LBij:UBij,4,Nbrec(ng),2)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ad_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: ad_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(inout) :: ad_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
#  endif
#  ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
#  else
      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
#  endif
      real(r8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,3)
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t_obc(LBij:UBij,N(ng),4,            &
     &                                    Nbrec(ng),2,NT(ng))
      real(r8), intent(inout) :: tl_u_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
      real(r8), intent(inout) :: tl_v_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
#   endif
      real(r8), intent(inout) :: tl_ubar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: tl_vbar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: tl_zeta_obc(LBij:UBij,4,Nbrec(ng),2)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: tl_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: tl_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(inout) :: tl_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
#  endif
#  ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: tl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: tl_v(LBi:UBi,LBj:UBj,N(ng),2)
#  else
      real(r8), intent(inout) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: tl_vbar(LBi:UBi,LBj:UBj,3)
#  endif
      real(r8), intent(inout) :: tl_zeta(LBi:UBi,LBj:UBj,3)
# endif
!
!  Local variable declarations.
!
      integer :: Lwrk1, Lwrk2, i, j, lstr, ndefERR, outLoop, rec
      integer :: neof
# ifdef SOLVE3D
      integer :: itrc, k
# endif
      real(r8) :: fac, fac1, fac2, Rerrtol

      real(r8), dimension(0:NstateVar(ng)) :: dot
      real(r8), dimension(NpostI) :: Reigval
      real(r8), dimension(NpostI) :: Rerrval

      character (len=256) :: ncname

# include "set_bounds.h"
!
      SourceFile='binv_norm.F, binv_norm_tile'
!
      Rerrtol=1.0e-5
!
!-----------------------------------------------------------------------
!  Compute adjont model initial conditions from the inv(B) using the 
!  EOFs. *** NOTE WE ARE USING THE ERR file name as storage for the
!  EOFs.
!-----------------------------------------------------------------------
!
!  Determine if single or multiple Lanczos vector NetCDF files.
!
      CALL netcdf_get_ivar (ng, iADM, TRIM(ERR(ng)%name), 'ndefADJ',    &
     &                      ndefERR)
      IF (exit_flag.ne.NoError) RETURN
!
!  Read the Ritz values.
!
      CALL netcdf_get_fvar (ng, iADM, ERR(ng)%name, 'Ritz',             &
     &                        Reigval)
      IF (exit_flag.ne. NoError) RETURN
      CALL netcdf_get_fvar (ng, iADM, ERR(ng)%name, 'Ritz_error',       &
     &                        Rerrval)
      IF (exit_flag.ne. NoError) RETURN

!
      Lwrk1=1
      Lwrk2=2
!
!  Initialize nl_var state arrays: nl_var(Lini) = fac
!
      fac=0.0_r8

      CALL state_initialize (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj, LBij, UBij,            &
     &                       Lwrk2, fac,                                &
# ifdef MASKING
     &                       rmask, umask, vmask,                       &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                       nl_t_obc, nl_u_obc, nl_v_obc,              &
#  endif
     &                       nl_ubar_obc, nl_vbar_obc,                  &
     &                       nl_zeta_obc,                               &
# endif
# ifdef ADJUST_WSTRESS
     &                       nl_ustr, nl_vstr,                          &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                       nl_tflux,                                  &
# endif
# ifdef SOLVE3D
     &                       nl_t, nl_u, nl_v,                          &
# else
     &                       nl_ubar, nl_vbar,                          &
# endif
     &                       nl_zeta)
!
      DO neof=1,NpostI/2
!
!  Determine EOF file to read.
!
        IF (ndefERR.gt.0) THEN
          lstr=LEN_TRIM(ERR(ng)%name)
          WRITE (ncname,10) ERR(ng)%name(1:lstr-8), neof
 10       FORMAT (a,'_',i4.4,'.nc')
        ELSE
          ncname=ERR(ng)%name
        END IF
!
!  Only use the EOFs that have converged sufficiently.
!
!AMM        IF (ABS(Rerrval(neof)).lt.Rerrtol) THEN
!
!  Read in the EOFs into the nl_var arrays
!
        CALL read_state (ng, tile, iTLM,                                &
     &                   LBi, UBi, LBj, UBj, LBij, UBij,                &
     &                   Lwrk1, neof,                                   &
     &                   ndefERR, ERR(ng)%ncid, TRIM(ncname),           &
# ifdef MASKING
     &                   rmask, umask, vmask,                           &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                   nl_t_obc, nl_u_obc, nl_v_obc,                  &
#  endif
     &                   nl_ubar_obc, nl_vbar_obc,                      &
     &                   nl_zeta_obc,                                   &
# endif
# ifdef ADJUST_WSTRESS
     &                   nl_ustr, nl_vstr,                              &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                   nl_tflux,                                      &
# endif
# ifdef SOLVE3D
     &                   nl_t, nl_u, nl_v,                              &
# else
     &                   nl_ubar, nl_vbar,                              &
# endif
     &                   nl_zeta)
        IF (exit_flag.ne.NoError) RETURN
!
!  Compute dot product between the EOFs and the TL model solution.
!
        CALL state_dotprod (ng, tile, iTLM,                             &
     &                      LBi, UBi, LBj, UBj, LBij, UBij,             &
     &                      NstateVar(ng), dot(0:),                     &
# ifdef MASKING
     &                      rmask, umask, vmask,                        &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                      nl_t_obc(:,:,:,:,Lwrk1,:),                  &
     &                      tl_t_obc(:,:,:,:,Lwrk,:),                   &
     &                      nl_u_obc(:,:,:,:,Lwrk1),                    &
     &                      tl_u_obc(:,:,:,:,Lwrk),                     &
     &                      nl_v_obc(:,:,:,:,Lwrk1),                    &
     &                      tl_v_obc(:,:,:,:,Lwrk),                     &
#  endif
     &                      nl_ubar_obc(:,:,:,Lwrk1),                   &
     &                      tl_ubar_obc(:,:,:,Lwrk),                    &
     &                      nl_vbar_obc(:,:,:,Lwrk1),                   &
     &                      tl_vbar_obc(:,:,:,Lwrk),                    &
     &                      nl_zeta_obc(:,:,:,Lwrk1),                   &
     &                      tl_zeta_obc(:,:,:,Lwrk),                    &
# endif
# ifdef ADJUST_WSTRESS
     &                      nl_ustr(:,:,:,Lwrk1), tl_ustr(:,:,:,Lwrk),  &
     &                      nl_vstr(:,:,:,Lwrk1), tl_vstr(:,:,:,Lwrk),  &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                      nl_tflux(:,:,:,Lwrk1,:),                    &
     &                      tl_tflux(:,:,:,Lwrk,:),                     &
# endif
# ifdef SOLVE3D
     &                      nl_t(:,:,:,Lwrk1,:), tl_t(:,:,:,Ninp,:),    &
     &                      nl_u(:,:,:,Lwrk1), tl_u(:,:,:,Ninp),        &
     &                      nl_v(:,:,:,Lwrk1), tl_v(:,:,:,Ninp),        &
# else
     &                      nl_ubar(:,:,Lwrk1), tl_ubar(:,:,Kinp),      &
     &                      nl_vbar(:,:,Lwrk1), tl_vbar(:,:,Kinp),      &
# endif
     &                      nl_zeta(:,:,Lwrk1), tl_zeta(:,:,Kinp))
!
!  Form the weighted sum of the EOFs.
!
!    nl_var(Lwrk2) = fac1 * nl_var(Lwrk2) + fac2 * nl_var(Lwrk1)
!
        fac1=1.0_r8
        fac2=dot(0)/Reigval(neof)

        CALL state_addition (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj, LBij, UBij,            &
     &                       Lwrk2, Lwrk1, Lwrk2, fac1, fac2,           &
# ifdef MASKING
     &                       rmask, umask, vmask,                       &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                       nl_t_obc, nl_t_obc,                        &
     &                       nl_u_obc, nl_u_obc,                        &
     &                       nl_v_obc, nl_v_obc,                        &
#  endif
     &                       nl_ubar_obc, nl_ubar_obc,                  &
     &                       nl_vbar_obc, nl_vbar_obc,                  &
     &                       nl_zeta_obc, nl_zeta_obc,                  &
# endif
# ifdef ADJUST_WSTRESS
     &                       nl_ustr, nl_ustr,                          &
     &                       nl_vstr, nl_vstr,                          &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                       nl_tflux, nl_tflux,                        &
# endif
# ifdef SOLVE3D
     &                       nl_t, nl_t,                                &
     &                       nl_u, nl_u,                                &
     &                       nl_v, nl_v,                                &
# else
     &                       nl_ubar, nl_ubar,                          &
     &                       nl_vbar, nl_vbar,                          &
# endif
     &                       nl_zeta, nl_zeta)
!AMM      END IF
      END DO
!
!  Compute the final time norm.
!
        CALL state_dotprod (ng, tile, iTLM,                             &
     &                      LBi, UBi, LBj, UBj, LBij, UBij,             &
     &                      NstateVar(ng), dot(0:),                     &
# ifdef MASKING
     &                      rmask, umask, vmask,                        &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                      nl_t_obc(:,:,:,:,Lwrk2,:),                  &
     &                      tl_t_obc(:,:,:,:,Lwrk,:),                   &
     &                      nl_u_obc(:,:,:,:,Lwrk2),                    &
     &                      tl_u_obc(:,:,:,:,Lwrk),                     &
     &                      nl_v_obc(:,:,:,:,Lwrk2),                    &
     &                      tl_v_obc(:,:,:,:,Lwrk),                     &
#  endif
     &                      nl_ubar_obc(:,:,:,Lwrk2),                   &
     &                      tl_ubar_obc(:,:,:,Lwrk),                    &
     &                      nl_vbar_obc(:,:,:,Lwrk2),                   &
     &                      tl_vbar_obc(:,:,:,Lwrk),                    &
     &                      nl_zeta_obc(:,:,:,Lwrk2),                   &
     &                      tl_zeta_obc(:,:,:,Lwrk),                    &
# endif
# ifdef ADJUST_WSTRESS
     &                      nl_ustr(:,:,:,Lwrk2), tl_ustr(:,:,:,Lwrk),  &
     &                      nl_vstr(:,:,:,Lwrk2), tl_vstr(:,:,:,Lwrk),  &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                      nl_tflux(:,:,:,Lwrk2,:),                    &
     &                      tl_tflux(:,:,:,Lwrk,:),                     &
# endif
# ifdef SOLVE3D
     &                      nl_t(:,:,:,Lwrk2,:), tl_t(:,:,:,Ninp,:),    &
     &                      nl_u(:,:,:,Lwrk2), tl_u(:,:,:,Ninp),        &
     &                      nl_v(:,:,:,Lwrk2), tl_v(:,:,:,Ninp),        &
# else
     &                      nl_ubar(:,:,Lwrk2), tl_ubar(:,:,Kinp),      &
     &                      nl_vbar(:,:,Lwrk2), tl_vbar(:,:,Kinp),      &
# endif
     &                      nl_zeta(:,:,Lwrk2), tl_zeta(:,:,Kinp))
!
      StateNorm=dot(0)
!
!
!-----------------------------------------------------------------------
!  Initialize adjoint state variables with nl_var arrays.
!-----------------------------------------------------------------------

#  ifdef FULL_GRID
#   define I_R_RANGE IstrT,IendT
#   define I_U_RANGE IstrP,IendT
#   define J_R_RANGE JstrT,JendT
#   define J_V_RANGE JstrP,JendT
#  else
#   define I_R_RANGE Istr,Iend
#   define I_U_RANGE IstrU,Iend
#   define J_R_RANGE Jstr,Jend
#   define J_V_RANGE JstrV,Jend
#  endif
!
!  Free-surface.
!
      DO j=J_R_RANGE
        DO i=I_R_RANGE
          ad_zeta(i,j,Kout)=nl_zeta(i,j,Lwrk2)
#  ifdef MASKING
          ad_zeta(i,j,Kout)=ad_zeta(i,j,Kout)*rmask(i,j)
#  endif
        END DO
      END DO

#  ifndef SOLVE3D
!
!  2D u-momentum component.
!
      DO j=J_R_RANGE
        DO i=I_U_RANGE
          ad_ubar(i,j,Kout)=nl_ubar(i,j,Lwrk2)
#   ifdef MASKING
          ad_ubar(i,j,Kout)=ad_ubar(i,j,Kout)*umask(i,j)
#   endif
        END DO
      END DO
!
!  2D v-momentum component.
!
      DO j=J_V_RANGE
        DO i=I_R_RANGE
          ad_vbar(i,j,Kout)=nl_vbar(i,j,Lwrk2)
#   ifdef MASKING
          ad_vbar(i,j,Kout)=ad_vbar(i,j,Kout)*vmask(i,j)
#   endif
        END DO
      END DO
#  else
!
!  3D u-momentum component.
!
      DO k=1,N(ng)
        DO j=J_R_RANGE
          DO i=I_U_RANGE
            ad_u(i,j,k,Nout)=nl_u(i,j,k,Lwrk2)
#   ifdef MASKING
            ad_u(i,j,k,Nout)=ad_u(i,j,k,Nout)*umask(i,j)
#   endif
          END DO
        END DO
      END DO
!
!  3D v-momentum component.
!
      DO k=1,N(ng)
        DO j=J_V_RANGE
          DO i=I_R_RANGE
            ad_v(i,j,k,Nout)=nl_v(i,j,k,Lwrk2)
#   ifdef MASKING
            ad_v(i,j,k,Nout)=ad_v(i,j,k,Nout)*vmask(i,j)
#   endif
          END DO
        END DO
      END DO
!
!  Tracers.
!
      DO itrc=1,NT(ng)
        DO k=1,N(ng)
          DO j=J_R_RANGE
            DO i=I_R_RANGE
              ad_t(i,j,k,Nout,itrc)=nl_t(i,j,k,Lwrk2,itrc)
#   ifdef MASKING
              ad_t(i,j,k,Nout,itrc)=ad_t(i,j,k,Nout,itrc)*rmask(i,j)
#   endif
            END DO
          END DO
        END DO
      END DO
#  endif

#  undef I_R_RANGE
#  undef I_U_RANGE
#  undef J_R_RANGE
#  undef J_V_RANGE

!
      RETURN
      END SUBROUTINE binv_norm_tile
#endif
      END MODULE binv_norm_mod

