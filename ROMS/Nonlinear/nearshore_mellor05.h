#undef GEO_ROTATION
      SUBROUTINE radiation_stress (ng, tile)
!
!svn $Id: nearshore_mellor05.h 1748 2018-02-10 03:25:17Z arango $
!***********************************************************************
!  Copyright (c) 2002-2018 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!**************************************************   John C. Warner ***
!                                                                      !
!  This routine computes the radiation stress contributions to the     !
!  momentum equations using Mellor (2005) formulation.                 !
!                                                                      !
!  References:                                                         !
!                                                                      !
!  Mellor, G. L., 2003: The three-dimensional current and surface      !
!    wave equations, Journal of Physical Oceanography 33, 1978-1989.   !
!                                                                      !
!  Mellor, G. L., 2005: Some consequences of the three-dimensional     !
!    currents and surface wave equations, Journal of Physical          !
!    Oceanography 35, 2291-2298.                                       !
!                                                                      !
!***********************************************************************
!
      USE mod_forces
      USE mod_grid
      USE mod_mixing
      USE mod_ocean
      USE mod_stepping
      USE mod_coupling
#if defined DIAGNOSTICS_UV
      USE mod_diags
#endif
!
      integer, intent(in) :: ng, tile

#include "tile.h"

#ifdef PROFILE
      CALL wclock_on (ng, iNLM, 21, __LINE__, __FILE__)
#endif
      CALL radiation_stress_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, N(ng),            &
     &                            IminS, ImaxS, JminS, JmaxS,           &
#ifdef SOLVE3D
     &                            nrhs(ng),                             &
#endif
#ifdef MASKING
     &                            GRID(ng) % pmask,                     &
     &                            GRID(ng) % rmask,                     &
     &                            GRID(ng) % umask,                     &
     &                            GRID(ng) % vmask,                     &
#endif
#ifdef WET_DRY
     &                            GRID(ng) % umask_wet,                 &
     &                            GRID(ng) % vmask_wet,                 &
#endif
     &                            GRID(ng) % om_u,                      &
     &                            GRID(ng) % om_v,                      &
     &                            GRID(ng) % on_u,                      &
     &                            GRID(ng) % on_v,                      &
     &                            GRID(ng) % pm,                        &
     &                            GRID(ng) % pn,                        &
     &                            GRID(ng) % angler,                    &
#if defined CURVGRID
     &                            GRID(ng) % dndx,                      &
     &                            GRID(ng) % dmde,                      &
#endif
     &                            GRID(ng) % h,                         &
     &                            OCEAN(ng) % zeta,                     &
#ifdef SOLVE3D
     &                            GRID(ng) % Hz,                        &
     &                            GRID(ng) % z_r,                       &
     &                            GRID(ng) % z_w,                       &
#endif
     &                            FORCES(ng) % Hwave,                   &
     &                            FORCES(ng) % Dwave,                   &
     &                            FORCES(ng) % Lwave,                   &
#ifdef SVENDSEN_ROLLER
     &                            FORCES(ng) % Wave_break,              &
#endif
#ifdef SOLVE3D
# ifdef DIAGNOSTICS_UV
     &                            DIAGS(ng) % DiaRU,                    &
     &                            DIAGS(ng) % DiaRV,                    &
# endif
     &                            MIXING(ng) % Sxx,                     &
     &                            MIXING(ng) % Sxy,                     &
     &                            MIXING(ng) % Syy,                     &
     &                            MIXING(ng) % Szx,                     &
     &                            MIXING(ng) % Szy,                     &
     &                            MIXING(ng) % rustr3d,                 &
     &                            MIXING(ng) % rvstr3d,                 &
     &                            OCEAN(ng) % rulag3d,                  &
     &                            OCEAN(ng) % rvlag3d,                  &
     &                            OCEAN(ng) % u_stokes,                 &
     &                            OCEAN(ng) % v_stokes,                 &
#endif
     &                            MIXING(ng) % Sxx_bar,                 &
     &                            MIXING(ng) % Sxy_bar,                 &
     &                            MIXING(ng) % Syy_bar,                 &
     &                            MIXING(ng) % rustr2d,                 &
     &                            MIXING(ng) % rvstr2d,                 &
     &                            OCEAN(ng) % rulag2d,                  &
     &                            OCEAN(ng) % rvlag2d,                  &
     &                            OCEAN(ng) % ubar_stokes,              &
     &                            OCEAN(ng) % vbar_stokes)
#ifdef PROFILE
      CALL wclock_off (ng, iNLM, 21, __LINE__, __FILE__)
#endif
      RETURN
      END SUBROUTINE radiation_stress
!
!***********************************************************************
      SUBROUTINE radiation_stress_tile (ng, tile,                       &
     &                                  LBi, UBi, LBj, UBj, UBk,        &
     &                                  IminS, ImaxS, JminS, JmaxS,     &
#ifdef SOLVE3D
     &                                  nrhs,                           &
#endif
#ifdef MASKING
     &                                  pmask, rmask, umask, vmask,     &
#endif
#ifdef WET_DRY
     &                                  umask_wet, vmask_wet,           &
#endif
     &                                  om_u, om_v, on_u, on_v,         &
     &                                  pm, pn,                         &
     &                                  angler,                         &
#if defined CURVGRID
     &                                  dndx, dmde,                     &
#endif
     &                                  h, zeta,                        &
#ifdef SOLVE3D
     &                                  Hz, z_r, z_w,                   &
#endif
     &                                  Hwave, Dwave, Lwave,            &
#ifdef SVENDSEN_ROLLER
     &                                  Wave_break,                     &
#endif
#ifdef SOLVE3D
# ifdef DIAGNOSTICS_UV
     &                                  DiaRU, DiaRV,                   &
# endif
     &                                  Sxx, Sxy, Syy, Szx, Szy,        &
     &                                  rustr3d,  rvstr3d,              &
     &                                  rulag3d, rvlag3d,               &
     &                                  u_stokes, v_stokes,             &
#endif
     &                                  Sxx_bar, Sxy_bar, Syy_bar,      &
     &                                  rustr2d, rvstr2d,               &
     &                                  rulag2d, rvlag2d,               &
     &                                  ubar_stokes, vbar_stokes)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE exchange_2d_mod
      USE exchange_3d_mod
#ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d, mp_exchange3d
#endif
      USE bc_2d_mod
      USE bc_3d_mod
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, UBk
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
#ifdef SOLVE3D
      integer, intent(in) :: nrhs
#endif

#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: pmask(LBi:,LBj:)
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
# endif
# ifdef WET_DRY
      real(r8), intent(in) :: umask_wet(LBi:,LBj:)
      real(r8), intent(in) :: vmask_wet(LBi:,LBj:)
# endif
      real(r8), intent(in) :: om_u(LBi:,LBj:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
      real(r8), intent(in) :: on_v(LBi:,LBj:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: angler(LBi:,LBj:)
# if defined CURVGRID
      real(r8), intent(in) :: dndx(LBi:,LBj:)
      real(r8), intent(in) :: dmde(LBi:,LBj:)
# endif
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)
# ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
# endif
      real(r8), intent(in) :: Hwave(LBi:,LBj:)
      real(r8), intent(in) :: Dwave(LBi:,LBj:)
      real(r8), intent(in) :: Lwave(LBi:,LBj:)
# ifdef SVENDSEN_ROLLER
      real(r8), intent(in) :: Wave_break(LBi:,LBj:)
# endif
# ifdef SOLVE3D
#  ifdef DIAGNOSTICS_UV
      real(r8), intent(inout) :: DiaRU(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: DiaRV(LBi:,LBj:,:,:,:)
#  endif
      real(r8), intent(inout) :: Sxx(LBi:,LBj:,:)
      real(r8), intent(inout) :: Sxy(LBi:,LBj:,:)
      real(r8), intent(inout) :: Syy(LBi:,LBj:,:)
      real(r8), intent(inout) :: Szx(LBi:,LBj:,:)
      real(r8), intent(inout) :: Szy(LBi:,LBj:,:)
      real(r8), intent(inout) :: rustr3d(LBi:,LBj:,:)
      real(r8), intent(inout) :: rvstr3d(LBi:,LBj:,:)
      real(r8), intent(inout) :: rulag3d(LBi:,LBj:,:)
      real(r8), intent(inout) :: rvlag3d(LBi:,LBj:,:)
      real(r8), intent(inout) :: u_stokes(LBi:,LBj:,:)
      real(r8), intent(inout) :: v_stokes(LBi:,LBj:,:)
# endif
      real(r8), intent(inout) :: Sxx_bar(LBi:,LBj:)
      real(r8), intent(inout) :: Sxy_bar(LBi:,LBj:)
      real(r8), intent(inout) :: Syy_bar(LBi:,LBj:)
      real(r8), intent(inout) :: rustr2d(LBi:,LBj:)
      real(r8), intent(inout) :: rvstr2d(LBi:,LBj:)
      real(r8), intent(inout) :: rulag2d(LBi:,LBj:)
      real(r8), intent(inout) :: rvlag2d(LBi:,LBj:)
      real(r8), intent(inout) :: ubar_stokes(LBi:,LBj:)
      real(r8), intent(inout) :: vbar_stokes(LBi:,LBj:)

#else

# ifdef MASKING
      real(r8), intent(in) :: pmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
# endif
# ifdef WET_DRY
      real(r8), intent(in) :: umask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask_wet(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: om_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: angler(LBi:UBi,LBj:UBj)
# if defined CURVGRID
      real(r8), intent(in) :: dndx(LBi:UBI,LBj:UBj)
      real(r8), intent(in) :: dmde(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: zeta(LBi:UBi,LBj:UBj,3)
# ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:UBk)
# endif
      real(r8), intent(in) :: Hwave(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Dwave(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Lwave(LBi:UBi,LBj:UBj)
# ifdef SVENDSEN_ROLLER
      real(r8), intent(in) :: Wave_break(LBi:UBi,LBj:UBj)
# endif
# ifdef SOLVE3D
#  ifdef DIAGNOSTICS_UV
      real(r8), intent(inout) :: DiaRU(LBi:UBi,LBj:UBj,N(ng),2,NDrhs)
      real(r8), intent(inout) :: DiaRV(LBi:UBi,LBj:UBj,N(ng),2,NDrhs)
#  endif
      real(r8), intent(inout) :: Sxx(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: Sxy(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: Syy(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: Szx(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: Szy(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: rustr3d(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: rvstr3d(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: rulag3d(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: rvlag3d(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: u_stokes(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: v_stokes(LBi:UBi,LBj:UBj,N(ng))
# endif
      real(r8), intent(inout) :: Sxx_bar(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: Sxy_bar(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: Syy_bar(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: rustr2d(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: rvstr2d(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: rulag2d(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: rvlag2d(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ubar_stokes(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: vbar_stokes(LBi:UBi,LBj:UBj)
#endif
!
!  Local variable declarations.
!
      integer :: i, j, k

      real(r8) :: cff, cff1, cff2, cff3, cff4, cff5, cff6
      real(r8) :: fac1, fac2, FCCr, FCSr, FSSr

      real(r8), parameter :: eps = 1.0E-14_r8
      real(r8), parameter :: kDmax = 5.0_r8
      real(r8), parameter :: Lwave_min = 1.0_r8

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Dstp
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: kD
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: wavec
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: wavecgoc
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: waven
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: owaven
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: wavenx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: waveny
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: waveE
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: waveEr
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: gamw
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: UFe
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Sxxl
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Sxyl
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Syyl
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Sxyl_psi

#ifdef SOLVE3D
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: CF
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FC
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FCC
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FCS
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FSS

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: rollA
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: oroller
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ocosh
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: osinh
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: o2sinh
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: z_psi
# ifdef CURVGRID
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: UFx
# endif

#endif

#include "set_bounds.h"

      fac1=1.0_r8/dt(ng)

#ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Compute radiation stress components Sxx, Sxy (= Syx), and Syy.
!  Compute stokes velocities u_stokes and v_stokes.
!  Preliminary step.
!-----------------------------------------------------------------------
!
      DO j=Jstr-1,Jend+1
        DO i=Istr-1,Iend+1
!
!  Compute total depth.
!
          Dstp(i,j)=zeta(i,j,1)+h(i,j)
!
!  Compute wave numbers and wave energy.
!
          waven (i,j)=2.0_r8*pi/MAX(Lwave(i,j),Lwave_min)
          owaven(i,j)=1.0_r8/waven(i,j)
          cff=1.5_r8*pi-Dwave(i,j)-angler(i,j)
          wavenx(i,j)=waven(i,j)*COS(cff)
          waveny(i,j)=waven(i,j)*SIN(cff)
          waveE(i,j)=0.0625_r8*g*Hwave(i,j)*Hwave(i,j)
!
!  Compute wave celerity and group velocity.
!
          kD(i,j)=MIN(waven(i,j)*Dstp(i,j)+eps,kDmax)
          wavec(i,j)=SQRT(g*owaven(i,j)*TANH(kD(i,j)))
!
!  Compute metrics for vertical structure functions.
!
          ocosh(i,j)=1.0_r8/COSH(kD(i,j))
          osinh(i,j)=1.0_r8/SINH(kD(i,j))
          o2sinh(i,j)=1.0_r8/SINH(2.0_r8*kD(i,j))

# ifdef SVENDSEN_ROLLER
!
!  Compute metrics for roller contribution.
!
          oroller(i,j)=0.0_r8
          gamw(i,j)=MIN(Dstp(i,j)/(Hwave(i,j)+eps),5.0_r8)
          DO k=1,N(ng)
            cff2=SCALARS(ng)%Cs_r(k)
!! tanh s^4 functional
            oroller(i,j)=oroller(i,j)+                                  &
     &                   Hz(i,j,k)*                                     &
     &                   (1.0_r8-TANH((2.0_r8*cff2*gamw(i,j))**4))
!! s^4 functional
!!             oroller(i,j)=oroller(i,j)+                               &
!!   &                      Hz(i,j,k)*(1.0_r8-cff2**4)
          END DO
          oroller(i,j)=1.0_r8/(oroller(i,j)+eps)

#  ifdef MONO_ROLLER
!
!  Here Wave_break is really wave_area.
!
          cff1=Wave_break(i,j)/MAX(Lwave(i,j),Lwave_min)
#  else
          cff1=0.0424_r8*Hwave(i,j)*Wave_break(i,j)
#  endif
          rollA(i,j)=cff1*wavec(i,j)*wavec(i,j)
# endif
!
!  Initialize depth independent arrays for summation.
!
        END DO
        IF ((j.ge.Jstr).and.(j.le.Jend)) THEN
          DO i=Istr,Iend
            Sxx_bar(i,j)=0.0_r8
            Sxy_bar(i,j)=0.0_r8
            Syy_bar(i,j)=0.0_r8
          END DO
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  Compute diagonal [Sxx,Syy] and off-diagonal [Sxy,Syx] components
!  of radiation stresses.
!-----------------------------------------------------------------------
!
      K_LOOP : DO k=1,N(ng)
        fac2=1.0_r8+SCALARS(ng)%Cs_r(k)
        DO j=Jstr-1,Jend+1
          DO i=Istr-1,Iend+1
            cff1=waven(i,j)*waveE(i,j)
            FCCr=COSH(kD(i,j)*fac2)*ocosh(i,j)
            FCSr=COSH(kD(i,j)*fac2)*osinh(i,j)
            FSSr=SINH(kD(i,j)*fac2)*osinh(i,j)
            waveEr(i,j)=cff1*FCSr*FCCr
# ifdef SVENDSEN_ROLLER
!! tanh s^4 functional
            cff3=(1.0_r8-TANH((2.0_r8*(fac2-1.0_r8)*gamw(i,j))**4))
!! s^4 functional
!!          cff3=(1.0_r8-(fac2-1.0_r8)**4)
            waveEr(i,j)=waveEr(i,j)+                                    &
     &                 cff3*rollA(i,j)*oroller(i,j)
# endif
!
!  Compute local copy of radiation stresses at RHO-points.
!
            cff4=waveE(i,j)*waven(i,j)*FCSr*(FCCr-FSSr)
            cff5=owaven(i,j)*owaven(i,j)
            Sxxl(i,j)=cff4+                                             &
     &                waveEr(i,j)*wavenx(i,j)*wavenx(i,j)*cff5
            Syyl(i,j)=cff4+                                             &
     &                waveEr(i,j)*waveny(i,j)*waveny(i,j)*cff5
            Sxyl(i,j)=waveEr(i,j)*wavenx(i,j)*waveny(i,j)*cff5
          END DO
!
!  Compute depth integrated terms and copy stress terms
!  from local arrays to state variables.
!
          IF ((j.ge.Jstr).and.(j.le.Jend)) THEN
            DO i=Istr,Iend
              Sxx(i,j,k)=Sxxl(i,j)
              Sxy(i,j,k)=Sxyl(i,j)
              Syy(i,j,k)=Syyl(i,j)
              Sxx_bar(i,j)=Sxx_bar(i,j)+Sxx(i,j,k)*Hz(i,j,k)
              Sxy_bar(i,j)=Sxy_bar(i,j)+Sxy(i,j,k)*Hz(i,j,k)
              Syy_bar(i,j)=Syy_bar(i,j)+Syy(i,j,k)*Hz(i,j,k)
            END DO
          END IF
        END DO
!
!  Average Sxy to PSI-points.
!
        DO j=Jstr,Jend+1
          DO i=Istr,Iend+1
            Sxyl_psi(i,j)=0.25_r8*(Sxyl(i-1,j-1)+Sxyl(i,j-1)+           &
     &                             Sxyl(i-1,j  )+Sxyl(i,j  ))
          END DO
        END DO
!
!  Compute Radiation Stress component for U-momentum.
!
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            cff=(Sxxl(i  ,j)*Hz(i  ,j,k)-                               &
     &           Sxxl(i-1,j)*Hz(i-1,j,k))*on_u(i,j)
            rustr3d(i,j,k)=cff
# ifdef DIAGNOSTICS_UV
            DiaRU(i,j,k,nrhs,M3hrad)=-cff
# endif
          END DO
        END DO
!
!  Add in cross stress term Sxy averaged at PSI-points.
!
        DO j=Jstr,Jend+1
          DO i=Istr,Iend+1
            UFe(i,j)=0.25_r8*(Hz(i,j  ,k)+Hz(i-1,j  ,k)+                &
     &                        Hz(i,j-1,k)+Hz(i-1,j-1,k))*               &
     &               Sxyl_psi(i,j)
          END DO
        END DO
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            cff=(UFe(i,j+1)-UFe(i,j))*om_u(i,j)
            rustr3d(i,j,k)=rustr3d(i,j,k)+cff
# ifdef DIAGNOSTICS_UV
            DiaRU(i,j,k,nrhs,M3hrad)=DiaRU(i,j,k,nrhs,M3hrad)-cff
# endif
          END DO
        END DO
!
!  Compute U-stokes velocity.
!
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            cff=fac1*om_u(i,j)*on_u(i,j)
            cff2=(waveE(i-1,j)+waveE(i,j))
            cff3=(kD(i-1,j)+kD(i,j))

# if defined SVENDSEN_ROLLER
#  ifdef MONO_ROLLER
!
!  Here Wave_break is really wave_area.
!
            cff4=1.0_r8/MAX(Lwave(i-1,j)+Lwave(i,j),Lwave_min)
            cff2=cff2+                                                  &
     &           g*cff4*(Dstp(i-1,j)+Dstp(i,j))*                        &
     &           (Wave_break(i-1,j)+Wave_break(i,j))
#  else
            cff2=cff2+                                                  &
     &           0.25_r8*0.0424_r8*g*(Hwave(i-1,j)+Hwave(i,j))*         &
     &           (Wave_break(i-1,j)+Wave_break(i,j))*                   &
     &           (Dstp(i-1,j)+Dstp(i,j))
#  endif
# endif
!
!  Store old value to compute tendency term.
!
            rulag3d(i,j,k)=u_stokes(i,j,k)
            u_stokes(i,j,k)=cff2*                                       &
     &                      (wavenx(i-1,j)+wavenx(i,j))/                &
     &                      (wavec (i-1,j)+wavec (i,j))*                &
     &                      COSH(cff3*fac2)*0.5_r8*                     &
     &                      (o2sinh(i-1,j)+o2sinh(i,j))
# ifdef MASKING
            u_stokes(i,j,k)=u_stokes(i,j,k)*umask(i,j)
# endif
# ifdef WET_DRY
            u_stokes(i,j,k)=u_stokes(i,j,k)*umask_wet(i,j)
# endif
!
!  Finalize computation of stokes tendency term.
!
            rulag3d(i,j,k)=0.5_r8*cff*                                  &
     &                     (Hz(i,j,k)+Hz(i-1,j,k))*                     &
     &                     (u_stokes(i,j,k)-rulag3d(i,j,k))
          END DO
        END DO
!
!  Compute Radiation Stress component to V-momentum.
!
        DO j=JstrV,Jend
          DO i=Istr,Iend
            cff=(Syyl(i,j  )*Hz(i,j  ,k)-                               &
     &           Syyl(i,j-1)*Hz(i,j-1,k))*om_v(i,j)
            rvstr3d(i,j,k)=cff
# ifdef DIAGNOSTICS_UV
            DiaRV(i,j,k,nrhs,M3hrad)=-cff
# endif
          END DO
        END DO
!
!  Add in cross stress term Sxy averaged at PSI-points.
!
        DO j=JstrV,Jend
          DO i=Istr,Iend
            cff=(UFe(i+1,j)-UFe(i,j))*on_v(i,j)
            rvstr3d(i,j,k)=rvstr3d(i,j,k)+cff
# ifdef DIAGNOSTICS_UV
            DiaRV(i,j,k,nrhs,M3hrad)=DiaRV(i,j,k,nrhs,M3hrad)-cff
# endif
          END DO
        END DO
!
!  Compute V-stokes velocity.
!
        DO j=JstrV,Jend
          DO i=Istr,Iend
            cff=fac1*om_v(i,j)*on_v(i,j)
            cff2=(waveE(i,j-1)+waveE(i,j))
            cff3=(kD(i,j-1)+kD(i,j))
# if defined SVENDSEN_ROLLER
#  ifdef MONO_ROLLER
!
!  Here Wave_break is really wave_area.
!
            cff4=1.0_r8/MAX(Lwave(i,j-1)+Lwave(i,j),Lwave_min)
            cff2=cff2+                                                  &
     &           g*cff4*(Dstp(i,j-1)+Dstp(i,j))*                        &
     &           (Wave_break(i,j-1)+Wave_break(i,j))
#  else
            cff2=cff2+                                                  &
     &           0.25_r8*0.0424_r8*g*(Hwave(i,j-1)+Hwave(i,j))*         &
     &           (Wave_break(i,j-1)+Wave_break(i,j))*                   &
     &           (Dstp(i,j-1)+Dstp(i,j))
#  endif
# endif
!
!  Store old value to compute tendency term.
!
             rvlag3d(i,j,k)=v_stokes(i,j,k)
             v_stokes(i,j,k)=cff2*                                      &
     &                       (waveny(i,j-1)+waveny(i,j))/               &
     &                       (wavec (i,j-1)+wavec (i,j))*               &
     &                       COSH(cff3*fac2)*                           &
     &                       0.5_r8*(o2sinh(i,j-1)+o2sinh(i,j))
# ifdef MASKING
              v_stokes(i,j,k)=v_stokes(i,j,k)*vmask(i,j)
# endif
# ifdef WET_DRY
              v_stokes(i,j,k)=v_stokes(i,j,k)*vmask_wet(i,j)
# endif
!
!  Finalize computation of stokes tendency term.
!
              rvlag3d(i,j,k)=0.5_r8*cff*                                &
     &                       (Hz(i,j,k)+Hz(i,j-1,k))*                   &
     &                       (v_stokes(i,j,k)-rvlag3d(i,j,k))
          END DO
        END DO

# if defined CURVGRID
!
!-----------------------------------------------------------------------
!  Add in curvilinear transformation terms.
!-----------------------------------------------------------------------
!
        DO j=Jstr,Jend
          DO i=IstrU-1,Iend
            UFx(i,j)=Hz(i,j,k)*                                         &
     &               (Sxyl(i,j)*dmde(i,j)-                              &
     &                Syyl(i,j)*dndx(i,j))
          END DO
          DO i=IstrU,Iend
            cff=0.5_r8*(UFx(i,j)+UFx(i-1,j))
            rustr3d(i,j,k)=rustr3d(i,j,k)+cff
#  ifdef DIAGNOSTICS_UV
            DiaRU(i,j,k,nrhs,M3hrad)=DiaRU(i,j,k,nrhs,M3hrad)-cff
#  endif
          END DO
        END DO
        DO i=Istr,Iend
          DO j=JstrV-1,Jend
            UFe(i,j)=Hz(i,j,k)*                                         &
     &               (Sxyl(i,j)*dndx(i,j)-                              &
     &                Sxxl(i,j)*dmde(i,j))
          END DO
          DO j=JstrV,Jend
            cff=0.5_r8*(UFe(i,j)+UFe(i,j-1))
            rvstr3d(i,j,k)=rvstr3d(i,j,k)+cff
#  ifdef DIAGNOSTICS_UV
            DiaRV(i,j,k,nrhs,M3hrad)=DiaRV(i,j,k,nrhs,M3hrad)-cff
#  endif
          END DO
        END DO
# endif
      END DO K_LOOP
# if defined GEO_ROTATION
!
!-----------------------------------------------------------------------
!  Add in vertical rotation along geopotential terms.
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend+1
        DO i=Istr,Iend+1
          z_psi(i,j)=0.25_r8*(z_r(i  ,j-1,k)+z_r(i  ,j,k)+              &
     &                        z_r(i-1,j-1,k)+z_r(i-1,j,k))
        END DO
      END DO
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          UFe(i,j)=0.5_r8*(kD(i,j)+kD(i-1,j))
          ocosh(i,j)=1.0_r8/COSH(UFe(i,j))
          osinh(i,j)=1.0_r8/SINH(UFe(i,j))
        END DO
        DO i=IstrU,Iend
          DO k=0,N(ng)
            cff1=0.5_r8*(waven(i  ,j)*waveE(i  ,j)+                     &
     &                   waven(i-1,j)*waveE(i-1,j))
            cff2=1.0_r8+SCALARS(ng)%Cs_w(k)
            FCCr=COSH(UFe(i,j)*cff2)*ocosh(i,j)
            FCSr=COSH(UFe(i,j)*cff2)*osinh(i,j)
            FSSr=SINH(UFe(i,j)*cff2)*osinh(i,j)
            waveEr(i,j)=cff1*FCSr*FCCr
#  ifdef SVENDSEN_ROLLER
            cff3=SCALARS(ng)%Cs_w(k)
            waveEr(i,j)=waveEr(i,j)+1.25_r8*(rollA(i,j)+rollA(i-1,j))*  &
     &                  (1.0_r8-cff3**4)/(Dstp(i,j)+Dstp(i-1,j))
#  endif
!
!  Compute radiation stresses at U-points.
!
            cff4=0.25_r8*FCSr*(FCCr-FSSr)*                              &
     &           (waveE(i,j)+waveE(i-1,j))*                             &
     &           (waven(i,j)+waven(i-1,j))
            cff5=0.5_r8*(owaven(i  ,j)*owaven(i  ,j)+                   &
     &                   owaven(i-1,j)*owaven(i-1,j))
            FCC(i,k)=cff4+0.5_r8*                                       &
     &               waveEr(i,j)*(wavenx(i,j)*wavenx(i,j)+              &
     &                            wavenx(i-1,j)*wavenx(i-1,j))*cff5
            FCS(i,k)=waveEr(i,j)*cff5*0.5_r8*(wavenx(i,j)*waveny(i,j)+  &
     &                                    wavenx(i-1,j)*waveny(i-1,j))
          END DO
        END DO
        DO i=IstrU,Iend
          DO k=1,N(ng)
            cff5=(z_r(i,j,k)-z_r(i-1,j,k))*                             &
     &           on_u(i,j)*(FCC(i,k)-FCC(i,k-1))+                       &
     &           (z_psi(i,j+1)-z_psi(i,j))*                             &
     &           om_u(i,j)*(FCS(i,k)-FCS(i,k-1))
            rustr3d(i,j,k)=rustr3d(i,j,k)-cff5
#  ifdef DIAGNOSTICS_UV
            DiaRU(i,j,k,nrhs,M3hrad)=DiaRU(i,j,k,nrhs,M3hrad)+cff5
#  endif
          END DO
        END DO
        IF (j.ge.JstrV) THEN
          DO i=Istr,Iend
            UFe(i,j)=0.5_r8*(kD(i,j)+kD(i,j-1))
            ocosh(i,j)=1.0_r8/COSH(UFe(i,j))
            osinh(i,j)=1.0_r8/SINH(UFe(i,j))
          END DO
          DO i=Istr,Iend
            DO k=0,N(ng)
              cff1=0.5_r8*(waven(i,j  )*waveE(i,j  )+                   &
     &                     waven(i,j-1)*waveE(i,j-1))
              cff2=1.0_r8+SCALARS(ng)%Cs_w(k)
              FCCr=COSH(UFe(i,j)*cff2)*ocosh(i,j)
              FCSr=COSH(UFe(i,j)*cff2)*osinh(i,j)
              FSSr=SINH(UFe(i,j)*cff2)*osinh(i,j)
              waveEr(i,j)=cff1*FCSr*FCCr
#  ifdef SVENDSEN_ROLLER
              cff3=SCALARS(ng)%Cs_w(k)
              waveEr(i,j)=waveEr(i,j)+
     &                    1.25_r8*(rollA(i,j)+rollA(i,j-1))*            &
     &                    (1.0_r8-cff3**4)/(Dstp(i,j)+Dstp(i,j-1))
#  endif
!
!  Compute radiation stresses at V-points.
!
              cff4=0.25_r8*FCSr*(FCCr-FSSr)*                            &
     &             (waveE(i,j)+waveE(i,j-1))*                           &
     &             (waven(i,j)+waven(i,j-1))
              cff5=0.5_r8*(owaven(i,j  )*owaven(i,j  )+                 &
     &                     owaven(i,j-1)*owaven(i,j-1))
              FCC(i,k)=cff4+0.5_r8*                                     &
     &                 waveEr(i,j)*(waveny(i,j  )*waveny(i,j  )+        &
     &                              waveny(i,j-1)*waveny(i,j-1))*cff5
              FCS(i,k)=waveEr(i,j)*cff5*                                &
     &                 0.5_r8*(wavenx(i  ,j)*waveny(i  ,j)+             &
     &                         wavenx(i-1,j)*waveny(i-1,j))
            END DO
          END DO
          DO i=Istr,Iend
            DO k=1,N(ng)
              cff5=(z_r(i,j,k)-z_r(i-1,j,k))*                           &
     &             on_u(i,j)*(FCS(i,k)-FCS(i,k-1))+                     &
     &             (z_psi(i,j+1)-z_psi(i,j))*                           &
     &             om_u(i,j)*(FCC(i,k)-FCC(i,k-1))
              rvstr3d(i,j,k)=rvstr3d(i,j,k)-cff5
#  ifdef DIAGNOSTICS_UV
              DiaRV(i,j,k,nrhs,M3hrad)=DiaRV(i,j,k,nrhs,M3hrad)+cff5
#  endif
            END DO
          END DO
        END IF
      END DO
# endif
!
!-----------------------------------------------------------------------
!  Determination of vertical stress terms.
!  Notice reuse variable waveEr.
!-----------------------------------------------------------------------
!
!  Component for U-momentum.
!
      J_LOOP : DO j=Jstr,Jend
        DO i=IstrU,Iend
          waveEr(i,j)=0.5_r8*(kD(i,j)+kD(i-1,j))
          ocosh(i,j)=1.0_r8/COSH(waveEr(i,j))
          osinh(i,j)=1.0_r8/SINH(waveEr(i,j))
        END DO
        DO i=IstrU,Iend
          DO k=0,N(ng)
              cff2=1.0_r8+SCALARS(ng)%Cs_w(k)
              FCC(i,k)=COSH(waveEr(i,j)*cff2)*ocosh(i,j)
              FCS(i,k)=COSH(waveEr(i,j)*cff2)*osinh(i,j)
              FSS(i,k)=SINH(waveEr(i,j)*cff2)*osinh(i,j)
          END DO
          cff1=waveE(i,j)-waveE(i-1,j)
          cff3=kD(i,j)-kD(i-1,j)
          cff4=1.0_r8/TANH(waveEr(i,j))
          DO k=0,N(ng)
            cff2=1.0_r8+SCALARS(ng)%Cs_w(k)
            CF(i,k)=0.25_r8*(pm(i,j)+pm(i-1,j))*                        &
     &              (FSS(i,k)*cff1+                                     &
     &               (waveE(i,j)+waveE(i-1,j))*                         &
     &               cff3*(FCS(i,k)*cff2-FSS(i,k)*cff4))
          END DO
        END DO
        DO k=0,N(ng)
          DO i=IstrU,Iend
            FC(i,k)=(FSS(i,k)-FCC(i,k))*CF(i,k)
          END DO
        END DO
        DO k=1,N(ng)
          DO i=IstrU,Iend
            cff=(FC(i,k)-FC(i,k-1))
            Szx(i,j,k)=cff
            rustr3d(i,j,k)=rustr3d(i,j,k)+                              &
     &                     cff*om_u(i,j)*on_u(i,j)
!
!  Convert units to m2/s2 for output purposes.
!
            cff1=0.25_r8*(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))
            rustr3d(i,j,k)=rustr3d(i,j,k)*cff1
# ifdef DIAGNOSTICS_UV
            DiaRU(i,j,k,nrhs,M3vrad)=-cff*om_u(i,j)*on_u(i,j)
# endif
          END DO
        END DO
!
!  Component for V-momentum.
!
        IF (j.ge.JstrV) THEN
          DO i=Istr,Iend
            waveEr(i,j)=0.5_r8*(kD(i,j)+kD(i,j-1))
            ocosh(i,j)=1.0_r8/COSH(waveEr(i,j))
            osinh(i,j)=1.0_r8/SINH(waveEr(i,j))
          END DO
          DO i=Istr,Iend
            DO k=0,N(ng)
              cff2=1.0_r8+SCALARS(ng)%Cs_w(k)
              FCC(i,k)=COSH(waveEr(i,j)*cff2)*ocosh(i,j)
              FCS(i,k)=COSH(waveEr(i,j)*cff2)*osinh(i,j)
              FSS(i,k)=SINH(waveEr(i,j)*cff2)*osinh(i,j)
            END DO
            cff1=waveE(i,j)-waveE(i,j-1)
            cff3=kD(i,j)-kD(i,j-1)
            cff4=1.0_r8/TANH(waveEr(i,j))
            DO k=0,N(ng)
              cff2=1.0_r8+SCALARS(ng)%Cs_w(k)
              CF(i,k)=0.25_r8*(pn(i,j)+pn(i,j-1))*                      &
     &                (FSS(i,k)*cff1+                                   &
     &                 (waveE(i,j)+waveE(i,j-1))*                       &
     &                 cff3*(FCS(i,k)*cff2-FSS(i,k)*cff4))
            END DO
          END DO
          DO k=0,N(ng)
            DO i=Istr,Iend
              FC(i,k)=(FSS(i,k)-FCC(i,k))*CF(i,k)
            END DO
          END DO
          DO k=1,N(ng)
            DO i=Istr,Iend
              cff=(FC(i,k)-FC(i,k-1))
              Szy(i,j,k)=cff
              rvstr3d(i,j,k)=rvstr3d(i,j,k)+cff*om_v(i,j)*on_v(i,j)
!
!  Convert units to m2/s2 for output purposes.
!
              cff1=0.25_r8*(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))
              rvstr3d(i,j,k)=rvstr3d(i,j,k)*cff1
# ifdef DIAGNOSTICS_UV
              DiaRV(i,j,k,nrhs,M3vrad)=-cff*om_v(i,j)*on_v(i,j)
# endif
            END DO
          END DO
        END IF
      END DO J_LOOP
!
!  For a 3D application, compute associated 2D fields by taking the
!  vertical integral of 3D fields.
!
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          cff5=0.5_r8*(Hz(i-1,j,1)+Hz(i,j,1))
          ubar_stokes(i,j)=cff5*u_stokes(i,j,1)
          rustr2d(i,j)=cff5*rustr3d(i,j,1)
          rulag2d(i,j)=cff5*rulag3d(i,j,1)
          DO k=2,N(ng)
            cff5=0.5_r8*(Hz(i-1,j,k)+Hz(i,j,k))
            ubar_stokes(i,j)=ubar_stokes(i,j)+cff5*u_stokes(i,j,k)
            rustr2d(i,j)=rustr2d(i,j)+cff5*rustr3d(i,j,k)
            rulag2d(i,j)=rulag2d(i,j)+cff5*rulag3d(i,j,k)
          END DO
          cff3=0.25_r8*(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))
          cff4=2.0_r8/(Dstp(i-1,j)+Dstp(i,j))
          ubar_stokes(i,j)=ubar_stokes(i,j)*cff4
          rustr2d(i,j)=rustr2d(i,j)*cff3*cff4
          rulag2d(i,j)=rulag2d(i,j)*cff4
        END DO
      END DO
      DO i=Istr,Iend
        DO j=JstrV,Jend
          cff5=0.5_r8*(Hz(i,j-1,1)+Hz(i,j,1))
          vbar_stokes(i,j)=cff5*v_stokes(i,j,1)
          rvstr2d(i,j)=cff5*rvstr3d(i,j,1)
          rvlag2d(i,j)=cff5*rvlag3d(i,j,1)
          DO k=2,N(ng)
            cff5=0.5_r8*(Hz(i,j-1,k)+Hz(i,j,k))
            vbar_stokes(i,j)=vbar_stokes(i,j)+cff5*v_stokes(i,j,k)
            rvstr2d(i,j)=rvstr2d(i,j)+cff5*rvstr3d(i,j,k)
            rvlag2d(i,j)=rvlag2d(i,j)+cff5*rvlag3d(i,j,k)
          END DO
          cff3=0.25_r8*(pm(i,j-1)+pm(i,j))*(pn(i,j-1)+pn(i,j))
          cff4=2.0_r8/(Dstp(i,j-1)+Dstp(i,j))
          vbar_stokes(i,j)=vbar_stokes(i,j)*cff4
          rvstr2d(i,j)=rvstr2d(i,j)*cff3*cff4
          rvlag2d(i,j)=rvlag2d(i,j)*cff4
        END DO
      END DO
#else
!
!   Compute the 2D fields of [Sxx_bar,Syy_bar] and
!   off-diagonal [Sxy_bar,Syx_bar] components of
!   radiation stresses, and 2D stokes veolicities. These terms
!   are only computed from the forcings for non-3D applications.
!   For a 3D application, these 2D fileds are computed as the
!   vertical integrals of the 3D fields (see above).
!
      DO j=Jstr-1,Jend+1
        DO i=Istr-1,Iend+1
!
!  Compute total depth.
!
          Dstp(i,j)=zeta(i,j,1)+h(i,j)
!
!  Compute wave numbers and wave energy.
!
          waven(i,j)=2.0_r8*pi/MAX(Lwave(i,j),Lwave_min)
          owaven(i,j)=1.0_r8/waven(i,j)
          cff=1.5_r8*pi-Dwave(i,j)-angler(i,j)
          wavenx(i,j)=waven(i,j)*DCOS(cff)
          waveny(i,j)=waven(i,j)*DSIN(cff)
!
          waveE(i,j) =0.0625_r8*g*Hwave(i,j)*Hwave(i,j)
!
!  Compute wave celerity and group velocity.
!
          kD(i,j)=MIN(waven(i,j)*Dstp(i,j)+eps,kDmax)
          wavec(i,j)=SQRT(g*owaven(i,j)*TANH(kD(i,j)))
          wavecgoc(i,j)=0.5_r8+                                         &
     &                  kD(i,j)/SINH(2.0_r8*kD(i,j))
!
          cff1=owaven(i,j)*owaven(i,j)
          Sxxl(i,j)=waveE(i,j)*                                         &
     &              (wavecgoc(i,j)*(wavenx(i,j)*wavenx(i,j)*cff1+       &
     &                              1.0_r8)-                            &
     &               0.5_r8)
          Sxyl(i,j)=waveE(i,j)*                                         &
     &              wavecgoc(i,j)*wavenx(i,j)*waveny(i,j)*cff1
          Syyl(i,j)=waveE(i,j)*                                         &
     &              (wavecgoc(i,j)*(waveny(i,j)*waveny(i,j)*cff1+       &
     &                              1.0_r8)-                            &
     &               0.5_r8)

# if defined SVENDSEN_ROLLER
!
!  Add roller contribution.
!
#  ifdef MONO_ROLLER
!
!  Here Wave_break is really wave_area.
!
          cff2=Wave_break(i,j)
#  else
          cff2=0.0424_r8*Hwave(i,j)*Lwave(i,j)*Wave_break(i,j)
#  endif
          cff3=1.0_r8/MAX(Lwave(i,j),Lwave_min)
          waveEr(i,j)=wavec(i,j)*wavec(i,j)*cff2*cff3
!
          Sxxl(i,j)=Sxxl(i,j)+                                          &
     &              waveEr(i,j)*wavenx(i,j)*wavenx(i,j)*cff1
          Sxyl(i,j)=Sxyl(i,j)+                                          &
     &              waveEr(i,j)*wavenx(i,j)*waveny(i,j)*cff1
          Syyl(i,j)=Syyl(i,j)+                                          &
     &              waveEr(i,j)*waveny(i,j)*waveny(i,j)*cff1
# endif
        END DO
      END DO
      DO j=Jstr,Jend+1
        DO i=Istr,Iend+1
          Sxy_psi(i,j)=0.25_r8*(Sxyl(i-1,j-1)+Sxyl(i,j-1)+              &
     &                          Sxyl(i-1,j  )+Sxyl(i,j  ))
        END DO
      END DO
!
!  Compute Radiation Stress contribution to ubar. These stresses
!  have dimensions of m4/s2.
!
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          cff=(Sxxl(i  ,j)-                                             &
     &         Sxxl(i-1,j))*on_u(i,j)+                                  &
     &        (Sxy_psi(i,j+1)-                                          &
     &         Sxy_psi(i,j  ))*om_u(i,j)
          rustr2d(i,j)=cff
!
!  Convert units to m2/s2 for output purposes.
!
          cff=0.25_r8*(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))
          rustr2d(i,j)=rustr2d(i,j)*cff
        END DO
      END DO
!
!  Compute 2D U-stokes velocity.
!
      DO j=JstrR,JendR
        DO i=Istr,IendR
          cff3=2.0_r8/(Dstp(i-1,j)+Dstp(i,j))
          cff5=(wavenx(i-1,j)+wavenx(i,j))/                             &
     &         ((waven(i-1,j)+waven(i,j))*                              &
     &          (wavec(i-1,j)+wavec(i,j)))
!
!  Compute ubar_stokes from wave energy and store old value to compute
!  tendency term.
!
          rulag2d(i,j)=ubar_stokes(i,j)
          ubar_stokes(i,j)=cff3*cff5*(waveE(i-1,j)+waveE(i,j))

# if defined SVENDSEN_ROLLER
#  ifdef MONO_ROLLER
!
!  Here Wave_break is really wave_area.
!
          cff1=0.5_r8*(Wave_break(i-1,j)+Wave_break(i,j))
#  else
          cff1=0.0053033_r8*(Hwave(i-1,j)+Hwave(i,j))*                  &
     &                      (Lwave(i-1,j)+Lwave(i,j))*                  &
     &                      (Wave_break(i-1,j)+Wave_break(i,j))
#  endif
          cff2=4.0_r8/Max(Lwave(i-1,j)+Lwave(i,j),Lwave_min)
          ubar_stokes(i,j)=ubar_stokes(i,j)+g*cff1*cff2*cff5
# endif
!
!  Finalize computation of stokes tendency term.
!
          cff=fac1*om_u(i,j)*on_u(i,j)
          rulag2d(i,j)=(ubar_stokes(i,j)-rulag2d(i,j))*cff
        END DO
      END DO
!
!  Compute Radiation Stress contribution to vbar. These stresses
!  have dimensions of m4/s2.
!
      DO j=JstrV,Jend
        DO i=Istr,Iend
          cff=(Syyl(i,j  )-                                             &
     &         Syyl(i,j-1))*on_v(i,j)+                                  &
     &        (Sxy_psi(i+1,j)-                                          &
     &         Sxy_psi(i  ,j))*on_u(i,j)
          rvstr2d(i,j)=cff
!
!  Convert units to m2/s2 for output purposes.
!
          cff=0.25_r8*(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))
          rvstr2d(i,j)=rvstr2d(i,j)*cff
        END DO
      END DO
!
!  Compute 2D V-stokes velocity.
!
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          cff3=2.0_r8/(Dstp(i,j-1)+Dstp(i,j))
          cff5=(waveny(i,j-1)+waveny(i,j))/                             &
     &         ((waven(i,j-1)+waven(i,j))*                              &
     &          (wavec(i,j-1)+wavec(i,j)))
!
!  Compute vbar_stokes from wave energy and store old value to compute
!  tendency term.
!
          rvlag2d(i,j)=vbar_stokes(i,j)
          vbar_stokes(i,j)=cff3*cff5*(waveE(i,j-1)+waveE(i,j))

# if defined SVENDSEN_ROLLER
#  ifdef MONO_ROLLER
!
!  Here Wave_break is really wave_area.
!
          cff1=0.5_r8*(Wave_break(i,j-1)+Wave_break(i,j))
#  else
          cff1=0.0053033_r8*(Hwave(i,j-1)+Hwave(i,j))*                  &
     &                      (Lwave(i,j-1)+Lwave(i,j))*                  &
     &                      (Wave_break(i,j-1)+Wave_break(i,j))
#  endif
          cff2=4.0_r8/Max(Lwave(i,j-1)+Lwave(i,j),Lwave_min)
          vbar_stokes(i,j)=vbar_stokes(i,j)+g*cff1*cff2*cff5
# endif
!
!  Finalize computation of stokes tendency term.
!
          cff=fac1*om_v(i,j)*on_v(i,j)
          rvlag2d(i,j)=(vbar_stokes(i,j)-rvlag2d(i,j))*cff
        END DO
      END DO
!
!  Fill global arrays with local values.
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          Sxx_bar(i,j)=Sxxl(i,j)
          Sxy_bar(i,j)=Sxyl(i,j)
          Syy_bar(i,j)=Syyl(i,j)
        END DO
      END DO
#endif
!
!  Apply boundary conditions.
!
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Sxx_bar)
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Sxy_bar)
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Syy_bar)
      CALL bc_u2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  ubar_stokes)
      CALL bc_v2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  vbar_stokes)
      CALL bc_u2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  rustr2d)
      CALL bc_v2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  rvstr2d)
#ifdef SOLVE3D
      CALL bc_r3d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, 1, N(ng),                   &
     &                  Sxx)
      CALL bc_r3d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, 1, N(ng),                   &
     &                  Sxy)
      CALL bc_r3d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, 1, N(ng),                   &
     &                  Syy)
      CALL bc_r3d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, 1, N(ng),                   &
     &                  Szx)
      CALL bc_r3d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, 1, N(ng),                   &
     &                  Szy)
      CALL bc_u3d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, 1, N(ng),                   &
     &                  u_stokes)
      CALL bc_v3d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, 1, N(ng),                   &
     &                  v_stokes)
      CALL bc_u3d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, 1, N(ng),                   &
     &                  rustr3d)
      CALL bc_v3d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, 1, N(ng),                   &
     &                  rvstr3d)
#endif
#ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 4,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    rustr2d, rvstr2d, rulag2d, rvlag2d)
      CALL mp_exchange2d (ng, tile, iNLM, 3,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Sxx_bar, Sxy_bar, Syy_bar)
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    ubar_stokes, vbar_stokes)
# ifdef SOLVE3D
      CALL mp_exchange3d (ng, tile, iNLM, 4,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    rustr3d,  rvstr3d, rulag3d,  rvlag3d)
      CALL mp_exchange3d (ng, tile, iNLM, 3,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Sxx,  Sxy, Syy)
      CALL mp_exchange3d (ng, tile, iNLM, 4,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Szx, Szy,                                     &
     &                    u_stokes, v_stokes)
# endif
#endif

      RETURN
      END SUBROUTINE radiation_stress_tile
