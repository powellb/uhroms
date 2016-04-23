#include "cppdefs.h"
      MODULE isotope_stflux_mod
#if defined ISOTOPE_FLUXES && defined BULK_FLUXES
!     
!svn $Id$
!================================================ Brian Powell, 2014 ===
!  Copyright (c) 2002-2015 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Thus module computes the surface fluxes of isotopes                 !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC  :: isotope_stflux
!
      CONTAINS

!
!***********************************************************************
      SUBROUTINE isotope_stflux (ng, tile, model)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_forces
      USE mod_ocean
      USE mod_stepping
!
! Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model

!
!  Local variable declarations.
!
#include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, model, 6)
# endif
      CALL isotope_stflux_tile (ng, tile, model,                        &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          nstp(ng),                               &
#  ifdef MASKING
     &                          GRID(ng) % rmask,                       &
#  endif
     &                          FORCES(ng) % stflx,                     &
#  ifdef TL_IOMS
     &                          FORCES(ng) % tl_stflx,                  &
#  endif
     &                          FORCES(ng) % Pair,                      &
     &                          FORCES(ng) % Hair,                      &
     &                          FORCES(ng) % Tair,                      &
     &                          FORCES(ng) % sustr,                     &
     &                          FORCES(ng) % svstr,                     &
#  ifdef EMINUSP
     &                          FORCES(ng) % evap,                      &
#  endif
     &                          FORCES(ng) % rain,                      &
     &                          FORCES(ng) % o16frac,                   &
     &                          FORCES(ng) % o18frac,                   &
     &                          OCEAN(ng) % t,                          &
     &                          OCEAN(ng) % rho)
      
# ifdef PROFILE
      CALL wclock_off (ng, model, 6)
# endif
      RETURN
      END SUBROUTINE isotope_stflux

!
!***********************************************************************
      SUBROUTINE isotope_stflux_tile (ng, tile, model,                  &
     &                            LBi, UBi, LBj, UBj, LBk, UBk,         &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            nstp,                                 &
# ifdef MASKING
     &                            rmask,                                &
# endif
     &                            stflx,                                &
# ifdef TL_IOMS
     &                            tl_stflx,                             &
# endif
     &                            Pair,                                 &
     &                            Hair,                                 &
     &                            Tair,                                 &
     &                            sustr,                                &
     &                            svstr,                                &
# ifdef EMINUSP
     &                            evap,                                 &
# endif
     &                            rain,                                 &
     &                            o16frac,                              &
     &                            o18frac,                              &
     &                            t,                                    &
     &                            rho)
!***********************************************************************
!
      USE mod_biology
      USE mod_param
      USE mod_scalars
      USE shapiro_mod
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
!
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#  endif
      real(r8), intent(inout) :: stflx(LBi:,LBj:,:)
#  ifdef TL_IOMS
      real(r8), intent(inout) :: tl_stflx(LBi:,LBj:,:)
#  endif
      real(r8), intent(in) :: Pair(LBi:,LBj:)
      real(r8), intent(in) :: Hair(LBi:,LBj:)
      real(r8), intent(in) :: Tair(LBi:,LBj:)
      real(r8), intent(in) :: sustr(LBi:,LBj:)
      real(r8), intent(in) :: svstr(LBi:,LBj:)
#  ifdef EMINUSP
      real(r8), intent(in) :: evap(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: rain(LBi:,LBj:)
      real(r8), intent(in) :: o16frac(LBi:,LBj:)
      real(r8), intent(in) :: o18frac(LBi:,LBj:)
      real(r8), intent(in) :: t(LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: rho(LBi:,LBj:,:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(inout) :: stflx(LBi:UBi,LBj:UBj,NT(ng))
#  ifdef TL_IOMS
      real(r8), intent(inout) :: tl_stflx(LBi:UBi,LBj:UBj,NT(ng))
#  endif
      real(r8), intent(in) :: Pair(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Hair(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Tair(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: svstr(LBi:UBi,LBj:UBj)
#  ifdef EMINUSP
      real(r8), intent(in) :: evap(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: rain(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: o16frac(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: o18frac(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(in) :: rho(LBi:UBi,LBj:UBj,N(ng))
# endif
!
!  Local variable declarations.
!
      integer :: i, j, itrc
      real(r8) :: precip, pevap, ustar, z0, rhoair, vmu, sc, reno
      real(r8) :: cff1, cff2, cff3, kn18, alphoce, alphair
      real(r8) :: qsea, qatm, delq, Hlv

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute isotopic surface fluxes.
!
!  Compute using model from Merlivat and Jouzel, "Global Climatic
!  Interpretation of the Deuterium-Oxygen 18 Relationship for
!  Precipitation", JGR, vol. 84, 5029--5033, 1979.
!-----------------------------------------------------------------------
!
      DO j=JstrR,JendR
         DO i=IstrR,IendR
            IF (Pair(i,j).gt.0) THEN
!
!     Compute saturation vapor pressure at ocean temperature (mb),
!     using Teten formula.
!
               cff1=(1.0007_r8 + 3.46E-6_r8*Pair(i,j))*                 &
     &              6.1121_r8*EXP(17.502_r8*                            &
     &              t(i,j,UBk,nstp,itemp) /                             &
     &              (240.97_r8 + (t(i,j,UBk,nstp,itemp))))
!
!     Compute Qsea (kg/kg) from saturation vapor pressure.
!
               qsea=0.62197_r8*(cff1/(Pair(i,j) - 0.378_r8*cff1))
!
!     Compute the density of surface air
!
               cff2=cff1*Hair(i,j) / Pair(i,j)
               rhoair=Pair(i,j) * 100.0_r8 /                            &
     &              (blk_Rgas*(Tair(i,j)+273.16_r8))*(1+0.61_r8*cff2)
!
!     Compute the frictional velocity
!
               cff2=0.5_r8 * (sustr(i,j) + sustr(i+1,j))
               cff3=0.5_r8 * (svstr(i,j) + svstr(i,j+1))
               cff1=SQRT(cff2*cff2 + cff3*cff3)
               ustar=SQRT(cff1)
!
!     Compute the roughness length equation: u*^2 / (81.1 * g)
!
               z0=ustar**2 * 0.00125693_r8
!
!     kinematic viscosity, sc, and reynold's number
!
               vmu=1.7E-5_r8/rhoair
               sc=vmu*42372.881_r8
               reno=ustar*z0/vmu
!
!     Compute turbulent resistance based on the Reynold's number
!
               IF (reno.LT.0.13) THEN
                  cff1=0.6666667_r8
                  cff2=(2.5_r8 * LOG(ustar*blk_ZW(ng)/(30.0_r8*vmu))) / &
     &                 (13.6_r8 * sc**cff1)
               ELSE
                  cff1=0.5_r8
                  cff2=(2.5_r8 * LOG(blk_ZW(ng)/z0) - 5.0_r8) /         &
     &                 (7.3_r8*(reno**0.25_r8)*(sc**cff1))
               END IF
!     
!     The kinetic fractionation coefficient
!     
               cff3=1.0323115_r8**cff1
               kn18=(cff3-1.0_r8)/(cff3+cff2)
!     
!     Compute saturation vapor pressure at air temperature (mb),
!     using Teten formula.
!     
               cff1=(1.0007_r8+3.46E-6_r8*Pair(i,j))*                   &
     &              6.1121_r8*EXP(17.502_r8*Tair(i,j) /                 &
     &              (240.97_r8+Tair(i,j)))
!     
!     Compute Qatm (kg/kg) from saturation vapor pressure.
!     
               qatm=0.62197_r8*(cff1/(Pair(i,j) -                       &
     &              0.378_r8*cff1))*Hair(i,j)
!     
!     Compute vapor deficit
!     
               delq=MAX(1.0E-5_r8,ABS(qsea-qatm))
               delq=SIGN(delq, qsea-qatm)
# ifdef EMINUSP
               pevap=evap(i,j)
# else
!     
!     Compute total evaporation from salt flux, (E-P) :
!     E = (stflx*(freshwater density)+P)
!     
               pevap=(stflx(i,j,isalt)*1000.0_r8) + rain(i,j)
# endif
!     
!     Compute the evaporation for O16
!     
               cff1=t(i,j,UBk,nstp,i16O) / (rho(i,j,UBk) + rho0)
               cff3=((qsea*cff1) - qatm * o16frac(i,j)) / delq
!     
!     Compute the precipitation for O16
!     
               precip=rain(i,j) * o16frac(i,j)
!     
!     Set the fluxes for O16
!     
               stflx(i,j,i16O)=precip - (pevap*cff3)
# ifdef TL_IOMS
               tl_stflx(i,j,i16O)=0.0_r8
# endif
!     
!     Compute the equilibrium fractionation coefficients
!     
               cff1=Tair(i,j)+273.16_r8
               alphair=EXP(1.137E3_r8/(cff1*cff1) -                     &
     &              0.4156_r8/cff1 - 2.0667E-3_r8 )
               cff1=t(i,j,UBk,nstp,itemp)+273.16_r8
               alphoce=EXP(1.137E3_r8/(cff1*cff1) -                     &
     &              0.4156_r8/cff1 - 2.0667E-3_r8 )
!     
!     Compute the evaporation for O18
!     
               cff1=t(i,j,UBk,nstp,i18O) / (rho(i,j,UBk) + rho0)
               cff3=((qsea*cff1/alphoce) - (qatm *                      &
     &              o18frac(i,j)/alphair)) / delq
!     
!     Compute the precipitation for O18
!     
               precip=rain(i,j) * o18frac(i,j)
!     
!     Set the fluxes for O18
!     
               stflx(i,j,i18O)=precip-(pevap * (1.0_r8 - kn18) * cff3)
# ifdef TL_IOMS
               tl_stflx(i,j,i18O)=0.0_r8
# endif
            END IF
         END DO
      END DO
!
!     Apply the Shapiro Filter to the fluxes of O16 and O18
!
      CALL shapiro2d_tile(ng, tile, model, LBi, UBi, LBj, UBj,          &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
# ifdef MASKING
     &                      rmask,                                      &
# endif
     &                      stflx(:,:,i16O))
      CALL shapiro2d_tile(ng, tile, model, LBi, UBi, LBj, UBj,          &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
# ifdef MASKING
     &                      rmask,                                      &
# endif
     &                      stflx(:,:,i18O))
!
!  Exchange boundary data.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
         DO itrc=NAT+1,NT(ng)
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              stflx(:,:,itrc))
# ifdef TL_IOMS
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              tl_stflx(:,:,itrc))
# endif
         END DO
      END IF
      
# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    stflx(:,:,i16O),                              &
     &                    stflx(:,:,i18O))
#  ifdef TL_IOMS
      CALL mp_exchange2d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_stflx(:,:,i16O),                           &
     &                    tl_stflx(:,:,i18O))
#  endif
# endif

      RETURN
      END SUBROUTINE isotope_stflux_tile
#endif
      END MODULE isotope_stflux_mod

