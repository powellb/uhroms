      SUBROUTINE isotope_stflux (ng, tile, model)
!     
!svn $Id$
!================================================ Brian Powell, 2014 ===
!  Copyright (c) 2002-2015 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Compute the surface fluxes of isotopes                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_forces
      USE mod_ncparam
!
! Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model

#include "tile.h"
!
      CALL isotope_stflux_tile (ng, tile, model,                        &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
#ifdef TL_IOMS
     &                      FORCES(ng) % tl_stflx,                      &
#endif
     &                      FORCES(ng) % stflx)

      RETURN
      END SUBROUTINE isotope_stflux
!
!***********************************************************************
      SUBROUTINE isotope_stflux_tile (ng, tile, model,                  &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
#ifdef TL_IOMS
     &                            tl_stflx,                             &
#endif
     &                            stflx)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_biology
      USE mod_forces
      USE mod_ocean
      USE mod_grid
      USE mod_stepping
      USE shapiro_mod

      USE exchange_2d_mod, ONLY : exchange_r2d_tile
#ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
#ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: stflx(LBi:,LBj:,:)
# ifdef TL_IOMS
      real(r8), intent(inout) :: tl_stflx(LBi:,LBj:,:)
# endif
#else
      real(r8), intent(inout) :: stflx(LBi:UBi,LBj:UBj,NT(ng))
# ifdef TL_IOMS
      real(r8), intent(inout) :: tl_stflx(LBi:UBi,LBj:UBj,NT(ng))
# endif
#endif
!
!  Local variable declarations.
!
      integer :: i, j
      real(r8) :: precip, evap, ustar, z0, rhoair, vmu, sc, reno
      real(r8) :: cff1, cff2, cff3, kn18, alphoce, alphair
      real(r8) :: qsea, qatm, delq, Hlv

#include "set_bounds.h"
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
            IF (FORCES(ng)%Pair(i,j).gt.0) THEN
!
!     Compute saturation vapor pressure at ocean temperature (mb),
!     using Teten formula.
!
               cff1=(1.0007_r8+3.46E-6_r8*FORCES(ng)%Pair(i,j))*        &
     &              6.1121_r8*EXP(17.502_r8*                            &
     &              OCEAN(ng)%t(i,j,N(ng),nstp(ng),itemp) /             &
     &              (240.97_r8+(OCEAN(ng)%t(i,j,N(ng),nstp(ng),itemp))))
!
!     Compute Qsea (kg/kg) from saturation vapor pressure.
!
               qsea=0.62197_r8*(cff1/(FORCES(ng)%Pair(i,j) -            &
     &              0.378_r8*cff1))
!
!     Compute the density of surface air
!
               cff2=cff1*FORCES(ng)%Hair(i,j) / FORCES(ng)%Pair(i,j)
               rhoair=FORCES(ng)%Pair(i,j) * 100.0_r8 /                 &
     &              (blk_Rgas*(FORCES(ng)%Tair(i,j)+273.16_r8))*        &
     &              (1+0.61_r8*cff2)
!
!     Compute the frictional velocity
!
               cff2=0.5_r8 * (FORCES(ng)%sustr(i,j) +                   &
     &              FORCES(ng)%sustr(MIN(i+1,Lm(ng)),j))
               cff3=0.5_r8 * (FORCES(ng)%svstr(i,j) +                   &
     &              FORCES(ng)%svstr(i,MIN(j+1,Mm(ng))))
               cff1=SQRT(cff2*cff2 + cff3*cff3)
               ustar=SQRT(cff1)
!
!     Compute the roughness length equation: u*^2 / (81.1 * g)
!
               z0=ustar**2*0.00125693_r8
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
               cff1=(1.0007_r8+3.46E-6_r8*FORCES(ng)%Pair(i,j))*        &
     &              6.1121_r8*EXP(17.502_r8*FORCES(ng)%Tair(i,j) /      &
     &              (240.97_r8+(FORCES(ng)%Tair(i,j))))
!     
!     Compute Qatm (kg/kg) from saturation vapor pressure.
!     
               qatm=0.62197_r8*(cff1/(FORCES(ng)%Pair(i,j) -            &
     &              0.378_r8*cff1))*FORCES(ng)%Hair(i,j)
!     
!     Compute vapor deficit
!     
               delq=MAX(1.0E-5_r8,ABS(qsea-qatm))
               delq=SIGN(delq, qsea-qatm)
#ifdef EMINUSP
               evap=FORCES(ng)%evap(i,j)
#else
!     
!     Compute total evaporation from salt flux, (E-P) :
!     E = (stflx*(freshwater density)+P)
!     
               evap=(stflx(i,j,isalt)*1000.0_r8) + FORCES(ng)%rain(i,j)
#endif
!     
!     Compute the evaporation for O16
!     
               cff1=OCEAN(ng)%t(i,j,N(ng),nstp(ng),i16O) /              &
     &              (OCEAN(ng)%rho(i,j,N(ng)) + rho0)
               cff3=((qsea*cff1) - qatm * FORCES(ng)%o16frac(i,j)) /    &
     &              (delq)
!     
!     Compute the precipitation for O16
!     
               precip=FORCES(ng)%rain(i,j) *                            &
     &              FORCES(ng)%o16frac(i,j)
!     
!     Set the fluxes for O16
!     
               stflx(i,j,i16O)=precip - (evap*cff3)
#ifdef TL_IOMS
               tl_stflx(i,j,i16O)=0.0_r8
#endif
!     
!     Compute the equilibrium fractionation coefficients
!     
               cff1=FORCES(ng)%Tair(i,j)+273.16_r8
               alphair=EXP(1.137E3_r8/(cff1*cff1) -                     &
     &              0.4156_r8/cff1 - 2.0667E-3_r8 )
               cff1=OCEAN(ng)%t(i,j,N(ng),nstp(ng),itemp)+273.16_r8
               alphoce=EXP(1.137E3_r8/(cff1*cff1) -                     &
     &              0.4156_r8/cff1 - 2.0667E-3_r8 )
!     
!     Compute the evaporation for O18
!     
               cff1=OCEAN(ng)%t(i,j,N(ng),nstp(ng),i18O) /              &
     &              (OCEAN(ng)%rho(i,j,N(ng)) + rho0)
               cff3=((qsea*cff1/alphoce) - (qatm *                      &
     &              FORCES(ng)%o18frac(i,j)/alphair)) /                 &
     &              (delq)
!     
!     Compute the precipitation for O18
!     
               precip=FORCES(ng)%rain(i,j) *                            &
     &              FORCES(ng)%o18frac(i,j)
!     
!     Set the fluxes for O18
!     
               stflx(i,j,i18O)=precip-(evap * (1.0_r8 - kn18) * cff3)
#ifdef TL_IOMS
               tl_stflx(i,j,i18O)=0.0_r8
#endif
            END IF
         END DO
      END DO
!
!     Apply the Shapiro Filter to the fluxes of O16 and O18
!
      CALL shapiro2d_tile(ng, tile, model, LBi, UBi, LBj, UBj,          &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      GRID(ng)%rmask, stflx(:,:,i16O))
      CALL shapiro2d_tile(ng, tile, model, LBi, UBi, LBj, UBj,          &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      GRID(ng)%rmask, stflx(:,:,i18O))
!
!  Exchange boundary data.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          stflx(:,:,itrc))
#ifdef TL_IOMS
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          tl_stflx(:,:,itrc))
#endif
      END IF

#ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    stflx(:,:,itrc))
# ifdef TL_IOMS
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_stflx(:,:,itrc))
# endif
#endif
      RETURN
      END SUBROUTINE isotope_stflux_tile
