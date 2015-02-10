      SUBROUTINE biology (ng,tile)
!
!svn $Id$
!************************************************** Brian Powell ***
!  Copyright (c) 2002-2011 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!***********************************************************************
!                                                                      !
!  Enterococcus Model.                                                 !
!                                                                      !
!***********************************************************************
!
      USE mod_param
      USE mod_forces
      USE mod_grid
      USE mod_ncparam
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
#include "tile.h"
!
!  Set header file name.
!
#ifdef DISTRIBUTE
      IF (Lbiofile(iNLM)) THEN
#else
      IF (Lbiofile(iNLM).and.(tile.eq.0)) THEN
#endif
        Lbiofile(iNLM)=.FALSE.
        BIONAME(iNLM)=__FILE__
      END IF
!
#ifdef PROFILE
      CALL wclock_on (ng, iNLM, 15)
#endif
      CALL biology_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj, N(ng), NT(ng),             &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   nstp(ng), nnew(ng),                            &
#ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#endif
     &                   GRID(ng) % Hz,                                 &
     &                   GRID(ng) % z_r,                                &
     &                   GRID(ng) % z_w,                                &
     &                   FORCES(ng) % srflx,                            &
     &                   OCEAN(ng) % t)

#ifdef PROFILE
      CALL wclock_off (ng, iNLM, 15)
#endif
      RETURN
      END SUBROUTINE biology
!
!-----------------------------------------------------------------------
      SUBROUTINE biology_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj, UBk, UBt,            &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         nstp, nnew,                              &
#ifdef MASKING
     &                         rmask,                                   &
#endif
     &                         Hz, z_r, z_w,                            &
     &                         srflx,                                   &
     &                         t)
!-----------------------------------------------------------------------
!
      USE mod_param
      USE mod_biology
      USE mod_ncparam
      USE mod_scalars
      USE nrutil, ONLY : gasdev
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, UBk, UBt
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp, nnew

#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
# endif
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: srflx(LBi:,LBj:)
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
#else
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:UBk)
      real(r8), intent(in) :: srflx(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,UBk,:,UBt)
#endif
!
!  Local variable declarations.
!
      integer, parameter :: Nsink = 3

      integer :: Iter, i, ibio, isink, itime, itrc, iTrcMax, j, k, ks, ii
      integer :: idxA, idxB, idxAlag, idxBlag

      real(r8), parameter :: MinVal = 1.0e-6_r8

      real(r8) :: Att, ExpAtt, Itop, PARUV, PARBlue
      real(r8) :: cff, cff1, cff2, cff3, cff4, dtdays
      real(r8) :: cffL, cffR, cu, dltL, dltR
      real(r8) :: delt, dels
      real(r8), dimension(N(ng)) :: zMeanA, zMeanB, zStdA, zStdB
      integer, dimension(N(ng)) :: countA, countB
      
      integer, dimension(Nsink) :: idsink
      real(r8), dimension(Nsink) :: Wbio

      integer, dimension(IminS:ImaxS,N(ng)) :: ksource

      real(r8), dimension(IminS:ImaxS) :: PARsurUV, PARsurBlue

      real(r8), dimension(NT(ng),2) :: BioTrc

      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: Bio
      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: Bio_bak

      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FC

      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv2
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv3
      real(r8), dimension(IminS:ImaxS,N(ng)) :: UVLight, BlueLight
      real(r8), dimension(IminS:ImaxS,N(ng)) :: WL
      real(r8), dimension(IminS:ImaxS,N(ng)) :: WR
      real(r8), dimension(IminS:ImaxS,N(ng)) :: bL
      real(r8), dimension(IminS:ImaxS,N(ng)) :: bR
      real(r8), dimension(IminS:ImaxS,N(ng)) :: qc

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Add biological Source/Sink terms.
!-----------------------------------------------------------------------
!
!  Set time-stepping size (days) according to the number of iterations.
!
      dtdays=dt(ng)*sec2day/REAL(BioIter(ng),r8)
!
!  Set vertical sinking indentification vector.
!
      idsink(1)=iEntero                  ! Enterococcus
      idsink(2)=iVulA                    ! Vulnificus A
      idsink(3)=iVulB                    ! Vulnificus B
!
!  Set vertical sinking velocity vector in the same order as the
!  identification vector, IDSINK.
!
      Wbio(1)=wEntero(ng)                ! Enterococcus
      Wbio(2)=wVulA(ng)                  ! Vulnificus A
      Wbio(3)=wVulB(ng)                  ! Vulnificus B

!
!  Set up the growth/mortality statistics for each layer
!
      idxA=MOD(iic(ng)-ntstart(ng),nVulA_lag(ng))+1
      idxAlag=idxA+1
      IF (idxAlag.GT.nVulA_lag(ng)) THEN
        idxAlag=1
      END IF
      idxB=MOD(iic(ng)-ntstart(ng),nVulB_lag(ng))+1
      idxBlag=idxB+1
      IF (idxBlag.GT.nVulB_lag(ng)) THEN
        idxBlag=1
      END IF
      DO k=1,N(ng)
        IF (iic(ng)-ntstart(ng).LT.nVulA_lag(ng)) THEN
          zMeanA(k)=zVulA(ng)*dtdays
          zStdA(k)=zMeanA(k)*0.025_r8
        ELSE
          zMeanA(k)=zVulA_avg(idxAlag,k,ng)
          zStdA(k)=zVulA_std(idxAlag,k,ng)
        END IF
        IF (iic(ng)-ntstart(ng).LT.nVulA_lag(ng)) THEN
          zMeanB(k)=zVulB(ng)*dtdays
          zStdB(k)=zMeanB(k)*0.025_r8
        ELSE
          zMeanB(k)=zVulB_avg(idxBlag,k,ng)
          zStdB(k)=zVulB_std(idxBlag,k,ng)
        END IF
        countA(k)=0
        countB(k)=0
        zVulA_avg(idxA,k,ng)=0.0_r8
        zVulA_std(idxA,k,ng)=0.0_r8
        zVulB_avg(idxB,k,ng)=0.0_r8
        zVulB_std(idxB,k,ng)=0.0_r8
      END DO

!
!  Compute inverse thickness to avoid repeated divisions.
!
      J_LOOP : DO j=Jstr,Jend
        DO k=1,N(ng)
          DO i=Istr,Iend
            Hz_inv(i,k)=1.0_r8/Hz(i,j,k)
          END DO
        END DO
        DO k=1,N(ng)-1
          DO i=Istr,Iend
            Hz_inv2(i,k)=1.0_r8/(Hz(i,j,k)+Hz(i,j,k+1))
          END DO
        END DO
        DO k=2,N(ng)-1
          DO i=Istr,Iend
            Hz_inv3(i,k)=1.0_r8/(Hz(i,j,k-1)+Hz(i,j,k)+Hz(i,j,k+1))
          END DO
        END DO
!
!  Restrict biological tracer to be positive definite. If a negative
!  concentration is detected, set to a minimum concentration
!
        DO k=1,N(ng)
          DO i=Istr,Iend
!
!  At input, all tracers (index nnew) from predictor step have
!  transport units (m Tunits) since we do not have yet the new
!  values for zeta and Hz. These are known after the 2D barotropic
!  time-stepping.
!
            DO itrc=1,NBT
              ibio=idbio(itrc)
              BioTrc(ibio,nstp)=t(i,j,k,nstp,ibio)
              BioTrc(ibio,nnew)=t(i,j,k,nnew,ibio)*Hz_inv(i,k)
            END DO
!
!  Impose positive definite concentrations.
!
            cff2=0.0_r8
            DO itime=1,2
              cff1=0.0_r8
              iTrcMax=idbio(1)
              DO itrc=1,NBT
                ibio=idbio(itrc)
                cff1=cff1+MAX(0.0_r8,MinVal-BioTrc(ibio,itime))
                IF (BioTrc(ibio,itime).gt.BioTrc(iTrcMax,itime)) THEN
                  iTrcMax=ibio
                END IF
                BioTrc(ibio,itime)=MAX(MinVal,BioTrc(ibio,itime))
              END DO
              IF (BioTrc(iTrcMax,itime).gt.cff1) THEN
                BioTrc(iTrcMax,itime)=BioTrc(iTrcMax,itime)-cff1
              END IF
            END DO
!
!  Load biological tracers into local arrays.
!
            DO itrc=1,NBT
              ibio=idbio(itrc)
              Bio_bak(i,k,ibio)=BioTrc(ibio,nstp)
              Bio(i,k,ibio)=BioTrc(ibio,nstp)
            END DO
          END DO
        END DO
!
!  Calculate surface Photosynthetically Available Radiation (PAR).  The
!  net shortwave radiation is scaled back to Watts/m2 and multiplied by
!  the fraction that is photosynthetically available, PARfrac.
!
        DO i=Istr,Iend
#ifdef CONST_PAR
          PARsurUV(i)=158.075_r8
#else
          PARsurUV(i)=PARfracUV(ng)*srflx(i,j)*rho0*Cp
          PARsurBlue(i)=PARfracBlue(ng)*srflx(i,j)*rho0*Cp
#endif
        END DO
!
!=======================================================================
!  Start internal iterations to achieve convergence of the nonlinear
!  backward-implicit solution.
!=======================================================================
!
!  During the iterative procedure a series of fractional time steps are
!  performed in a chained mode (splitting by different biological
!  conversion processes) in sequence of the main food chain.  In all
!  stages the concentration of the component being consumed is treated
!  in a fully implicit manner, so the algorithm guarantees non-negative
!  values, no matter how strong the concentration of active consuming
!  component (Phytoplankton or Zooplankton).  The overall algorithm,
!  as well as any stage of it, is formulated in conservative form
!  (except explicit sinking) in sense that the sum of concentration of
!  all components is conserved.
!
!  In the implicit algorithm, we have for example (N: nutrient,
!                                                  P: phytoplankton),
!
!     N(new) = N(old) - uptake * P(old)     uptake = mu * N / (Kn + N)
!                                                    {Michaelis-Menten}
!  below, we set
!                                           The N in the numerator of
!     cff = mu * P(old) / (Kn + N(old))     uptake is treated implicitly
!                                           as N(new)
!
!  so the time-stepping of the equations becomes:
!
!     N(new) = N(old) / (1 + cff)     (1) when substracting a sink term,
!                                         consuming, divide by (1 + cff)
!  and
!
!     P(new) = P(old) + cff * N(new)  (2) when adding a source term,
!                                         growing, add (cff * source)
!
!  Notice that if you substitute (1) in (2), you will get:
!
!     P(new) = P(old) + cff * N(old) / (1 + cff)    (3)
!
!  If you add (1) and (3), you get
!
!     N(new) + P(new) = N(old) + P(old)
!
!  implying conservation regardless how "cff" is computed. Therefore,
!  this scheme is unconditionally stable regardless of the conversion
!  rate. It does not generate negative values since the constituent
!  to be consumed is always treated implicitly. It is also biased
!  toward damping oscillations.
!
!  The iterative loop below is to iterate toward an universal Backward-
!  Euler treatment of all terms. So if there are oscillations in the
!  system, they are only physical oscillations. These iterations,
!  however, do not improve the accuaracy of the solution.
!
        ITER_LOOP: DO Iter=1,BioIter(ng)
!
!  Compute light attenuation as function of depth.
!
          DO i=Istr,Iend
!
!  Make sure it is daytime
!
            PARUV=PARsurUV(i)
            PARBlue=PARsurBlue(i)
            IF (PARsurUV(i).gt.0.0_r8.and.PARsurBlue(i).gt.0.0_r8) THEN
              DO k=N(ng),1,-1
!
!  Compute average light attenuation for each grid cell. Here, AttSW is
!  the light attenuation due to seawater.
!
!  UV first
                Att=AttSWUV(ng)*(z_w(i,j,k)-z_w(i,j,k-1))
                ExpAtt=EXP(-Att)
                Itop=PARUV
                PARUV=Itop*(1.0_r8-ExpAtt)/Att    ! average at cell center
                UVLight(i,k)=PARUV
                PARUV=Itop*ExpAtt
!  Blue second
                Att=AttSWBlue(ng)*(z_w(i,j,k)-z_w(i,j,k-1))
                ExpAtt=EXP(-Att)
                Itop=PARBlue
                PARBlue=Itop*(1.0_r8-ExpAtt)/Att    ! average at cell center
                BlueLight(i,k)=PARBlue
                PARBlue=Itop*ExpAtt
!
!  Light attenuation at the bottom of the grid cell. It is the starting
!  PAR value for the next (deeper) vertical grid cell.
!
              END DO
            ELSE                                       ! night time
              DO k=1,N(ng)
                UVLight(i,k)=0.0_r8
                BlueLight(i,k)=0.0_r8
              END DO
            END IF
          END DO

!
!  Integrate the Microbial model
!
          DO k=1,N(ng)
            DO i=Istr,Iend
!
!  Enterococcus growth to blue-light exposure
!
              cff1=dtdays*BlueLight(i,k)*Ent_GrowthBlue(ng)*            &
     &             Bio(i,k,iEntero)
              Bio(i,k,iEntero)=Bio(i,k,iEntero)+cff1
!
!  Enterococcus mortality to UV exposure
!
              cff1=dtdays*UVLight(i,k)*Ent_DecayUV(ng)*                 &
     &             Bio(i,k,iEntero)
              Bio(i,k,iEntero)=Bio(i,k,iEntero)/(1.0_r8+cff1)
!
!  Vibrio Vulnificus A growth. First, compute the growth rate
!
              cff1=0.0_r8
              DO ii=1,nVulAWeights(ng)
                delt=vulAtemp(ii) - t(i,j,k,nstp,itemp)
                dels=vulAsalt(ii) - t(i,j,k,nstp,isalt)
                cff1=cff1 + vulAwght(ii)*                               &
     &              SQRT(0.25*(delt*delt+dels*dels)+1)
              END DO
!
!  Apply the growth-rate
!
              cff1=dtdays*cff1
              Bio(i,k,iVulA)=Bio(i,k,iVulA)+cff1*Bio(i,k,iVulA)
!
! Build the statistics of the growth-rates relative to the decay
!
              IF (Bio(i,k,iVulA).GT.0.0_r8) THEN
                countA(k) = countA(k) + 1
                zVulA_avg(idxA,k,ng)=zVulA_avg(idxA,k,ng)+cff1
                zVulA_std(idxA,k,ng)=zVulA_std(idxA,k,ng)+              &
     &              (cff1-zMeanA(k))*(cff1-zMeanA(k))
              END IF
!
!  Vibrio Vulnificus B growth. First, compute the growth rate
!
              cff1=0.0_r8
              DO ii=1,nVulBWeights(ng)
                delt=vulAtemp(ii) - t(i,j,k,nstp,itemp)
                dels=vulAsalt(ii) - t(i,j,k,nstp,isalt)
                cff1=cff1 + vulBwght(ii)*                               &
     &              SQRT(0.25*(delt*delt+dels*dels)+1)
              END DO
!
!  Apply the growth-rate
!
              cff1=dtdays*cff1
              Bio(i,k,iVulB)=Bio(i,k,iVulB)+cff1*Bio(i,k,iVulB)
!
! Build the statistics of the growth-rates relative to the decay
!
              IF (Bio(i,k,iVulB).GT.0.0_r8) THEN
                countB = countB + 1
                zVulB_avg(idxB,k,ng)=zVulB_avg(idxB,k,ng)+cff1
                zVulB_std(idxB,k,ng)=zVulB_std(idxB,k,ng)+              &
     &              (cff1-zMeanB(k))*(cff1-zMeanB(k))
              END IF
!
!  Vibrio Vulnificus A mortality.
!
              CALL gasdev(cff3)
              cff2=ABS(zMeanA(k) + zStdA(k)*cff3)
              cff1=1.0_r8+cff2
              Bio(i,k,iVulA)=Bio(i,k,iVulA)/cff1

!
!  Vibrio Vulnificus B mortality.
!
              CALL gasdev(cff3)
              cff2=ABS(zMeanB(k) + zStdB(k)*cff3)
              cff1=1.0_r8+cff2
              Bio(i,k,iVulB)=Bio(i,k,iVulB)/cff1
            END DO
          END DO

!
!-----------------------------------------------------------------------
!  Vertical sinking terms: Enterococcus
!-----------------------------------------------------------------------
!
!  Reconstruct vertical profile of selected biological constituents
!  "Bio(:,:,isink)" in terms of a set of parabolic segments within each
!  grid box. Then, compute semi-Lagrangian flux due to sinking.
!
          SINK_LOOP: DO isink=1,Nsink
            ibio=idsink(isink)
!
!  Copy concentration of biological particulates into scratch array
!  "qc" (q-central, restrict it to be positive) which is hereafter
!  interpreted as a set of grid-box averaged values for biogeochemical
!  constituent concentration.
!
            DO k=1,N(ng)
              DO i=Istr,Iend
                qc(i,k)=Bio(i,k,ibio)
              END DO
            END DO
!
            DO k=N(ng)-1,1,-1
              DO i=Istr,Iend
                FC(i,k)=(qc(i,k+1)-qc(i,k))*Hz_inv2(i,k)
              END DO
            END DO
            DO k=2,N(ng)-1
              DO i=Istr,Iend
                dltR=Hz(i,j,k)*FC(i,k)
                dltL=Hz(i,j,k)*FC(i,k-1)
                cff=Hz(i,j,k-1)+2.0_r8*Hz(i,j,k)+Hz(i,j,k+1)
                cffR=cff*FC(i,k)
                cffL=cff*FC(i,k-1)
!
!  Apply PPM monotonicity constraint to prevent oscillations within the
!  grid box.
!
                IF ((dltR*dltL).le.0.0_r8) THEN
                  dltR=0.0_r8
                  dltL=0.0_r8
                ELSE IF (ABS(dltR).gt.ABS(cffL)) THEN
                  dltR=cffL
                ELSE IF (ABS(dltL).gt.ABS(cffR)) THEN
                  dltL=cffR
                END IF
!
!  Compute right and left side values (bR,bL) of parabolic segments
!  within grid box Hz(k); (WR,WL) are measures of quadratic variations.
!
!  NOTE: Although each parabolic segment is monotonic within its grid
!        box, monotonicity of the whole profile is not guaranteed,
!        because bL(k+1)-bR(k) may still have different sign than
!        qc(i,k+1)-qc(i,k).  This possibility is excluded,
!        after bL and bR are reconciled using WENO procedure.
!
                cff=(dltR-dltL)*Hz_inv3(i,k)
                dltR=dltR-cff*Hz(i,j,k+1)
                dltL=dltL+cff*Hz(i,j,k-1)
                bR(i,k)=qc(i,k)+dltR
                bL(i,k)=qc(i,k)-dltL
                WR(i,k)=(2.0_r8*dltR-dltL)**2
                WL(i,k)=(dltR-2.0_r8*dltL)**2
              END DO
            END DO
            cff=1.0E-14_r8
            DO k=2,N(ng)-2
              DO i=Istr,Iend
                dltL=MAX(cff,WL(i,k  ))
                dltR=MAX(cff,WR(i,k+1))
                bR(i,k)=(dltR*bR(i,k)+dltL*bL(i,k+1))/(dltR+dltL)
                bL(i,k+1)=bR(i,k)
              END DO
            END DO
            DO i=Istr,Iend
              FC(i,N(ng))=0.0_r8            ! NO-flux boundary condition
#if defined LINEAR_CONTINUATION
              bL(i,N(ng))=bR(i,N(ng)-1)
              bR(i,N(ng))=2.0_r8*qc(i,N(ng))-bL(i,N(ng))
#elif defined NEUMANN
              bL(i,N(ng))=bR(i,N(ng)-1)
              bR(i,N(ng))=1.5_r8*qc(i,N(ng))-0.5_r8*bL(i,N(ng))
#else
              bR(i,N(ng))=qc(i,N(ng))       ! default strictly monotonic
              bL(i,N(ng))=qc(i,N(ng))       ! conditions
              bR(i,N(ng)-1)=qc(i,N(ng))
#endif
#if defined LINEAR_CONTINUATION
              bR(i,1)=bL(i,2)
              bL(i,1)=2.0_r8*qc(i,1)-bR(i,1)
#elif defined NEUMANN
              bR(i,1)=bL(i,2)
              bL(i,1)=1.5_r8*qc(i,1)-0.5_r8*bR(i,1)
#else
              bL(i,2)=qc(i,1)               ! bottom grid boxes are
              bR(i,1)=qc(i,1)               ! re-assumed to be
              bL(i,1)=qc(i,1)               ! piecewise constant.
#endif
            END DO
!
!  Apply monotonicity constraint again, since the reconciled interfacial
!  values may cause a non-monotonic behavior of the parabolic segments
!  inside the grid box.
!
            DO k=1,N(ng)
              DO i=Istr,Iend
                dltR=bR(i,k)-qc(i,k)
                dltL=qc(i,k)-bL(i,k)
                cffR=2.0_r8*dltR
                cffL=2.0_r8*dltL
                IF ((dltR*dltL).lt.0.0_r8) THEN
                  dltR=0.0_r8
                  dltL=0.0_r8
                ELSE IF (ABS(dltR).gt.ABS(cffL)) THEN
                  dltR=cffL
                ELSE IF (ABS(dltL).gt.ABS(cffR)) THEN
                  dltL=cffR
                END IF
                bR(i,k)=qc(i,k)+dltR
                bL(i,k)=qc(i,k)-dltL
              END DO
            END DO
!
!  After this moment reconstruction is considered complete. The next
!  stage is to compute vertical advective fluxes, FC. It is expected
!  that sinking may occurs relatively fast, the algorithm is designed
!  to be free of CFL criterion, which is achieved by allowing
!  integration bounds for semi-Lagrangian advective flux to use as
!  many grid boxes in upstream direction as necessary.
!
!  In the two code segments below, WL is the z-coordinate of the
!  departure point for grid box interface z_w with the same indices;
!  FC is the finite volume flux; ksource(:,k) is index of vertical
!  grid box which contains the departure point (restricted by N(ng)).
!  During the search: also add in content of whole grid boxes
!  participating in FC.
!
            cff=dtdays*ABS(Wbio(isink))
            DO k=1,N(ng)
              DO i=Istr,Iend
                FC(i,k-1)=0.0_r8
                WL(i,k)=z_w(i,j,k-1)+cff
                WR(i,k)=Hz(i,j,k)*qc(i,k)
                ksource(i,k)=k
              END DO
            END DO
            DO k=1,N(ng)
              DO ks=k,N(ng)-1
                DO i=Istr,Iend
                  IF (WL(i,k).gt.z_w(i,j,ks)) THEN
                    ksource(i,k)=ks+1
                    FC(i,k-1)=FC(i,k-1)+WR(i,ks)
                  END IF
                END DO
              END DO
            END DO
!
!  Finalize computation of flux: add fractional part.
!
            DO k=1,N(ng)
              DO i=Istr,Iend
                ks=ksource(i,k)
                cu=MIN(1.0_r8,(WL(i,k)-z_w(i,j,ks-1))*Hz_inv(i,ks))
                FC(i,k-1)=FC(i,k-1)+                                    &
     &                    Hz(i,j,ks)*cu*                                &
     &                    (bL(i,ks)+                                    &
     &                     cu*(0.5_r8*(bR(i,ks)-bL(i,ks))-              &
     &                         (1.5_r8-cu)*                             &
     &                         (bR(i,ks)+bL(i,ks)-                      &
     &                          2.0_r8*qc(i,ks))))
              END DO
            END DO
            DO k=1,N(ng)
              DO i=Istr,Iend
                Bio(i,k,ibio)=qc(i,k)+(FC(i,k)-FC(i,k-1))*Hz_inv(i,k)
              END DO
            END DO

          END DO SINK_LOOP
        END DO ITER_LOOP
!
!-----------------------------------------------------------------------
!  Update global tracer variables (m Tunits).
!-----------------------------------------------------------------------
!
        DO itrc=1,NBT
          ibio=idbio(itrc)
          DO k=1,N(ng)
            DO i=Istr,Iend
              t(i,j,k,nnew,ibio)=MAX(t(i,j,k,nnew,ibio)+                &
     &                               (Bio(i,k,ibio)-Bio_bak(i,k,ibio))* &
     &                               Hz(i,j,k),                         &
     &                               0.0_r8)
#ifdef TS_MPDATA
              t(i,j,k,3,ibio)=t(i,j,k,nnew,ibio)*Hz_inv(i,k)
#endif
            END DO
          END DO
        END DO

      END DO J_LOOP

!
! Update the mortality statistics
!
      DO k=1,N(ng)
        zVulA_avg(idxA,k,ng)=zVulA_avg(idxA,k,ng)/countA(k)
        zVulA_std(idxA,k,ng)=MIN(0.025_r8*zVulA_avg(idxA,k,ng),         &
     &                  SQRT(zVulA_std(idxA,k,ng)/countA(k))/2.0_r8)
        zVulB_avg(idxB,k,ng)=zVulB_avg(idxB,k,ng)/countB(k)
        zVulB_std(idxB,k,ng)=MIN(0.025_r8*zVulB_avg(idxB,k,ng),         &
     &                  SQRT(zVulB_std(idxB,k,ng)/countB(k))/2.0_r8)
      END DO
      
      RETURN
      END SUBROUTINE biology_tile
