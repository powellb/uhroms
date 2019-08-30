#include "cppdefs.h"
#if (defined FOUR_DVAR || defined VERIFICATION) && defined OBSERVATIONS
      SUBROUTINE extract_transport (ng, tile, model)
!
!svn $Id: extract_transport.F 1804 2019-07-09 20:17:48Z arango $
!================================================== Brian S. Powell ====
!  Copyright (c) 2002-2019 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!   This routine extracts the lines of transport and calculates        !
!   the total transport for each line.                                 !
!                                                                      !
!   The transport line angle is stored in obs_meta, and every grid
!   location that goes into the calculation of the transport must be
!   an "observation", such that if the transport is from the upper 10
!   layers, then every k value in the upper 10 layers for a given i, j
!   location is specified as an observation.
!
!   The observation value for every transport observation should be the
!   total observation
!   Currently, only a single transport line is supported. We need an
!   additional vector to provide unique IDs if we want to support more
!   lines.
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_ocean
      USE mod_scalars
      USE mod_stepping

!
!         Use the u and v radial vectors for the u and v transport
!         vectors.
!
          DO iobs=Mstr,Mend
            uradial(iobs)=IniVal
            vradial(iobs)=IniVal
          END DO

!
!  Get the u-component for transport calculations.
!
          DO k=1,N(ng)
            DO j=Jstr-1,Jend+1
              DO i=IstrU-1,Iend+1
                GRID(ng)%z_v(i,j,k)=0.5_r8*(GRID(ng)%z_r(i-1,j,k)+      &
     &                                      GRID(ng)%z_r(i  ,j,k))
              END DO
            END DO
          END DO
          CALL extract_obs3d (ng, 1, Lm(ng)+1, 0, Mm(ng)+1,             &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        ObsState2Type(isTransport),               &
     &                        Mobs, Mstr, Mend,                         &
     &                        uXmin(ng)+0.5_r8, uXmax(ng)+0.5_r8,       &
     &                        uYmin(ng), uYmax(ng),                     &
     &                        time(ng), dt(ng),                         &
     &                        ObsType,  ObsVetting,                     &
     &                        Tobs, Xobs+0.5_r8, Yobs, Zobs,            &
     &                        OCEAN(ng)%u(:,:,:,NOUT),                  &
     &                        GRID(ng)%z_v,                             &
#   ifdef MASKING
     &                        GRID(ng)%umask,                           &
#   endif
     &                        uradial)
#   ifdef BGQC
            CALL extract_obs3d (ng, 1, Lm(ng)+1, 0, Mm(ng)+1,           &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          ObsState2Type(isTransport),             &
     &                          Mobs, Mstr, Mend,                       &
     &                          uXmin(ng)+0.5_r8, rXmax(ng)+0.5_r8,     &
     &                          uYmin(ng), uYmax(ng),                   &
     &                          time(ng), dt(ng),                       &
     &                          ObsType, ObsVetting,                    &
     &                          Tobs, Xobs+0.5_r8, Yobs, Zobs,          &
     &                          OCEAN(ng)%e_u(:,:,:,1),                 &
     &                          GRID(ng)%z_v,                           &
#    ifdef MASKING
     &                          GRID(ng)%umask,                         &
#    endif
     &                          uBgErr)
#   endif
!
!  Get the v-component for transport calculations
!
          DO k=1,N(ng)
            DO j=JstrV-1,Jend+1
              DO i=Istr-1,Iend+1
                GRID(ng)%z_v(i,j,k)=0.5_r8*(GRID(ng)%z_r(i,j-1,k)+      &
     &                                      GRID(ng)%z_r(i,j  ,k))
              END DO
            END DO
          END DO
          CALL extract_obs3d (ng, 0, Lm(ng)+1, 1, Mm(ng)+1,             &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        ObsState2Type(isTransport),               &
     &                        Mobs, Mstr, Mend,                         &
     &                        vXmin(ng), vXmax(ng),                     &
     &                        vYmin(ng)+0.5_r8, vYmax(ng)+0.5_r8,       &
     &                        time(ng), dt(ng),                         &
     &                        ObsType, ObsVetting,                      &
     &                        Tobs, Xobs, Yobs+0.5_r8, Zobs,            &
     &                        OCEAN(ng)%v(:,:,:,NOUT),                  &
     &                        GRID(ng)%z_v,                             &
#   ifdef MASKING
     &                        GRID(ng)%vmask,                           &
#   endif
     &                        vradial)
#   ifdef BGQC
          CALL extract_obs3d (ng, 0, Lm(ng)+1, 1, Mm(ng)+1,             &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        ObsState2Type(isTransport),               &
     &                        Mobs, Mstr, Mend,                         &
     &                        vXmin(ng), vXmax(ng),                     &
     &                        vYmin(ng), vYmax(ng),                     &
     &                        time(ng), dt(ng),                         &
     &                        ObsType, ObsVetting,                      &
     &                        Tobs, Xobs, Yobs+0.5_r8, Zobs,            &
     &                        OCEAN(ng)%e_v(:,:,:,1),                   &
     &                        GRID(ng)%z_v,                             &
#    ifdef MASKING
     &                        GRID(ng)%vmask,                           &
#    endif
     &                        vBgErr)
#   endif
!
!  Get the pm-component for transport calculations.
!
          CALL extract_obs2d (ng, 0, Lm(ng)+1, 0, Mm(ng)+1,             &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        ObsState2Type(isTransport),               &
     &                        Mobs, Mstr, Mend,                         &
     &                        rXmin(ng), rXmax(ng),                     &
     &                        rYmin(ng), rYmax(ng),                     &
     &                        time(ng), dt(ng),                         &
     &                        ObsType, ObsVetting,                      &
     &                        Tobs, Xobs, Yobs,                         &
     &                        GRID(ng)%pm,                              &
#  ifdef MASKING
     &                        GRID(ng)%rmask,                           &
#  endif
     &                        pm)

!
!  Get the pn-component for transport calculations.
!
          CALL extract_obs2d (ng, 0, Lm(ng)+1, 0, Mm(ng)+1,             &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        ObsState2Type(isTransport),               &
     &                        Mobs, Mstr, Mend,                         &
     &                        rXmin(ng), rXmax(ng),                     &
     &                        rYmin(ng), rYmax(ng),                     &
     &                        time(ng), dt(ng),                         &
     &                        ObsType, ObsVetting,                      &
     &                        Tobs, Xobs, Yobs,                         &
     &                        GRID(ng)%pm,                              &
#  ifdef MASKING
     &                        GRID(ng)%rmask,                           &
#  endif
     &                        pn)

!
!  Get the thickness-component for transport calculations.
!
            CALL extract_obs3d (ng, 0, Lm(ng)+1, 0, Mm(ng)+1,           &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          ObsState2Type(isTransport),             &
     &                          Mobs, Mstr, Mend,                       &
     &                          rXmin(ng), rXmax(ng),                   &
     &                          rYmin(ng), rYmax(ng),                   &
     &                          time(ng), dt(ng),                       &
     &                          ObsType, ObsVetting,                    &
     &                          Tobs, Xobs, Yobs, Zobs,                 &
     &                          GRID(ng)%Hz,                            &
     &                          GRID(ng)%z_r,                           &
#   ifdef MASKING
     &                          GRID(ng)%rmask,                         &
#   endif
     &                          dz)

!
!  Compute NLM transport.
!
          DO iobs=Mstr,Mend
            IF (ObsType(iobs).eq.ObsState2Type(isTransport)) THEN
              NLmodVal(iobs)=uradial(iobs)*SIN(ObsMeta(iobs)*           &
     &                                       dz(iobs)/pn(iobs) +        &
     &                       vradial(iobs)*COS(ObsMeta(iobs)*           &
     &                                       dz(iobs)/pm(iobs)
#   ifdef BGQC
              BgErr(iobs)=MAX(uBgErr(iobs), vBgErr(iobs))
#   endif
            END IF
          END DO

          END SUBROUTINE extract_transport
#else
      SUBROUTINE extract_transport
      RETURN
      END SUBROUTINE extract_transport
#endif
