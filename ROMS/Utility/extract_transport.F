#include "cppdefs.h"
#if (defined FOUR_DVAR || defined VERIFICATION) && defined OBSERVATIONS
      SUBROUTINE extract_transport (ng, Imin, Imax, Jmin, Jmax,         &
     &                          LBi, UBi, LBj, UBj, LBk, UBk,           &
     &                          ifield, Mobs, NobsSTR, NobsEND,         &
     &                          Xmin, Xmax, Ymin, Ymax,                 &
     &                          time, dt,                               &
     &                          ObsType, ObsVetting,                    &
     &                          Tobs, Xobs, Yobs, Zobs, Aobs)

!
!svn $Id: extract_transport.F 1804 2019-07-09 20:17:48Z arango $
!================================================== Brian S. Powell ====
!  Copyright (c) 2002-2019 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!   This routine extracts the paths of transport and calculates        !
!   the total transport for each path.                                 !
!                                                                      !
!   The transport path angle is stored in obs_meta, and every grid
!   location that goes into the calculation of the transport must be
!   an "observation", such that if the transport is from the upper 10
!   layers, then every k value in the upper 10 layers for a given i, j
!   location is specified as an observation.
!
!   The observation value for every transport observation should be the
!   total observation
!   Currently, only a single transport path is supported. We need an
!   additional vector to provide unique IDs if we want to support more
!   paths.
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
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Imin, Imax, Jmin, Jmax
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
      integer, intent(in) :: ifield, Mobs, NobsSTR, NobsEND
!
      real(r8), intent(in) :: Xmin, Xmax, Ymin, Ymax
      real(dp), intent(in) :: time, dt
!
#  ifdef ASSUMED_SHAPE
      integer, intent(in) :: ObsType(:)

      real(dp), intent(in) :: Tobs(:)
      real(r8), intent(in) :: Xobs(:)
      real(r8), intent(in) :: Yobs(:)
      real(r8), intent(inout) :: ObsVetting(:)
      real(r8), intent(inout) :: Zobs(:)
      real(r8), intent(inout) :: Aobs(:)
#  else
      integer, intent(in) :: ObsType(Mobs)

      real(dp), intent(in) :: Tobs(Mobs)
      real(r8), intent(in) :: Xobs(Mobs)
      real(r8), intent(in) :: Yobs(Mobs)
      real(r8), intent(inout) :: ObsVetting(Mobs)
      real(r8), intent(inout) :: Zobs(Mobs)
      real(r8), intent(inout) :: Aobs(Mobs)
#  endif
!
!  Local variable declarations.
!
      logical :: Lwrote
      integer :: i, j, k, ic, iobs, Rstr, Rend
      real(dp) :: TimeLB, TimeUB

      real(r8), parameter :: IniVal = 0.0_r8

      real(r8), allocatable, dimension(:) :: path_time
      real(r8), allocatable, dimension(:) :: path_type, path_vetting
      real(r8), allocatable, dimension(:) :: path_u, path_v, path_hz
      real(r8) :: transport
# ifdef BGQC
      real(r8), allocatable, dimension(:) :: uBgErr, vBgErr
# endif

!
!  Find the transport path observations that are available
!  for the current time.
!
      TimeLB=(time-0.5_dp*dt)/86400.0_dp
      TimeUB=(time+0.5_dp*dt)/86400.0_dp

      OBS_THREAD : DO iobs=NobsSTR,NobsEND
        OBS_TYPE : IF ((ObsType(iobs).eq.isTransport).and.             &
     &      ((TimeLB.le.Tobs(iobs)).and.(Tobs(iobs).lt.TimeUB))) THEN

!
!  Find the transport path points for this transport observation
!
          Rend=0
          Rstr=99999999
          DO ic=1,Npath(ng)
            IF (PathID(ic).eq.ObsMeta(iobs)) THEN
              Rstr=MIN(Rstr,ic)
              Rend=MAX(Rend,ic)
            END IF
          END DO
          Robs=Rend-Rstr+1
          IF (Rstr.gt.Rend) THEN
            ObsVetting(iobs)=0.0_r8
            CYCLE
          END IF

!
!  Allocate the arrays of the transport path data
!
          allocate ( path_time(Robs) )
          path_time(ic)=Tobs(iobs)

          allocate ( path_type(Robs) )
          path_type(ic)=ObsState2Type(isTransport)

          allocate ( path_vetting(Robs) )
          path_vetting(ic)=IniVal

          allocate ( path_u(Robs) )
          path_u(ic)=IniVal

          allocate ( path_v(Robs) )
          path_v(ic)=IniVal

          allocate ( path_hz(Robs) )
          path_hz(ic)=IniVal

# ifdef BGQC
          allocate ( uBgErr(Robs) )
          path_hz(ic)=IniVal

          allocate ( vBgErr(Robs) )
          path_hz(ic)=IniVal
# endif

!
!  Get the u-component for transport calculations.
!
          DO k=1,N(ng)
            DO j=Jstr-1,Jend+1
              DO i=IstrU-1,Iend+1
                GRID(ng)%z_v(i,j,k)=0.5_r8*(GRID(ng)%z_r(i-1,j,k)+
     &            &GRID(ng)%z_r(i  ,j,k))
              END DO
            END DO
          END DO
          CALL extract_obs3d (ng, 1, Lm(ng)+1, 0, Mm(ng)+1,             &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        ObsState2Type(isTransport),               &
     &                        Mobs, Mstr, Mend,                         &
     &                        uXmin(ng)+0.5_r8,                         &
     &                        uXmax(ng)+0.5_r8,                         &
     &                        uYmin(ng), uYmax(ng),                     &
     &                        time(ng), dt(ng), path_type,              &
     &                        path_vetting, path_time,                  &
     &                        PathX(Rstr:Rend)+0.5_r8,                  &
     &                        PathY(Rstr:Rend), PathZ(Rstr:Rend),       &
     &                        OCEAN(ng)%u(:,:,:,NOUT),                  &
     &                        GRID(ng)%z_v,                             &
#   ifdef MASKING
     &                        GRID(ng)%umask,                           &
#   endif
     &                        path_u)
#   ifdef BGQC
            CALL extract_obs3d (ng, 1, Lm(ng)+1, 0, Mm(ng)+1,           &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          ObsState2Type(isTransport),             &
     &                          Mobs, Mstr, Mend,                       &
     &                          uXmin(ng)+0.5_r8, rXmax(ng)+0.5_r8,     &
     &                          uYmin(ng), uYmax(ng),                   &
     &                          time(ng), dt(ng),                       &
     &                          path_type, path_vetting,                &
     &                          path_time, PathX(Rstr:Rend)+0.5_r8,     &
     &                          PathY(Rstr:Rend), PathZ(Rstr:Rend),     &
     &                          OCEAN(ng)%e_u(:,:,:,1),                 &
     &                          GRID(ng)%z_v,                           &
#    ifdef MASKING
     &                          GRID(ng)%umask,                         &
#    endif
     &                          uBgErr)
#   endif

!
!  Get the u-component for transport calculations.
!
          DO k=1,N(ng)
            DO j=Jstr-1,Jend+1
              DO i=IstrU-1,Iend+1
                GRID(ng)%z_v(i,j,k)=0.5_r8*(GRID(ng)%z_r(i-1,j,k)+
     &            &GRID(ng)%z_r(i  ,j,k))
              END DO
            END DO
          END DO
          CALL extract_obs3d (ng, 0, Lm(ng)+1, 1, Mm(ng)+1,             &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        ObsState2Type(isTransport),               &
     &                        Mobs, Mstr, Mend,                         &
     &                        uXmin(ng), uXmax(ng),                     &
     &                        uYmin(ng)+0.5_r8,                         &
     &                        uYmax(ng)+0.5_r8,                         &
     &                        time(ng), dt(ng), path_type,              &
     &                        path_vetting, path_time,                  &
     &                        PathX(Rstr:Rend),                         &
     &                        PathY(Rstr:Rend)+0.5_r8,                  &
     &                        PathZ(Rstr:Rend),                         &
     &                        OCEAN(ng)%v(:,:,:,NOUT),                  &
     &                        GRID(ng)%z_v,                             &
#   ifdef MASKING
     &                        GRID(ng)%vmask,                           &
#   endif
     &                        path_v)
#   ifdef BGQC
          CALL extract_obs3d (ng, 0, Lm(ng)+1, 1, Mm(ng)+1,             &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        ObsState2Type(isTransport),               &
     &                        Mobs, Mstr, Mend,                         &
     &                        uXmin(ng), uXmax(ng),                     &
     &                        uYmin(ng)+0.5_r8,                         &
     &                        uYmax(ng)+0.5_r8,                         &
     &                        time(ng), dt(ng), path_type,              &
     &                        path_vetting, path_time,                  &
     &                        PathX(Rstr:Rend),                         &
     &                        PathY(Rstr:Rend)+0.5_r8,                  &
     &                        PathZ(Rstr:Rend),                         &
     &                        OCEAN(ng)%e_v(:,:,:,NOUT),                &
     &                        GRID(ng)%z_v,                             &
#   ifdef MASKING
     &                        GRID(ng)%vmask,                           &
#   endif
     &                        vBgErr)
#   endif

!
!  Get the thickness-component for transport calculations.
!
          CALL extract_obs3d (ng, 0, Lm(ng)+1, 0, Mm(ng)+1,             &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        ObsState2Type(isTransport),               &
     &                        Mobs, Mstr, Mend,                         &
     &                        uXmin(ng), uXmax(ng),                     &
     &                        uYmin(ng), uYmax(ng),                     &
     &                        time(ng), dt(ng), path_type,              &
     &                        path_vetting, path_time,                  &
     &                        PathX(Rstr:Rend)+0.5_r8,                  &
     &                        PathY(Rstr:Rend), PathZ(Rstr:Rend),       &
     &                        GRID(ng)%Hz,                              &
     &                        GRID(ng)%z_r,                             &
#   ifdef MASKING
     &                        GRID(ng)%rmask,                           &
#   endif
     &                        path_hz)

# ifdef DISTRIBUTE
!
!  Collect the vectors from the tiles
!
          Ncollect=Robs
#  ifdef SPHERICAL
          CALL mp_collect (ng, iNLM, Ncollect, IniVal, path_time)
#  endif
          CALL mp_collect (ng, iNLM, Ncollect, IniVal, path_hz)
          CALL mp_collect (ng, iNLM, Ncollect, IniVal, path_u)
          CALL mp_collect (ng, iNLM, Ncollect, IniVal, path_v)
#  ifndef IS4DVAR_SENSITIVITY
          CALL mp_collect (ng, iNLM, Ncollect, IniVal, path_vetting)
#  endif
#  ifdef BGQC
          CALL mp_collect (ng, iNLM, Ncollect, IniVal, uBgErr)
          CALL mp_collect (ng, iNLM, Ncollect, IniVal, vBgErr)
#  endif
# endif
!
!  Compute NLM transport in the master thread after the fields are
!  collected.
!
          IF (Master) THEN
            Aobs(iobs)=IniVal
            ObsVetting(iobs)=0.0_r8
            DO ic=1,Robs
              IF path_vetting(ic).ne.0.0_r8 THEN
                Aobs(iobs)=Aobs(iobs)+                                  &
     &                     PathLength(Rstr-1+ic)*path_hz(ic)*           &           
     &                     (path_u(ic)*SIN(PathAngle(Rstr-1+ic))+       &
     &                      path_v(ic)*COS(PathAngle(Rstr-1+ic))
              END IF
            END IF
#   ifdef BGQC
            BgErr(iobs)=MAX(MAXVAL(uBgErr), MAXVAL(vBgErr))
#   endif
            IF (Aobs(iobs).ne.0.0_r8) THEN
              ObsVetting(iobs)=1.0_r8
            END IF
          END IF
!
!  Free the allocated arrays
!
          IF (allocated(path_time)) THEN
            deallocate (path_time)
          END IF

          IF (allocated(path_type)) THEN
            deallocate (path_type)
          END IF

          IF (allocated(path_vetting)) THEN
            deallocate (path_vetting)
          END IF

          IF (allocated(path_u)) THEN
            deallocate (path_u)
          END IF

          IF (allocated(path_v)) THEN
            deallocate (path_v)
          END IF

          IF (allocated(path_hz)) THEN
            deallocate (path_hz)
          END IF

# ifdef BGQC
          IF (allocated(uBgErr)) THEN
            deallocate (uBgErr)
          END IF

          IF (allocated(vBgErr)) THEN
            deallocate (vBgErr)
          END IF
# endif
        END IF OBS_TYPE
      END DO OBS_THREAD
      END SUBROUTINE extract_transport
#else
      SUBROUTINE extract_transport
      RETURN
      END SUBROUTINE extract_transport
#endif
