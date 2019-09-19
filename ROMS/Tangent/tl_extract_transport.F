#include "cppdefs.h"
      MODULE tl_extract_transport_mod
#if (defined FOUR_DVAR || defined VERIFICATION) && \
    (defined OBSERVATIONS && defined TLM_OBS)
!
!svn $Id: tl_extract_transport.F
!================================================== Brian S. Powell ====
!  Copyright (c) 2002-2019 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!   This module extracts the paths of transport obs and calculates     !
!   the tangent transport for each path. Transport can use the         !
!   barotropic velocity (isTransport2d) or the baroclinic velocity     !
!   (isTransport3d).                                                   !
!                                                                      !
!   This will then extract the appropriate velocity at each point      !
!   along with the height (full-column for 2D or cell height for 3D)   !
!   and compute:                                                       !
!     transport = SUM_i [( u_i * sin(path_angle_i) +                   !
!                          v_i * cos(path_angle_i) ) *                 !
!                        path_length(i) * height(i)]                   !
!                                                                      !
!     tl_transport = SUM_i [( tl_u_i * sin(path_angle_i) +             !
!                             tl_v_i * cos(path_angle_i) ) *           !
!                        path_length(i) * height(i)]                   !
!=======================================================================
!
      USE mod_kinds

      implicit none

      PUBLIC tl_extract_transport2d
# ifdef SOLVE3D
      PUBLIC tl_extract_transport3d
# endif

      CONTAINS
!
!***********************************************************************
      SUBROUTINE tl_extract_transport2d (ng, Imin, Imax, Jmin, Jmax,    &
     &                                 LBi, UBi, LBj, UBj,              &
     &                                 Mstr, Mend, Kout,                &
     &                                 Xmin, Xmax, Ymin, Ymax, Aobs)
!***********************************************************************
!
      USE mod_fourdvar
      USE mod_grid
      USE mod_ncparam
      USE mod_ocean
      USE mod_parallel
      USE mod_scalars
# ifdef DISTRIBUTE
      USE distribute_mod,  ONLY : mp_collect
# endif
      USE extract_obs_mod, ONLY : extract_obs2d
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Imin, Imax, Jmin, Jmax
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Mstr, Mend, Kout
!
      real(r8), intent(in) :: Xmin, Xmax, Ymin, Ymax
!
#  ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: Aobs(:)
#  else
      real(r8), intent(inout) :: Aobs(Mobs)
#  endif
!
!  Local variable declarations.
!
      integer :: ic, iobs, Rstr, Rend, Robs, Ncollect
      real(dp) :: TimeLB, TimeUB

      real(r8), parameter :: IniVal = 0.0_r8, Sv=1.0e-6_r8

      integer, allocatable, dimension(:) :: path_type
      real(r8), allocatable, dimension(:) :: path_time, path_h
      real(r8), allocatable, dimension(:) :: path_vetting
      real(r8), allocatable, dimension(:) :: path_tl_ubar, path_tl_vbar
!
!  Find the transport path observations that are available
!  for the current time.
!
      TimeLB=(time(ng)-0.5_dp*dt(ng))/86400.0_dp
      TimeUB=(time(ng)+0.5_dp*dt(ng))/86400.0_dp

      OBS_THREAD : DO iobs=Mstr,Mend
        OBS_TYPE : IF ((ObsType(iobs).eq.                               &
     &                  ObsState2Type(isTransport2d)).and.              &
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
          path_time=Tobs(iobs)

          allocate ( path_type(Robs) )
          path_type=ObsState2Type(isTransport2d)

          allocate ( path_vetting(Robs) )
          path_vetting=IniVal

          allocate ( path_tl_ubar(Robs) )
          path_tl_ubar=IniVal

          allocate ( path_tl_vbar(Robs) )
          path_tl_vbar=IniVal

          allocate ( path_h(Robs) )
          path_h=IniVal
!
!  Get the ubar-component for transport calculations.
!
          CALL extract_obs2d (ng, 1, Lm(ng)+1, 0, Mm(ng)+1,             &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        ObsState2Type(isTransport2d),             &
     &                        Robs, 1, Robs,                            &
     &                        uXmin(ng)+0.5_r8,                         &
     &                        uXmax(ng)+0.5_r8,                         &
     &                        uYmin(ng), uYmax(ng),                     &
     &                        time(ng), dt(ng), path_type,              &
     &                        path_vetting, path_time,                  &
     &                        PathX(Rstr:Rend)+0.5_r8,                  &
     &                        PathY(Rstr:Rend),                         & 
     &                        OCEAN(ng)%tl_ubar(:,:,Kout),              &
#   ifdef MASKING
     &                        GRID(ng)%umask,                           &
#   endif
     &                        path_tl_ubar)
!
!  Get the vbar-component for transport calculations.
!
          CALL extract_obs2d (ng, 0, Lm(ng)+1, 1, Mm(ng)+1,             &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        ObsState2Type(isTransport2d),             &
     &                        Robs, 1, Robs,                            &
     &                        uXmin(ng), uXmax(ng),                     &
     &                        uYmin(ng)+0.5_r8,                         &
     &                        uYmax(ng)+0.5_r8,                         &
     &                        time(ng), dt(ng), path_type,              &
     &                        path_vetting, path_time,                  &
     &                        PathX(Rstr:Rend),                         &
     &                        PathY(Rstr:Rend)+0.5_r8,                  &
     &                        OCEAN(ng)%tl_vbar(:,:,Kout),              &
#   ifdef MASKING
     &                        GRID(ng)%vmask,                           &
#   endif
     &                        path_tl_vbar)

!
!  Get the depth-component for transport calculations.
!
          CALL extract_obs2d (ng, 0, Lm(ng)+1, 0, Mm(ng)+1,             &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        ObsState2Type(isTransport2d),             &
     &                        Robs, 1, Robs,                            &
     &                        uXmin(ng), uXmax(ng),                     &
     &                        uYmin(ng), uYmax(ng),                     &
     &                        time(ng), dt(ng), path_type,              &
     &                        path_vetting, path_time,                  &
     &                        PathX(Rstr:Rend)+0.5_r8,                  &
     &                        PathY(Rstr:Rend),                         &
     &                        GRID(ng)%h,                               &
#   ifdef MASKING
     &                        GRID(ng)%rmask,                           &
#   endif
     &                        path_h)

# ifdef DISTRIBUTE
!
!  Collect the vectors from the tiles
!
          Ncollect=Robs
          CALL mp_collect (ng, iNLM, Ncollect, IniVal, path_h)
          CALL mp_collect (ng, iNLM, Ncollect, IniVal, path_tl_ubar)
          CALL mp_collect (ng, iNLM, Ncollect, IniVal, path_tl_vbar)
# endif
!
!  Compute NLM transport in the master thread after the fields are
!  collected.
!
          IF (Master) THEN
            Aobs(iobs)=IniVal
            DO ic=1,Robs
              IF (path_vetting(ic).ne.0.0_r8) THEN
                Aobs(iobs)=Aobs(iobs)+                                  &
     &                     PathLength(Rstr-1+ic)*path_h(ic)*Sv*         &
     &                     (path_tl_ubar(ic)*SIN(PathAngle(Rstr-1+ic))+ &
     &                      path_tl_vbar(ic)*COS(PathAngle(Rstr-1+ic)))
              END IF
            END DO
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

          IF (allocated(path_tl_ubar)) THEN
            deallocate (path_tl_ubar)
          END IF

          IF (allocated(path_tl_vbar)) THEN
            deallocate (path_tl_vbar)
          END IF

          IF (allocated(path_h)) THEN
            deallocate (path_h)
          END IF

        END IF OBS_TYPE
      END DO OBS_THREAD
      END SUBROUTINE tl_extract_transport2d

# ifdef SOLVE3D
!
!***********************************************************************
      SUBROUTINE tl_extract_transport3d (ng, tile, Imin, Imax,          &
     &                                 Jmin, Jmax, LBi, UBi, LBj, UBj,  &
     &                                 Mstr, Mend, Nout,                &
     &                                 Xmin, Xmax, Ymin, Ymax, Aobs)
!***********************************************************************
!
      USE mod_fourdvar
      USE mod_grid
      USE mod_ncparam
      USE mod_ocean
      USE mod_parallel
      USE mod_scalars
# ifdef DISTRIBUTE
      USE distribute_mod,  ONLY :  mp_collect
# endif
      USE extract_obs_mod, ONLY : extract_obs3d
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, Imin, Imax, Jmin, Jmax
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Mstr, Mend, Nout
!
      real(r8), intent(in) :: Xmin, Xmax, Ymin, Ymax
!
#  ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: Aobs(:)
#  else
      real(r8), intent(inout) :: Aobs(Mobs)
#  endif
!
!  Local variable declarations.
!
      integer :: i, j, k, ic, iobs, Ncollect, Robs, Rstr, Rend
      real(dp) :: TimeLB, TimeUB

      real(r8), parameter :: IniVal = 0.0_r8, Sv = 1.0e-6_r8

      integer, allocatable, dimension(:) :: path_type
      real(r8), allocatable, dimension(:) :: path_time
      real(r8), allocatable, dimension(:) :: path_vetting, path_hz
      real(r8), allocatable, dimension(:) :: path_tl_u, path_tl_v
#  include "set_bounds.h"
!
!  Find the transport path observations that are available
!  for the current time.
!
      TimeLB=(time(ng)-0.5_dp*dt(ng))/86400.0_dp
      TimeUB=(time(ng)+0.5_dp*dt(ng))/86400.0_dp

      OBS_THREAD : DO iobs=Mstr,Mend
        OBS_TYPE : IF ((ObsType(iobs).eq.                               &
     &                  ObsState2Type(isTransport3d)).and.              &
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
          path_time=Tobs(iobs)

          allocate ( path_type(Robs) )
          path_type=ObsState2Type(isTransport3d)

          allocate ( path_vetting(Robs) )
          path_vetting=IniVal

          allocate ( path_tl_u(Robs) )
          path_tl_u=IniVal

          allocate ( path_tl_v(Robs) )
          path_tl_v=IniVal

          allocate ( path_hz(Robs) )
          path_hz=IniVal
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
     &                        ObsState2Type(isTransport3d),             &
     &                        Robs, 1, Robs,                            &
     &                        uXmin(ng)+0.5_r8,                         &
     &                        uXmax(ng)+0.5_r8,                         &
     &                        uYmin(ng), uYmax(ng),                     &
     &                        time(ng), dt(ng), path_type,              &
     &                        path_vetting, path_time,                  &
     &                        PathX(Rstr:Rend)+0.5_r8,                  &
     &                        PathY(Rstr:Rend), PathZ(Rstr:Rend),       &
     &                        OCEAN(ng)%tl_u(:,:,:,Nout),               &
     &                        GRID(ng)%z_v,                             &
#   ifdef MASKING
     &                        GRID(ng)%umask,                           &
#   endif
     &                        path_tl_u)
!
!  Get the v-component for transport calculations.
!
          DO k=1,N(ng)
            DO j=Jstr-1,Jend+1
              DO i=IstrU-1,Iend+1
                GRID(ng)%z_v(i,j,k)=0.5_r8*(GRID(ng)%z_r(i-1,j,k)+      &
     &                                      GRID(ng)%z_r(i  ,j,k))
              END DO
            END DO
          END DO
          CALL extract_obs3d (ng, 0, Lm(ng)+1, 1, Mm(ng)+1,             &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        ObsState2Type(isTransport3d),             &
     &                        Robs, 1, Robs,                            &
     &                        uXmin(ng), uXmax(ng),                     &
     &                        uYmin(ng)+0.5_r8,                         &
     &                        uYmax(ng)+0.5_r8,                         &
     &                        time(ng), dt(ng), path_type,              &
     &                        path_vetting, path_time,                  &
     &                        PathX(Rstr:Rend),                         &
     &                        PathY(Rstr:Rend)+0.5_r8,                  &
     &                        PathZ(Rstr:Rend),                         &
     &                        OCEAN(ng)%tl_v(:,:,:,Nout),               &
     &                        GRID(ng)%z_v,                             &
#   ifdef MASKING
     &                        GRID(ng)%vmask,                           &
#   endif
     &                        path_tl_v)

!
!  Get the thickness-component for transport calculations.
!
          CALL extract_obs3d (ng, 0, Lm(ng)+1, 0, Mm(ng)+1,             &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        ObsState2Type(isTransport3d),             &
     &                        Robs, 1, Robs,                            &
     &                        uXmin(ng), uXmax(ng),                     &
     &                        uYmin(ng), uYmax(ng),                     &
     &                        time(ng), dt(ng), path_type,              &
     &                        path_vetting, path_time,                  &
     &                        PathX(Rstr:Rend)+0.5_r8,                  &
     &                        PathY(Rstr:Rend), PathZ(Rstr:Rend),       &
     &                        GRID(ng)%Hz,                              &
     &                        GRID(ng)%z_r,                             &
#  ifdef MASKING
     &                        GRID(ng)%rmask,                           &
#  endif
     &                        path_hz)

#  ifdef DISTRIBUTE
!
!  Collect the vectors from the tiles
!
          Ncollect=Robs
          CALL mp_collect (ng, iNLM, Ncollect, IniVal, path_hz)
          CALL mp_collect (ng, iNLM, Ncollect, IniVal, path_tl_u)
          CALL mp_collect (ng, iNLM, Ncollect, IniVal, path_tl_v)
#  endif
!
!  Compute NLM transport in the master thread after the fields are
!  collected.
!
          IF (Master) THEN
            Aobs(iobs)=IniVal
            DO ic=1,Robs
              IF (path_vetting(ic).ne.0.0_r8) THEN
                Aobs(iobs)=Aobs(iobs)+                                  &
     &                     PathLength(Rstr-1+ic)*path_hz(ic)*Sv*        &
     &                     (path_tl_u(ic)*SIN(PathAngle(Rstr-1+ic))+    &
     &                      path_tl_v(ic)*COS(PathAngle(Rstr-1+ic)))
              END IF
            END DO
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

          IF (allocated(path_tl_u)) THEN
            deallocate (path_tl_u)
          END IF

          IF (allocated(path_tl_v)) THEN
            deallocate (path_tl_v)
          END IF

          IF (allocated(path_hz)) THEN
            deallocate (path_hz)
          END IF

        END IF OBS_TYPE
      END DO OBS_THREAD
      END SUBROUTINE tl_extract_transport3d
# endif
#endif
      END MODULE tl_extract_transport_mod
