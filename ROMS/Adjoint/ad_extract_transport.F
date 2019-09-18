#include "cppdefs.h"
      MODULE ad_extract_transport_mod
#if defined ADJOINT && defined FOUR_DVAR && defined OBSERVATIONS
      USE mod_kinds

      implicit none

      PUBLIC :: ad_extract_transport_2d
# ifdef SOLVE3D
      PUBLIC :: ad_extract_transport_3d
# endif

      CONTAINS
!svn $Id: ad_extract_transport.F
!================================================== Brian S. Powell ====
!  Copyright (c) 2002-2019 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!   This module extracts the paths of transport obs and calculates     !
!   the adjoint transport for each path. Transport can use the         !
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
!                                                                      !
!                                                                      !
!     f_v_i = f_v_i + ad_transport * cos(path_angle_i) *               !
!                     path_length(i) * height(i)                       !
!                                                                      !
!     f_u_i = f_u_i + ad_transport * cos(path_angle_i) *               !
!                     path_length(i) * height(i)                       !
!                                                                      !
!   where,                                                             !
!     ad_transport = ADmodVal                                          !
!=======================================================================
!
      USE mod_kinds

      implicit none

      PUBLIC ad_extract_transport_2d
# ifdef SOLVE3D
      PUBLIC ad_extract_transport_3d
# endif

      CONTAINS
!
!***********************************************************************
      SUBROUTINE ad_extract_transport_2d (ng, Imin, Imax, Jmin, Jmax,   &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Mstr, Mend, Kout,                       &
     &                          Xmin, Xmax, Ymin, Ymax,                 &
     &                          f_ubar, f_vbar, ADmodVal)
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
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
      integer, intent(in) :: ifield, Mstr, Mend, Kout
!
      real(r8), intent(in) :: Xmin, Xmax, Ymin, Ymax
!
#  ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: f_ubar(LBi:,LBj:)
      real(r8), intent(inout) :: f_vbar(LBi:,LBj:)
      real(r8), intent(inout) :: ADmodVal(:)
#  else
      real(r8), intent(inout) :: f_ubar(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: f_vbar(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ADmodVal(Mobs)
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
      real(r8), allocatable, dimension(:) :: path_ad_ubar, path_ad_vbar
      real(r8) :: ad_transport, cff1
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

          allocate ( path_ad_ubar(Robs) )
          path_ad_ubar=IniVal

          allocate ( path_ad_vbar(Robs) )
          path_ad_vbar=IniVal

          allocate ( path_h(Robs) )
          path_h=IniVal
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
!
!  Calculate the adjoint transport
!
          ad_transport=ADmodVal(iobs)
          ADmodVal(iobs)=IniVal
          DO ic=1,Robs
            cff1=PathLength(Rstr-1+ic)*path_h(ic)
            path_ad_ubar(ic) = ad_transport *                           &
     &                          SIN(PathAngle(Rstr-1+ic) * cff1
            path_ad_vbar(ic) = ad_transport *                           &
     &                          COS(PathAngle(Rstr-1+ic) * cff1
          END DO
# ifdef DISTRIBUTE
!
!  Collect the vectors from the tiles
!
          Ncollect=Robs
          CALL mp_collect (ng, iNLM, Ncollect, IniVal, path_h)
          CALL mp_collect (ng, iNLM, Ncollect, IniVal, path_ad_ubar)
          CALL mp_collect (ng, iNLM, Ncollect, IniVal, path_ad_vbar)
# endif
!
!  Set the adjoint ubar-component from the transport
!
          CALL ad_extract_obs2d (ng, 1, Lm(ng)+1, 0, Mm(ng)+1,          &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           ObsState2Type(isTransport2d),          &
     &                           Robs, 1, Robs,                         &
     &                           uXmin(ng)+0.5_r8,                      &
     &                           uXmax(ng)+0.5_r8,                      &
     &                           uYmin(ng), uYmax(ng),                  &
     &                           time(ng), dt(ng), path_type,           &
     &                           path_vetting, path_time,               &
     &                           PathX(Rstr:Rend)+0.5_r8,               &
     &                           PathY(Rstr:Rend),                      & 
     &                           f_ubar,                                &
#   ifdef MASKING
     &                           GRID(ng)%umask,                        &
#   endif
     &                           path_ad_ubar)
!
!  Set the adjoint vbar-component from the transport
!
          CALL ad_extract_obs2d (ng, 0, Lm(ng)+1, 1, Mm(ng)+1,          &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           ObsState2Type(isTransport2d),          &
     &                           Robs, 1, Robs,                         &
     &                           uXmin(ng), uXmax(ng),                  &
     &                           uYmin(ng)+0.5_r8,                      &
     &                           uYmax(ng)+0.5_r8,                      &
     &                           time(ng), dt(ng), path_type,           &
     &                           path_vetting, path_time,               &
     &                           PathX(Rstr:Rend),                      &
     &                           PathY(Rstr:Rend)+0.5_r8,               &
     &                           f_vbar,                                &
#   ifdef MASKING
     &                           GRID(ng)%vmask,                        &
#   endif
     &                        path_ad_vbar)
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

          IF (allocated(path_ad_ubar)) THEN
            deallocate (path_ad_ubar)
          END IF

          IF (allocated(path_ad_vbar)) THEN
            deallocate (path_ad_vbar)
          END IF

          IF (allocated(path_h)) THEN
            deallocate (path_h)
          END IF

        END IF OBS_TYPE
      END DO OBS_THREAD
      END SUBROUTINE ad_extract_transport_2d

# ifdef SOLVE3D
!
!***********************************************************************
      SUBROUTINE ad_extract_transport_3d (ng, tile Imin, Imax, Jmin,    &
     &                                    Jmax, LBi, UBi, LBj, UBj,     &
     &                                    Mstr, Mend, Nout,             &
     &                                    Xmin, Xmax, Ymin, Ymax,       &
     &                                    f_u, f_v, ADmodVal)
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
      real(r8), intent(inout) :: f_u(LBi:,LBj:,:)
      real(r8), intent(inout) :: f_v(LBi:,LBj:,:)
      real(r8), intent(inout) :: ADmodVal(:)
#  else
      real(r8), intent(inout) :: f_u(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: f_v(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: ADmodVal(Mobs)
#  endif
!
!  Local variable declarations.
!
      integer :: i, j, k, ic, iobs, Ncollect, Robs, Rstr, Rend
      real(dp) :: TimeLB, TimeUB

      real(r8), parameter :: IniVal = 0.0_r8, Sv=1.0e-6_r8

      integer, allocatable, dimension(:) :: path_type
      real(r8), allocatable, dimension(:) :: path_time, path_hz
      real(r8), allocatable, dimension(:) :: path_vetting
      real(r8), allocatable, dimension(:) :: path_ad_u, path_ad_v
      real(r8) :: ad_transport, cff1
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

          allocate ( path_ad_u(Robs) )
          path_ad_u=IniVal

          allocate ( path_ad_v(Robs) )
          path_ad_v=IniVal

          allocate ( path_hz(Robs) )
          path_hz=IniVal
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
!
!  Compute NLM transport in the master thread after the fields are
!  collected.
!
          ad_transport=ADmodVal(ic)
          ADmodVal(ic)=IniVal
          DO ic=1,Robs
            cff1=PathLength(Rstr-1+ic)*path_h(ic)
            path_ad_u(ic) = ad_transport *                              &
     &                          SIN(PathAngle(Rstr-1+ic) * cff1
            path_ad_v(ic) = ad_transport *                              &
     &                          COS(PathAngle(Rstr-1+ic) * cff1
          END DO
#  ifdef DISTRIBUTE
!
!  Collect the vectors from the tiles
!
          Ncollect=Robs
          CALL mp_collect (ng, iNLM, Ncollect, IniVal, path_hz)
          CALL mp_collect (ng, iNLM, Ncollect, IniVal, path_ad_u)
          CALL mp_collect (ng, iNLM, Ncollect, IniVal, path_ad_v)
#  endif
!
!  Set the adjoint u-component from the transport
!
          DO k=1,N(ng)
            DO j=Jstr-1,Jend+1
              DO i=IstrU-1,Iend+1
                GRID(ng)%z_v(i,j,k)=0.5_r8*(GRID(ng)%z_r(i-1,j,k)+      &
     &                                      GRID(ng)%z_r(i  ,j,k))
              END DO
            END DO
          END DO
          CALL ad_extract_obs3d (ng, 1, Lm(ng)+1, 0, Mm(ng)+1,          &
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
     &                        f_u,                                      &
     &                        GRID(ng)%z_v,                             &
#   ifdef MASKING
     &                        GRID(ng)%umask,                           &
#   endif
     &                        path_ad_u)
!
!  Set the adjoint v-component from the transport
!
          DO k=1,N(ng)
            DO j=Jstr-1,Jend+1
              DO i=IstrU-1,Iend+1
                GRID(ng)%z_v(i,j,k)=0.5_r8*(GRID(ng)%z_r(i-1,j,k)+      &
     &                                      GRID(ng)%z_r(i  ,j,k))
              END DO
            END DO
          END DO
          CALL ad_extract_obs3d (ng, 0, Lm(ng)+1, 1, Mm(ng)+1,          &
     &                           LBi, UBi, LBj, UBj, 1, N(ng),          &
     &                           ObsState2Type(isTransport3d),          &
     &                           Robs, 1, Robs,                         &
     &                           uXmin(ng), uXmax(ng),                  &
     &                           uYmin(ng)+0.5_r8,                      &
     &                           uYmax(ng)+0.5_r8,                      &
     &                           time(ng), dt(ng), path_type,           &
     &                           path_vetting, path_time,               &
     &                           PathX(Rstr:Rend),                      &
     &                           PathY(Rstr:Rend)+0.5_r8,               &
     &                           PathZ(Rstr:Rend),                      &
     &                           f_v,                                   &
     &                           GRID(ng)%z_v,                          &
#   ifdef MASKING
     &                           GRID(ng)%vmask,                        &
#   endif
     &                           path_ad_v)

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

          IF (allocated(path_ad_u)) THEN
            deallocate (path_ad_u)
          END IF

          IF (allocated(path_ad_v)) THEN
            deallocate (path_ad_v)
          END IF

          IF (allocated(path_hz)) THEN
            deallocate (path_hz)
          END IF

        END IF OBS_TYPE
      END DO OBS_THREAD
      END SUBROUTINE ad_extract_transport_3d
# endif
#endif
      END MODULE ad_extract_transport_mod
