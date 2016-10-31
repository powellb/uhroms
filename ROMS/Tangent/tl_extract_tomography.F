#include "cppdefs.h"
      MODULE tl_extract_tomography_mod
#if defined SOLVE3D && defined TLM_OBS && \ 
   (defined FOUR_DVAR || defined VERIFICATION)
!
!svn $Id$
!================================================== Brian S. Powell ====
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
! This module contains the TLM code for extracting the                 !
! total travel time along an acoustic ray path.                        !
!=======================================================================
!
      USE mod_kinds
      USE mod_ocean
      USE mod_grid
      USE mod_ncparam
      USE mod_scalars
      USE mod_stepping
      USE mod_fourdvar
      USE mod_parallel
# ifdef DISTRIBUTE
      USE distribute_mod, ONLY :  mp_collect
# endif

      implicit none
!
!-----------------------------------------------------------------------
!  Set all of the constants
!-----------------------------------------------------------------------
!
      real(r8), parameter :: A1 = +1.44896e+3_r8
      real(r8), parameter :: A2 = +4.591_r8
      real(r8), parameter :: A3 = +5.304e-2_r8
      real(r8), parameter :: A4 = +2.374e-4_r8
      real(r8), parameter :: A5 = +1.34_r8
      real(r8), parameter :: A6 = +1.630e-2_r8
      real(r8), parameter :: A7 = +1.675e-7_r8
      real(r8), parameter :: A8 = +1.025e-2_r8
      real(r8), parameter :: A9 = +7.139e-13_r8
      real(r8), parameter :: B1 = +4.068063159e+13_r8
      real(r8), parameter :: B2 = +4.040830031e+13_r8

      PUBLIC :: tl_extract_tomography

      CONTAINS
!
!***********************************************************************
      SUBROUTINE tl_extract_tomography (ng, Imin, Imax, Jmin, Jmax,     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Mobs, NobsSTR, NobsEND,                 &
     &                          Xmin, Xmax, Ymin, Ymax,                 &
     &                          time, dt,                               &
     &                          ObsScale, Aobs)
!***********************************************************************
      USE extract_obs_mod, ONLY : extract_obs2d, extract_obs3d
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Imin, Imax, Jmin, Jmax
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Mobs, NobsSTR, NobsEND
!
      real(r8), intent(in) :: Xmin, Xmax, Ymin, Ymax
      real(r8), intent(in) :: time, dt
!
#  ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: ObsScale(:)
      real(r8), intent(inout) :: Aobs(:)
#  else
      real(r8), intent(inout) :: ObsScale(Mobs)
      real(r8), intent(inout) :: Aobs(Mobs)
#  endif
!
!  Local variable declarations.
!
      real(r8), parameter :: IniVal = 0.0_r8

      real(r8), allocatable, dimension(:) :: raytime, rayscale
      real(r8), allocatable, dimension(:) :: rayzgrid
      real(r8), allocatable, dimension(:) :: raytemp, raysalt
      real(r8), allocatable, dimension(:) :: raytl_temp, raytl_salt
#  ifdef SPHERICAL
      real(r8), allocatable, dimension(:) :: raylat
#  endif
      integer, allocatable, dimension(:) :: raytype
      real(r8) :: sound, tl_sound, cff1, cff2, cffc, cffs
      real(r8) :: TimeLB, TimeUB, tl_travel_time
      integer :: iobs, ic, idx, Robs, Rstr, Rend, Ncollect
      logical :: isValid
!
!-----------------------------------------------------------------------
!   Go through all of the travel time observations. If they are
!   valid for the current time, then process the ray
!-----------------------------------------------------------------------
!
      TimeLB=(time-0.5_r8*dt)/86400.0_r8
      TimeUB=(time+0.5_r8*dt)/86400.0_r8
!
      OBS_THREAD : DO iobs=NobsSTR,NobsEND
        OBS_TYPE : IF ((ObsType(iobs).eq.isTravelTime).and.           &
     &      ((TimeLB.le.Tobs(iobs)).and.(Tobs(iobs).lt.TimeUB))) THEN
!
!-----------------------------------------------------------------------
!   Find all of the points for this ray
!-----------------------------------------------------------------------
!
          Rend=0
          Rstr=999999999
          DO ic=1,Nray(ng)
            IF (RayID(ic).EQ.ObsMeta(iobs)) THEN
              Rstr=MIN(Rstr,ic)
              Rend=MAX(Rend,ic)
            END IF
          END DO
          Robs=Rend-Rstr+1
          IF (Robs.gt.Nray(ng).or.Rstr.gt.Rend) THEN
            ObsScale(iobs)=0.0_r8
            CYCLE
          END IF
!
!-----------------------------------------------------------------------
!  Allocate the arrays
!-----------------------------------------------------------------------
!
          allocate ( raytime(Robs) )
          raytime = 1.0_r8
          allocate ( rayscale(Robs) )
          rayscale = IniVal
          allocate ( rayzgrid(Robs) )
          rayzgrid = RayZ(Rstr:Rend)
          allocate ( raytemp(Robs) )
          raytemp = IniVal
          allocate ( raysalt(Robs) )
          raysalt = IniVal
          allocate ( raytl_temp(Robs) )
          raytl_temp = IniVal
          allocate ( raytl_salt(Robs) )
          raytl_salt = IniVal
          allocate ( raytype(Robs) )
          raytype = 1
#  ifdef SPHERICAL
          allocate ( raylat(Robs) )
          raylat = IniVal
!
!-----------------------------------------------------------------------
!  Grab the latitude for computing the Earth flattening correction
!-----------------------------------------------------------------------
!
          CALL extract_obs2d (ng, 0, Lm(ng)+1, 1, Mm(ng)+1,             &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        isTvar(itemp),                            &
     &                        Robs, Rstr, Rend,                         &
     &                        Xmin, Xmax, Ymin, Ymax,                   &
     &                        time, dt,                                 &
     &                        raytype*isTvar(itemp), rayscale,          &
     &                        raytime*Tobs(iobs),                       &
     &                        RayX(Rstr:Rend),                          &
     &                        RayY(Rstr:Rend),                          &
     &                        GRID(ng)%latr,                            &
#   ifdef MASKING
     &                        GRID(ng)%vmask,                           &
#   endif
     &                        raylat)
#  endif
!
!-----------------------------------------------------------------------
!  Grab the temperature along the raypath
!-----------------------------------------------------------------------
!
          CALL extract_obs3d (ng, 0, Lm(ng)+1, 0, Mm(ng)+1,             &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          isTvar(itemp),                          &
     &                          Robs, 1, Robs,                          &
     &                          Xmin, Xmax, Ymin, Ymax,                 &
     &                          time, dt,                               &
     &                          raytype*isTvar(itemp), rayscale,        &
     &                          raytime*Tobs(iobs),                     &
     &                          RayX(Rstr:Rend),                        &
     &                          RayY(Rstr:Rend),                        &
     &                          rayzgrid,                               &
     &                          OCEAN(ng)%t(:,:,:,NOUT,itemp),          &
     &                          GRID(ng)%z_r,                           &
#  ifdef MASKING
     &                          GRID(ng)%rmask,                         &
#  endif
     &                          raytemp)
!
!-----------------------------------------------------------------------
!  Grab the salt along the raypath
!-----------------------------------------------------------------------
!
          CALL extract_obs3d (ng, 0, Lm(ng)+1, 0, Mm(ng)+1,             &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          isTvar(isalt),                          &
     &                          Robs, Rstr, Rend,                       &
     &                          Xmin, Xmax, Ymin, Ymax,                 &
     &                          time, dt,                               &
     &                          raytype*isTvar(isalt), rayscale,        &
     &                          raytime*Tobs(iobs),                     &
     &                          RayX(Rstr:Rend),                        &
     &                          RayY(Rstr:Rend),                        &
     &                          rayzgrid,                               &
     &                          OCEAN(ng)%t(:,:,:,NOUT,isalt),          &
     &                          GRID(ng)%z_r,                           &
#  ifdef MASKING
     &                          GRID(ng)%rmask,                         &
#  endif
     &                          raysalt)
!
!-----------------------------------------------------------------------
!  Grab the TLM temperature along the raypath
!-----------------------------------------------------------------------
!
          CALL extract_obs3d (ng, 0, Lm(ng)+1, 0, Mm(ng)+1,             &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          isTvar(itemp),                          &
     &                          Robs, 1, Robs,                          &
     &                          Xmin, Xmax, Ymin, Ymax,                 &
     &                          time, dt,                               &
     &                          raytype*isTvar(itemp), rayscale,        &
     &                          raytime*Tobs(iobs),                     &
     &                          RayX(Rstr:Rend),                        &
     &                          RayY(Rstr:Rend),                        &
     &                          rayzgrid,                               &
     &                          OCEAN(ng)%tl_t(:,:,:,NOUT,itemp),       &
     &                          GRID(ng)%z_r,                           &
#  ifdef MASKING
     &                          GRID(ng)%rmask,                         &
#  endif
     &                          raytl_temp)
!
!-----------------------------------------------------------------------
!  Grab the TLM salt along the raypath
!-----------------------------------------------------------------------
!
          CALL extract_obs3d (ng, 0, Lm(ng)+1, 0, Mm(ng)+1,             &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          isTvar(isalt),                          &
     &                          Robs, Rstr, Rend,                       &
     &                          Xmin, Xmax, Ymin, Ymax,                 &
     &                          time, dt,                               &
     &                          raytype*isTvar(isalt), rayscale,        &
     &                          raytime*Tobs(iobs),                     &
     &                          RayX(Rstr:Rend),                        &
     &                          RayY(Rstr:Rend),                        &
     &                          rayzgrid,                               &
     &                          OCEAN(ng)%tl_t(:,:,:,NOUT,isalt),       &
     &                          GRID(ng)%z_r,                           &
#  ifdef MASKING
     &                          GRID(ng)%rmask,                         &
#  endif
     &                          raytl_salt)
#  ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Collect all extracted data.
!-----------------------------------------------------------------------
!
          Ncollect=Robs
#   ifdef SPHERICAL
          CALL mp_collect (ng, iTLM, Ncollect, IniVal, raylat)
#   endif
          CALL mp_collect (ng, iTLM, Ncollect, IniVal, raytemp)
          CALL mp_collect (ng, iTLM, Ncollect, IniVal, raysalt)
          CALL mp_collect (ng, iTLM, Ncollect, IniVal, raytl_temp)
          CALL mp_collect (ng, iTLM, Ncollect, IniVal, raytl_salt)
#   ifndef IS4DVAR_SENSITIVITY
          CALL mp_collect (ng, iTLM, Ncollect, IniVal, rayscale)
#   endif
#  endif
!
!-----------------------------------------------------------------------
!  Only the Master needs to integrate the entire ray
!-----------------------------------------------------------------------
!
          MASTER_THREAD : IF (Master) THEN
!
!-----------------------------------------------------------------------
!  Compute the sound along the raypath and the total travel time
!-----------------------------------------------------------------------
!
            tl_travel_time=0.0_r8
            ObsScale(iobs)=1.0_r8
            DO ic=Rstr,Rend
              idx=ic-Rstr+1
#  ifndef IS4DVAR_SENSITIVITY
!
!-----------------------------------------------------------------------
!  If any of the points fail, then this is a bad ray that shouldn't
!  be used
!-----------------------------------------------------------------------
!
              IF (rayscale(idx).EQ.0.0_r8) THEN
                ObsScale(iobs)=0.0_r8
                EXIT
              END IF
#  endif
#  ifdef SPHERICAL
!
!-----------------------------------------------------------------------
!  Compute the Earth flattening for the ray path.
!-----------------------------------------------------------------------
!
              cffc = COS(raylat(idx)*deg2rad)**2
              cffs = SIN(raylat(idx)*deg2rad)**2
              cff1 = B1 * SQRT(cffc + (B1*B1)/(B2*B2)*cffs) /           &
     &             SQRT(B1*cffc + B2*cffs)
              cff1 = 1.0_r8 / cff1
              cff2 = 1.0_r8 + RayZ(ic)*cff1*( 1.0_r8+RayZ(ic)*cff1 )
#  else        
              cff2 = 1.0_r8
#  endif
              tl_sound = cff2 * ( raytl_temp(idx) *                     &
     &                (A2 - 2.0_r8*A3*raytemp(idx) +                    &
     &                  3.0_r8*A4*raytemp(idx)*raytemp(idx) +           &
     &                  A8*(raysalt(idx)-35.0_r8) +                     &
     &                  A9*RayDepth(ic)**3) +                           &
     &                  raytl_salt(idx) * (A5 + A8*raytemp(idx)))
              sound = cff2 * (A1 + A2*raytemp(idx) -                    &
     &                A3*raytemp(idx)*raytemp(idx) +                    &
     &                A4*raytemp(idx)**3 + A5*(raysalt(idx)-35.0_r8) +  &
     &                A6*RayDepth(ic) + A7*RayDepth(ic)*RayDepth(ic) +  &
     &                A8*raytemp(idx)*(raysalt(idx)-35.0_r8) +          &
     &                A9*raytemp(idx)*RayDepth(ic)**3)
              tl_travel_time = tl_travel_time -                         &
     &                tl_sound*RayLength(ic)/(sound*sound)
            END DO
!
!-----------------------------------------------------------------------
!  Free the memory and store the result
!-----------------------------------------------------------------------
!
            deallocate ( raytime )
            deallocate ( rayscale )
            deallocate ( rayzgrid )
            deallocate ( raytemp )
            deallocate ( raysalt )
            deallocate ( raytl_temp )
            deallocate ( raytl_salt )
            deallocate ( raytype )
#  ifdef SPHERICAL
            deallocate ( raylat )
#  endif
            IF (ObsScale(iobs).eq.1.0_r8) THEN
              Aobs(iobs)=tl_travel_time
            END IF
          END IF MASTER_THREAD
        END IF OBS_TYPE
      END DO OBS_THREAD
      RETURN

      END SUBROUTINE tl_extract_tomography
#endif
      END MODULE tl_extract_tomography_mod
