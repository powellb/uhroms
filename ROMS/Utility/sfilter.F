#include "cppdefs.h"
      MODULE sfilter_mod
# ifdef HESSIAN_SV
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2012 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  These routines applies a shapiro filter to the TL model fields.     !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC :: sfilter
      PUBLIC :: ad_sfilter

      CONTAINS

      SUBROUTINE sfilter(ng, tile, model, Kinp, Ninp, Napp)
!
!=======================================================================
!                                                                      !
!  Apply an order 2 Shapiro filter to each of the TL model fields.      !
!                                                                      !
!=======================================================================
!
      USE mod_param
#  ifdef SOLVE3D
      USE mod_coupling
#  endif
      USE mod_grid
      USE mod_ocean

      USE bc_2d_mod, ONLY: bc_r2d_tile
      USE bc_2d_mod, ONLY: bc_u2d_tile
      USE bc_2d_mod, ONLY: bc_v2d_tile
      USE bc_3d_mod, ONLY: bc_r3d_tile
      USE bc_3d_mod, ONLY: bc_u3d_tile
      USE bc_3d_mod, ONLY: bc_v3d_tile

      USE shapiro_mod, ONLY: shapiro2d_tile
      USE shapiro_mod, ONLY: shapiro3d_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model, Kinp, Ninp, Napp
!
!  Local variable declarations.
!
      integer :: ii, itrc
#  include "tile.h"
!
#  ifdef PROFILE
      CALL wclock_on (ng, iTLM, 7)
#  endif
      DO ii=1,Napp

      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  OCEAN(ng) % tl_zeta(:,:,Kinp))
      CALL shapiro2d_tile (ng, tile, model,                             &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     OCEAN(ng) % tl_zeta(:,:,Kinp))
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  OCEAN(ng) % tl_zeta(:,:,Kinp))
#  ifdef SOLVE3D
      CALL bc_u3d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, 1, N(ng),                   &
     &                  OCEAN(ng) % tl_u(:,:,:,Ninp))
      CALL shapiro3d_tile (ng, tile, model,                             &
     &                     LBi, UBi, LBj, UBj, 1, N(ng),                &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
#  ifdef MASKING
     &                     GRID(ng) % umask,                            &
#  endif
     &                     OCEAN(ng) % tl_u(:,:,:,Ninp))
      CALL bc_u3d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, 1, N(ng),                   &
     &                  OCEAN(ng) % tl_u(:,:,:,Ninp))
      CALL bc_v3d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, 1, N(ng),                   &
     &                  OCEAN(ng) % tl_v(:,:,:,Ninp))
      CALL shapiro3d_tile (ng, tile, model,                             &
     &                     LBi, UBi, LBj, UBj, 1, N(ng),                &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
#  ifdef MASKING
     &                     GRID(ng) % vmask,                            &
#  endif
     &                     OCEAN(ng) % tl_v(:,:,:,Ninp))
      CALL bc_v3d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, 1, N(ng),                   &
     &                  OCEAN(ng) % tl_v(:,:,:,Ninp))
      DO itrc=1,NT(ng)
        CALL bc_r3d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,1, N(ng),                 &
     &                    OCEAN(ng) % tl_t(:,:,:,Ninp,itrc))
        CALL shapiro3d_tile (ng, tile, model,                           &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       IminS, ImaxS, JminS, JmaxS,                &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       OCEAN(ng) % tl_t(:,:,:,Ninp,itrc))
        CALL bc_r3d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    OCEAN(ng) % tl_t(:,:,:,Ninp,itrc))
      END DO
#  else
      CALL bc_u2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  OCEAN(ng) % tl_ubar(:,:,Kinp))
      CALL shapiro2d_tile (ng, tile, model,                             &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
#  ifdef MASKING
     &                     GRID(ng) % umask,                            &
#  endif
     &                     OCEAN(ng) % tl_ubar(:,:,Kinp))
      CALL bc_u2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  OCEAN(ng) % tl_ubar(:,:,Kinp))
      CALL bc_v2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  OCEAN(ng) % tl_vbar(:,:,Kinp))
      CALL shapiro2d_tile (ng, tile, model,                             &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
#  ifdef MASKING
     &                     GRID(ng) % vmask,                            &
#  endif
     &                     OCEAN(ng) % tl_vbar(:,:,Kinp))
      CALL bc_v2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  OCEAN(ng) % tl_vbar(:,:,Kinp))
#  endif

      END DO

#  ifdef PROFILE
      CALL wclock_off (ng, iTLM, 7)
#  endif

      END SUBROUTINE sfilter

      SUBROUTINE ad_sfilter(ng, tile, model, Kinp, Ninp, Napp)
!
!=======================================================================
!                                                                      !
!  Apply an order 2 Shapiro filter to each of the TL model fields.      !
!                                                                      !
!=======================================================================
!
      USE mod_param
#  ifdef SOLVE3D
      USE mod_coupling
#  endif
      USE mod_grid
      USE mod_ocean

      USE ad_bc_2d_mod, ONLY: ad_bc_r2d_tile
      USE ad_bc_2d_mod, ONLY: ad_bc_u2d_tile
      USE ad_bc_2d_mod, ONLY: ad_bc_v2d_tile
      USE ad_bc_3d_mod, ONLY: ad_bc_r3d_tile
      USE ad_bc_3d_mod, ONLY: ad_bc_u3d_tile
      USE ad_bc_3d_mod, ONLY: ad_bc_v3d_tile

      USE shapiro_mod, ONLY: ad_shapiro2d_tile
      USE shapiro_mod, ONLY: ad_shapiro3d_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model, Kinp, Ninp, Napp
!
!  Local variable declarations.
!
      integer :: ii, itrc
#  include "tile.h"
!
#  ifdef PROFILE
      CALL wclock_on (ng, iADM, 7)
#  endif
      DO ii=1,Napp

      CALL ad_bc_r2d_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     OCEAN(ng) % ad_zeta(:,:,Kinp))
      CALL ad_shapiro2d_tile (ng, tile, model,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
#  ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#  endif
     &                        OCEAN(ng) % ad_zeta(:,:,Kinp))
      CALL ad_bc_r2d_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     OCEAN(ng) % ad_zeta(:,:,Kinp))
#  ifdef SOLVE3D
      CALL ad_bc_u3d_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj, 1, N(ng),                &
     &                     OCEAN(ng) % tl_u(:,:,:,Ninp))
      CALL ad_shapiro3d_tile (ng, tile, model,                          &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        IminS, ImaxS, JminS, JmaxS,               &
#  ifdef MASKING
     &                        GRID(ng) % umask,                         &
#  endif
     &                        OCEAN(ng) % ad_u(:,:,:,Ninp))
      CALL ad_bc_u3d_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj, 1, N(ng),                &
     &                     OCEAN(ng) % ad_u(:,:,:,Ninp))
!
      CALL ad_bc_v3d_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj, 1, N(ng),                &
     &                     OCEAN(ng) % ad_v(:,:,:,Ninp))
      CALL ad_shapiro3d_tile (ng, tile, model,                          &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        IminS, ImaxS, JminS, JmaxS,               &
#  ifdef MASKING
     &                        GRID(ng) % vmask,                         &
#  endif
     &                        OCEAN(ng) % ad_v(:,:,:,Ninp))
      CALL ad_bc_v3d_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj, 1, N(ng),                &
     &                     OCEAN(ng) % ad_v(:,:,:,Ninp))
      DO itrc=1,NT(ng)
        CALL ad_bc_r3d_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       OCEAN(ng) % ad_t(:,:,:,Ninp,itrc))
        CALL ad_shapiro3d_tile (ng, tile, model,                        &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          IminS, ImaxS, JminS, JmaxS,             &
#  ifdef MASKING
     &                          GRID(ng) % rmask,                       &
#  endif
     &                          OCEAN(ng) % ad_t(:,:,:,Ninp,itrc))
        CALL ad_bc_r3d_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       OCEAN(ng) % ad_t(:,:,:,Ninp,itrc))
      END DO
#  else
      CALL ad_bc_u2d_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     OCEAN(ng) % ad_ubar(:,:,Kinp))
      CALL ad_shapiro2d_tile (ng, tile, model,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
#  ifdef MASKING
     &                        GRID(ng) % umask,                         &
#  endif
     &                        OCEAN(ng) % ad_ubar(:,:,Kinp))
      CALL ad_bc_u2d_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     OCEAN(ng) % ad_ubar(:,:,Kinp))
      CALL ad_bc_v2d_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     OCEAN(ng) % ad_vbar(:,:,Kinp))
      CALL ad_shapiro2d_tile (ng, tile, model,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
#  ifdef MASKING
     &                        GRID(ng) % vmask,                         &
#  endif
     &                        OCEAN(ng) % ad_vbar(:,:,Kinp))
      CALL ad_bc_v2d_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     OCEAN(ng) % ad_vbar(:,:,Kinp))
#  endif

      END DO

#  ifdef PROFILE
      CALL wclock_off (ng, iADM, 7)
#  endif

      END SUBROUTINE ad_sfilter

# endif
      END MODULE sfilter_mod
