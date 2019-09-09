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
!
!***********************************************************************
      SUBROUTINE ad_extract_transport_2d (ng, Imin, Imax, Jmin, Jmax,   &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ifield, Mobs, NobsSTR, NobsEND,      &
     &                             Xmin, Xmax, Ymin, Ymax,              &
     &                             time, dt,                            &
     &                             ObsType, ObsVetting,                 &
     &                             Tobs, Xobs, Yobs,                    &
     &                             ad_A,                                &
# ifdef MASKING
     &                             Amask,                               &
# endif
     &                             ad_Aobs)
!***********************************************************************
!

      END SUBROUTINE ad_extract_transport_2d

# ifdef SOLVE3D
!
!***********************************************************************
      SUBROUTINE ad_extract_transport3d (ng, Imin, Imax, Jmin, Jmax,    &
     &                             LBi, UBi, LBj, UBj, LBk, UBk,        &
     &                             ifield, Mobs, NobsSTR, NobsEND,      &
     &                             Xmin, Xmax, Ymin, Ymax,              &
     &                             time, dt,                            &
     &                             ObsType, ObsVetting,                 &
     &                             Tobs, Xobs, Yobs, Zobs,              &
     &                             ad_A, Adepth,                        &
#  ifdef MASKING
     &                             Amask,                               &
#  endif
     &                             ad_Aobs)
!***********************************************************************
!

      END SUBROUTINE ad_extract_transport_3d
# endif
#endif
      END MODULE ad_extract_transport_mod
