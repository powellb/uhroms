!
!svn $Id$
!================================================ Brian Powell, 2014 ===
!  Copyright (c) 2002-2015 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Parameters for Isotope ecosystem model:                             !
!                                                                      !
!  BioIni(i16O)   Initial O16 concentration                            !
!  BioIni(i18O)   Initial O18 concentration                            !
!  w16O           Sinking-Rate of O16 (enforce mixing if needed)       !
!  w18O           Sinking-Rate of O18 (enforce mixing if needed)       !
!=======================================================================
!
      USE mod_param
!
      implicit none
!
!  Set biological tracer identification indices.
!
      integer, allocatable :: idbio(:)  ! Biological tracers
      integer :: i16O, i18O  ! Isotope Concentrations
!
!  Biological parameters.
!
      integer, allocatable :: BioIter(:)

#ifdef ANA_BIOLOGY
      real(r8), allocatable :: BioIni(:,:)
#endif
      real(r8), allocatable :: w16O(:)   ! m/day
      real(r8), allocatable :: w18O(:)   ! m/day

#ifdef TANGENT
      real(r8), allocatable :: tl_w16O(:)
      real(r8), allocatable :: tl_w18O(:)
#endif
#ifdef ADJOINT
      real(r8), allocatable :: ad_w16O(:)
      real(r8), allocatable :: ad_w18O(:)
#endif

      CONTAINS

      SUBROUTINE initialize_biology
!
!=======================================================================
!                                                                      !
!  This routine sets several variables needed by the biology model.    !
!  It allocates and assigns biological tracers indices.                !
!                                                                      !
!=======================================================================
!
!  Local variable declarations
!
      integer :: i, ic
!
!-----------------------------------------------------------------------
!  Set number of biological tracers.
!-----------------------------------------------------------------------
!
      NBT=2
!
!-----------------------------------------------------------------------
!  Initialize tracer identification indices.
!-----------------------------------------------------------------------
!
!  Allocate biological tracer vector.
!
      IF (.not.allocated(idbio)) THEN
        allocate ( idbio(NBT) )
      END IF
      IF (.not.allocated(BioIter)) THEN
        allocate ( BioIter(Ngrids) )
      END IF
      IF (.not.allocated(w16O)) THEN
        allocate ( w16O(Ngrids) )
      END IF
      IF (.not.allocated(w18O)) THEN
        allocate ( w18O(Ngrids) )
      END IF
#ifdef TANGENT
      IF (.not.allocated(tl_w16O)) THEN
        allocate ( tl_w16O(Ngrids) )
      END IF
      IF (.not.allocated(tl_w18O)) THEN
        allocate ( tl_w18O(Ngrids) )
      END IF
#endif
#ifdef ADJOINT
      IF (.not.allocated(ad_w16O)) THEN
        allocate ( ad_w16O(Ngrids) )
      END IF
      IF (.not.allocated(ad_w18O)) THEN
        allocate ( ad_w18O(Ngrids) )
      END IF
#endif
!
!  Set identification indices.
!
      ic=NAT+NPT+NCS+NNS
      DO i=1,NBT
        idbio(i)=ic+i
      END DO
      i16O=ic+1
      i18O=ic+2

      RETURN
      END SUBROUTINE initialize_biology
