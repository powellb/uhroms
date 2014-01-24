!
!svn $Id$
!================================================ Brian Powell, 2014 ===
!  Copyright (c) 2002-2011 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Parameters for Enterococcus ecosystem model:                        !
!                                                                      !
!  AttSW     Light attenuation due to sea water, [1/m].                !
!  BioIter   Maximum number of iterations to achieve convergence of    !
!              the nonlinear solution.                                 !
!  BioIni    Initial concentration for analytical initial (uniform)    !
!              conditions.                                             !
!  Ent_Att   Decay of Enterococcus due to UV.                          !
!  wEntero   Enterococcus sinking rate, [m/day].                       !
!                                                                      !
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
      integer, dimension(Ngrids) :: BioIter

      real(r8), allocatable :: BioIni(:,:)
      real(r8), dimension(Ngrids) :: w16O, w18O      ! m/day
      
#ifdef TANGENT
      real(r8), dimension(Ngrids) :: tl_w16O, tl_w18O
#endif
#ifdef ADJOINT
      real(r8), dimension(Ngrids) :: ad_w16O, ad_w18O
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
