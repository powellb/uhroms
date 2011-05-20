!
!svn $Id$
!================================================== Hernan G. Arango ===
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
      integer :: iEntero                ! Enterococcus concentration
!
!  Biological parameters.
!
      integer, dimension(Ngrids) :: BioIter

#ifdef ANA_BIOLOGY
      real(r8), allocatable :: BioIni(:,:)
#endif
      real(r8), dimension(Ngrids) :: AttSWUV, AttSWBlue  ! 1/m
      real(r8), dimension(Ngrids) :: Ent_DecayUV, Ent_GrowthBlue ! nmol/day
      real(r8), dimension(Ngrids) :: PARfracUV, PARfracBlue ! nondimensional
      real(r8), dimension(Ngrids) :: wEntero         ! m/day
#ifdef TANGENT
      real(r8), dimension(Ngrids) :: tl_PARfracUV, tl_PARfracBlue
      real(r8), dimension(Ngrids) :: tl_wEntero
#endif
#ifdef ADJOINT
      real(r8), dimension(Ngrids) :: ad_PARfracUV, ad_PARfracBlue
      real(r8), dimension(Ngrids) :: ad_wEntero
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
      NBT=1
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
      iEntero=ic+1

      RETURN
      END SUBROUTINE initialize_biology
