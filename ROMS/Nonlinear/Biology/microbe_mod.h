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
      integer :: iEntero, iVulA, iVulB  ! Microbial Concentrations
!
!  Biological parameters.
!
      integer, dimension(Ngrids) :: BioIter

      real(r8), allocatable :: BioIni(:,:)
      real(r8), dimension(Ngrids) :: AttSWUV, AttSWBlue  ! 1/m
      real(r8), dimension(Ngrids) :: Ent_uvd, Ent_blug ! nmol/day
      real(r8), dimension(Ngrids) :: VulA_uvd, VulA_blug ! nmol/day
      real(r8), dimension(Ngrids) :: VulB_uvd, VulB_blug ! nmol/day
      real(r8), dimension(Ngrids) :: PARfracUV, PARfracBlue ! nondimensional
      real(r8), dimension(Ngrids) :: wEntero, wVulA, wVulB      ! m/day
      real(r8), dimension(Ngrids) :: zVulA, zVulB               ! nmol/day
      real(r8), allocatable :: VulA_pop(:,:,:), VulB_pop(:,:,:)
      real(r8), allocatable :: zVulA_win(:,:,:), zVulB_win(:,:,:)
      real(r8), allocatable :: zVulA_avg(:,:,:), zVulB_avg(:,:,:)
      real(r8), allocatable :: zVulA_std(:,:,:), zVulB_std(:,:,:)

      integer, dimension(Ngrids) :: nVulA_win, nVulB_win
      integer, dimension(Ngrids) :: nVulA_lag, nVulB_lag
      integer, dimension(Ngrids) :: nVulAWeights, nVulBWeights
      real(r8), dimension(40) :: vulAwght, vulAtemp, vulAsalt
      real(r8), dimension(40) :: vulBwght, vulBtemp, vulBsalt
      
#ifdef TANGENT
      real(r8), dimension(Ngrids) :: tl_PARfracUV, tl_PARfracBlue
      real(r8), dimension(Ngrids) :: tl_wEntero
      real(r8), dimension(Ngrids) :: tl_wVulA, tl_wVulB
      real(r8), dimension(Ngrids) :: tl_zVulA, tl_zVulB
#endif
#ifdef ADJOINT
      real(r8), dimension(Ngrids) :: ad_PARfracUV, ad_PARfracBlue
      real(r8), dimension(Ngrids) :: ad_wEntero
      real(r8), dimension(Ngrids) :: ad_wVulA, ad_wVul
      real(r8), dimension(Ngrids) :: ad_zVulA, ad_zVul
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
      NBT=3
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
      iVulA=ic+2
      iVulB=ic+3

      RETURN
      END SUBROUTINE initialize_biology
