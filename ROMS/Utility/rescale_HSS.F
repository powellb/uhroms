#include "cppdefs.h"

      MODULE rescale_HSS_mod

#if defined WEAK_CONSTRAINT && defined EVOLVED_LCZ
!
!=======================================================================
! This routine reads the dual Lanczos vectors evolves into primal space
! by the adjoint model and normalizes them according to cg_beta.
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC :: rescale_HSS

      CONTAINS
!
!***********************************************************************
      SUBROUTINE rescale_HSS (ng, tile, model, innLoop, outLoop)
!***********************************************************************
!
      USE mod_param
# ifdef ADJUST_BOUNDARY
      USE mod_boundary
# endif
# ifdef SOLVE3D
      USE mod_coupling
# endif
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      USE mod_forces
# endif
      USE mod_grid
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model, innLoop, outLoop
!
!  Local variable declarations.
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, model, 36)
# endif
      CALL rescale_HSS_tile (ng, tile, model,                           &
     &                     LBi, UBi, LBj, UBj, LBij, UBij,              &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     innLoop, outLoop,                            &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
     &                     GRID(ng) % umask,                            &
     &                     GRID(ng) % vmask,                            &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                     BOUNDARY(ng) % ad_t_obc,                     &
     &                     BOUNDARY(ng) % ad_u_obc,                     &
     &                     BOUNDARY(ng) % ad_v_obc,                     &
#  endif
     &                     BOUNDARY(ng) % ad_ubar_obc,                  &
     &                     BOUNDARY(ng) % ad_vbar_obc,                  &
     &                     BOUNDARY(ng) % ad_zeta_obc,                  &
# endif
# ifdef ADJUST_WSTRESS
     &                     FORCES(ng) % ad_ustr,                        &
     &                     FORCES(ng) % ad_vstr,                        &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                     FORCES(ng) % ad_tflux,                       &
#  endif
     &                     OCEAN(ng) % ad_t,                            &
     &                     OCEAN(ng) % ad_u,                            &
     &                     OCEAN(ng) % ad_v,                            &
# else
     &                     OCEAN(ng) % ad_ubar,                         &
     &                     OCEAN(ng) % ad_vbar,                         &
# endif
     &                     OCEAN(ng) % ad_zeta)
# ifdef PROFILE
      CALL wclock_on (ng, model, 36)
# endif
      RETURN
      END SUBROUTINE rescale_HSS
!
!***********************************************************************
      SUBROUTINE rescale_HSS_tile (ng, tile, model,                     &
     &                           LBi, UBi, LBj, UBj, LBij, UBij,        &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           innLoop, outLoop,                      &
# ifdef MASKING
     &                           rmask, umask, vmask,                   &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                           ad_t_obc, ad_u_obc, ad_v_obc,          &
#  endif
     &                           ad_ubar_obc, ad_vbar_obc,              &
     &                           ad_zeta_obc,                           &
# endif
# ifdef ADJUST_WSTRESS
     &                           ad_ustr, ad_vstr,                      &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                           ad_tflux,                              &
#  endif
     &                           ad_t, ad_u, ad_v,                      &
# else
     &                           ad_ubar, ad_vbar,                      &
# endif
     &                           ad_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_iounits
      USE mod_scalars
      USE mod_ncparam
      USE mod_netcdf

# ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcastf, mp_bcasti
# endif
      USE state_scale_mod, ONLY : state_scale
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBij, UBij
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: innLoop, outLoop
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t_obc(LBij:,:,:,:,:,:)
      real(r8), intent(inout) :: ad_u_obc(LBij:,:,:,:,:)
      real(r8), intent(inout) :: ad_v_obc(LBij:,:,:,:,:)
#   endif
      real(r8), intent(inout) :: ad_ubar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: ad_vbar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: ad_zeta_obc(LBij:,:,:,:)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ad_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_vstr(LBi:,LBj:,:,:)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: ad_tflux(LBi:,LBj:,:,:,:)
#   endif
      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_v(LBi:,LBj:,:,:)
#  else
      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
#  endif
      real(r8), intent(inout) :: ad_zeta(LBi:,LBj:,:)

# else

#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t_obc(LBij:UBij,N(ng),4,            &
     &                                    Nbrec(ng),2,NT(ng))
      real(r8), intent(inout) :: ad_u_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
      real(r8), intent(inout) :: ad_v_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
#   endif
      real(r8), intent(inout) :: ad_ubar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: ad_vbar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: ad_zeta_obc(LBij:UBij,4,Nbrec(ng),2)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ad_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: ad_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: ad_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
#   endif
      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
#  else
      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
#  endif
      real(r8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,3)
# endif
!
!  Local variable declarations.
!
      logical :: Ltrans

      integer :: L1 = 1
      integer :: L2 = 2

      integer :: rec, ndefHSS
      integer, save :: ncidsav

      real(r8) :: fac

      character (len=13) :: string
      character (len=256) :: ncname

# include "set_bounds.h"
!
        WRITE (ncname,10) TRIM(HSS(ng)%base), outLoop
 10     FORMAT (a,'_',i3.3,'.nc')

!AMM        ndefHSS=innLoop
        ndefHSS=0
!
        CALL netcdf_put_fvar (ng, model, HSS(ng)%name, 'cg_beta',       &
     &                        cg_beta, (/1,1/), (/Ninner+1,Nouter/),    &
     &                        ncid = HSS(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
!
        CALL netcdf_put_fvar (ng, model, HSS(ng)%name, 'cg_delta',      &
     &                        cg_delta, (/1,1/), (/Ninner,Nouter/),     &
     &                        ncid = HSS(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
!
        DO rec=1,innLoop
          CALL read_state (ng, tile, model,                             &
     &                     LBi, UBi, LBj, UBj, LBij, UBij,              &
     &                     L1, rec,                                     &
     &                     ndefHSS, HSS(ng)%ncid, TRIM(ncname),         &
# ifdef MASKING
     &                     rmask, umask, vmask,                         &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                     ad_t_obc, ad_u_obc, ad_v_obc,                &
#  endif
     &                     ad_ubar_obc, ad_vbar_obc,                    &
     &                     ad_zeta_obc,                                 &
# endif
# ifdef ADJUST_WSTRESS
     &                     ad_ustr, ad_vstr,                            &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                     ad_tflux,                                    &
#  endif
     &                     ad_t, ad_u, ad_v,                            &
# else
     &                     ad_ubar, ad_vbar,                            &
# endif
     &                     ad_zeta)
          IF (exit_flag.ne.NoError) RETURN

      IF (rec.eq.1) THEN
        fac=1.0_r8/cg_Gnorm_v(outLoop)
      ELSE
        fac=1.0_r8/cg_beta(rec,outLoop)
      END IF

      CALL state_scale (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, LBij, UBij,                 &
     &                  L1, L1, fac,                                    &
# ifdef MASKING
     &                  rmask, umask, vmask,                            &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                  ad_t_obc, ad_u_obc, ad_v_obc,                   &
#  endif
     &                  ad_ubar_obc, ad_vbar_obc,                       &
     &                  ad_zeta_obc,                                    &
# endif
# ifdef ADJUST_WSTRESS
     &                  ad_ustr, ad_vstr,                               &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                  ad_tflux,                                       &
#  endif
     &                  ad_t, ad_u, ad_v,                               &
# else
     &                  ad_ubar, ad_vbar,                               &
# endif
     &                  ad_zeta)

          IF (rec.gt.1) THEN
               HSS(ng)%ncid=ncidsav
          ELSE
               ncidsav=HSS(ng)%ncid
          END IF
          HSS(ng)%Rindex=rec-1
          CALL wrt_hessian (ng, L1, L1)

        END DO

      RETURN
      END SUBROUTINE rescale_HSS_tile


!
!***********************************************************************
      SUBROUTINE read_state (ng, tile, model,                           &
     &                       LBi, UBi, LBj, UBj, LBij, UBij,            &
     &                       Lwrk, rec,                                 &
     &                       ndef, ncfileid, ncname,                    &
# ifdef MASKING
     &                       rmask, umask, vmask,                       &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                       s_t_obc, s_u_obc, s_v_obc,                 &
#  endif
     &                       s_ubar_obc, s_vbar_obc,                    &
     &                       s_zeta_obc,                                &
# endif
# ifdef ADJUST_WSTRESS
     &                       s_ustr, s_vstr,                            &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                       s_tflux,                                   &
#  endif
     &                       s_t, s_u, s_v,                             &
# else
     &                       s_ubar, s_vbar,                            &
# endif
     &                       s_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
!
# ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti
# endif
# ifdef ADJUST_BOUNDARY
      USE nf_fread2d_bry_mod, ONLY : nf_fread2d_bry
#  ifdef SOLVE3D
      USE nf_fread3d_bry_mod, ONLY : nf_fread3d_bry
#  endif
# endif
      USE nf_fread2d_mod, ONLY : nf_fread2d
# ifdef SOLVE3D
      USE nf_fread3d_mod, ONLY : nf_fread3d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBij, UBij
      integer, intent(in) :: Lwrk, rec, ndef

      integer, intent(inout) :: ncfileid

      character (len=*), intent(in) :: ncname
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: s_t_obc(LBij:,:,:,:,:,:)
      real(r8), intent(inout) :: s_u_obc(LBij:,:,:,:,:)
      real(r8), intent(inout) :: s_v_obc(LBij:,:,:,:,:)
#   endif
      real(r8), intent(inout) :: s_ubar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: s_vbar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: s_zeta_obc(LBij:,:,:,:)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: s_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: s_vstr(LBi:,LBj:,:,:)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: s_tflux(LBi:,LBj:,:,:,:)
#   endif
      real(r8), intent(inout) :: s_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: s_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: s_v(LBi:,LBj:,:,:)
#  else
      real(r8), intent(inout) :: s_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: s_vbar(LBi:,LBj:,:)
#  endif
      real(r8), intent(inout) :: s_zeta(LBi:,LBj:,:)

# else

#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: s_t_obc(LBij:UBij,N(ng),4,             &
     &                                   Nbrec(ng),2,NT(ng))
      real(r8), intent(inout) :: s_u_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
      real(r8), intent(inout) :: s_v_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
#   endif
      real(r8), intent(inout) :: s_ubar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: s_vbar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: s_zeta_obc(LBij:UBij,4,Nbrec(ng),2)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: s_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: s_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: s_tflux(LBi:UBi,LBj:UBj,               &
     &                                   Nfrec(ng),2,NT(ng))
#   endif
      real(r8), intent(inout) :: s_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: s_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: s_v(LBi:UBi,LBj:UBj,N(ng),2)
#  else
      real(r8), intent(inout) :: s_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: s_vbar(LBi:UBi,LBj:UBj,3)
#  endif
      real(r8), intent(inout) :: s_zeta(LBi:UBi,LBj:UBj,3)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k
      integer :: ifield, it
      integer :: gtype, ncid, status, varid

      integer, dimension(4) :: Vsize

      real(r8) :: Fmin, Fmax, scale

# include "set_bounds.h"
!
      SourceFile='rescale_HSS.F, read_state'
!
!-----------------------------------------------------------------------
!  Read in requested model state record. Load data into state array
!  index Lwrk.
!-----------------------------------------------------------------------
!
!  Determine file and variables ids.
!
      IF (ndef.gt.0) THEN
        CALL netcdf_open (ng, model, ncname, 0, ncid)
        IF (exit_flag.ne.NoError) THEN
          WRITE (stdout,10) TRIM(ncname)
          RETURN
        END IF
        ncfileid=ncid
      ELSE
        ncid=ncfileid
      END IF

      DO i=1,4
        Vsize(i)=0
      END DO
!
!  Read in free-surface.
!
      gtype=r2dvar
      scale=1.0_r8
      CALL netcdf_inq_varid (ng, model, ncname, Vname(1,idFsur),        &
     &                       ncid, varid)
      IF (exit_flag.ne.NoError) RETURN

      status=nf_fread2d(ng, model, ncname, ncid,                        &
     &                  Vname(1,idFsur), varid,                         &
     &                  rec, gtype, Vsize,                              &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  scale, Fmin, Fmax,                              &
# ifdef MASKING
     &                  rmask,                                          &
# endif
     &                  s_zeta(:,:,Lwrk))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idFsur)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF

# ifdef ADJUST_BOUNDARY
!
!  Read in free-surface open boundaries.
!
      IF (ANY(Lobc(:,isFsur,ng))) THEN
        ifield=idSbry(isFsur)
        gtype=r2dvar
        scale=1.0_r8
        CALL netcdf_inq_varid (ng, model, ncname, Vname(1,ifield),      &
     &                         ncid, varid)
        IF (exit_flag.ne.NoError) RETURN

        status=nf_fread2d_bry (ng, model, ncname, ncid,                 &
     &                         Vname(1,ifield), varid,                  &
     &                         rec, gtype,                              &
     &                         LBij, UBij, Nbrec(ng),                   &
     &                         scale, Fmin, Fmax,                       &
     &                         s_zeta_obc(:,:,:,Lwrk))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,ifield)), rec, TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# ifndef SOLVE3D
!
!  Read in 2D U-momentum component.
!
      gtype=u2dvar
      scale=1.0_r8
      CALL netcdf_inq_varid (ng, model, ncname, Vname(1,idUbar),        &
     &                       ncid, varid)
      IF (exit_flag.ne.NoError) RETURN

      status=nf_fread2d(ng, model, ncname, ncid,                        &
     &                  Vname(1,idUbar), varid,                         &
     &                  rec, gtype, Vsize,                              &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  scale, Fmin, Fmax,                              &
#  ifdef MASKING
     &                  umask,                                          &
#  endif
     &                  s_ubar(:,:,Lwrk))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idUbar)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Read in 2D V-momentum component.
!
      gtype=v2dvar
      scale=1.0_r8
      CALL netcdf_inq_varid (ng, model, ncname, Vname(1,idVbar),        &
     &                       ncid, varid)
      IF (exit_flag.ne.NoError) RETURN

      status=nf_fread2d(ng, model, ncname, ncid,                        &
     &                  Vname(1,idVbar), varid,                         &
     &                  rec, gtype, Vsize,                              &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  scale, Fmin, Fmax,                              &
#  ifdef MASKING
     &                  vmask,                                          &
#  endif
     &                  s_vbar(:,:,Lwrk))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idVbar)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
# endif
# ifdef ADJUST_BOUNDARY
!
!  Read in 2D U-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isUbar,ng))) THEN
        ifield=idSbry(isUbar)
        gtype=u2dvar
        scale=1.0_r8
        CALL netcdf_inq_varid (ng, model, ncname, Vname(1,ifield),      &
     &                         ncid, varid)
        IF (exit_flag.ne.NoError) RETURN

        status=nf_fread2d_bry (ng, model, ncname, ncid,                 &
     &                         Vname(1,ifield), varid,                  &
     &                         rec, gtype,                              &
     &                         LBij, UBij, Nbrec(ng),                   &
     &                         scale, Fmin, Fmax,                       &
     &                         s_ubar_obc(:,:,:,Lwrk))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,ifield)), rec, TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Read in 2D V-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isVbar,ng))) THEN
        ifield=idSbry(isVbar)
        gtype=u2dvar
        scale=1.0_r8
        CALL netcdf_inq_varid (ng, model, ncname, Vname(1,ifield),      &
     &                         ncid, varid)
        IF (exit_flag.ne.NoError) RETURN

        status=nf_fread2d_bry (ng, model, ncname, ncid,                 &
     &                         Vname(1,ifield), varid,                  &
     &                         rec, gtype,                              &
     &                         LBij, UBij, Nbrec(ng),                   &
     &                         scale, Fmin, Fmax,                       &
     &                         s_vbar_obc(:,:,:,Lwrk))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,ifield)), rec, TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# ifdef ADJUST_WSTRESS
!
!  Read surface momentum stress.
!
      gtype=u3dvar
      scale=1.0_r8
      CALL netcdf_inq_varid (ng, model, ncname, Vname(1,idUsms),        &
     &                       ncid, varid)
      IF (exit_flag.ne.NoError) RETURN

      status=nf_fread3d(ng, model, ncname, ncid,                        &
     &                  Vname(1,idUsms), varid,                         &
     &                  rec, gtype, Vsize,                              &
     &                  LBi, UBi, LBj, UBj, 1, Nfrec(ng),               &
     &                  scale, Fmin, Fmax,                              &
#  ifdef MASKING
     &                  umask,                                          &
#  endif
     &                  s_ustr(:,:,:,Lwrk))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idUsms)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF

      gtype=v3dvar
      scale=1.0_r8
      CALL netcdf_inq_varid (ng, model, ncname, Vname(1,idVsms),        &
     &                       ncid, varid)
      IF (exit_flag.ne.NoError) RETURN

      status=nf_fread3d(ng, model, ncname, ncid,                        &
     &                  Vname(1,idVsms), varid,                         &
     &                  rec, gtype, Vsize,                              &
     &                  LBi, UBi, LBj, UBj, 1, Nfrec(ng),               &
     &                  scale, Fmin, Fmax,                              &
#  ifdef MASKING
     &                  vmask,                                          &
#  endif
     &                  s_vstr(:,:,:,Lwrk))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idVsms)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
# endif

# ifdef SOLVE3D
!
!  Read in 3D U-momentum component.
!
      gtype=u3dvar
      scale=1.0_r8
      CALL netcdf_inq_varid (ng, model, ncname, Vname(1,idUvel),        &
     &                       ncid, varid)
      IF (exit_flag.ne.NoError) RETURN

      status=nf_fread3d(ng, model, ncname, ncid,                        &
     &                  Vname(1,idUvel), varid,                         &
     &                  rec, gtype, Vsize,                              &
     &                  LBi, UBi, LBj, UBj, 1, N(ng),                   &
     &                  scale, Fmin, Fmax,                              &
#  ifdef MASKING
     &                  umask,                                          &
#  endif
     &                  s_u(:,:,:,Lwrk))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idUvel)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF

#  ifdef ADJUST_BOUNDARY
!
!  Read in 3D U-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isUvel,ng))) THEN
        ifield=idSbry(isUvel)
        gtype=u3dvar
        scale=1.0_r8
        CALL netcdf_inq_varid (ng, model, ncname, Vname(1,ifield),      &
     &                         ncid, varid)
        IF (exit_flag.ne.NoError) RETURN

        status=nf_fread3d_bry (ng, model, ncname, ncid,                 &
     &                         Vname(1,ifield), varid,                  &
     &                         rec, gtype,                              &
     &                         LBij, UBij, 1, N(ng), Nbrec(ng),         &
     &                         scale, Fmin, Fmax,                       &
     &                         s_u_obc(:,:,:,:,Lwrk))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,ifield)), rec, TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
!
!  Read in 3D U-momentum component.
!
      gtype=v3dvar
      scale=1.0_r8
      CALL netcdf_inq_varid (ng, model, ncname, Vname(1,idVvel),        &
     &                       ncid, varid)
      IF (exit_flag.ne.NoError) RETURN

      status=nf_fread3d(ng, model, ncname, ncid,                        &
     &                  Vname(1,idVvel), varid,                         &
     &                  rec, gtype, Vsize,                              &
     &                  LBi, UBi, LBj, UBj, 1, N(ng),                   &
     &                  scale, Fmin, Fmax,                              &
#  ifdef MASKING
     &                  vmask,                                          &
#  endif
     &                  s_v(:,:,:,Lwrk))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idVvel)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF

#  ifdef ADJUST_BOUNDARY
!
!  Read in 3D V-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isVvel,ng))) THEN
        ifield=idSbry(isVvel)
        gtype=u3dvar
        scale=1.0_r8
        CALL netcdf_inq_varid (ng, model, ncname, Vname(1,ifield),      &
     &                         ncid, varid)
        IF (exit_flag.ne.NoError) RETURN

        status=nf_fread3d_bry (ng, model, ncname, ncid,                 &
     &                         Vname(1,ifield), varid,                  &
     &                         rec, gtype,                              &
     &                         LBij, UBij, 1, N(ng), Nbrec(ng),         &
     &                         scale, Fmin, Fmax,                       &
     &                         s_v_obc(:,:,:,:,Lwrk))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,ifield)), rec, TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
!
!  Read in tracers.
!
      gtype=r3dvar
      scale=1.0_r8
      DO it=1,NT(ng)
        CALL netcdf_inq_varid (ng, model, ncname, Vname(1,idTvar(it)),  &
     &                         ncid, varid)
        IF (exit_flag.ne.NoError) RETURN

        status=nf_fread3d(ng, model, ncname, ncid,                      &
     &                    Vname(1,idTvar(it)), varid,                   &
     &                    rec, gtype, Vsize,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    scale, Fmin, Fmax,                            &
#  ifdef MASKING
     &                    rmask,                                        &
#  endif
     &                    s_t(:,:,:,Lwrk,it))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idTvar(it))), rec,           &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END DO

#  ifdef ADJUST_BOUNDARY
!
!  Read in tracers open boundaries.
!
      DO it=1,NT(ng)
        IF (ANY(Lobc(:,isTvar(it),ng))) THEN
          ifield=idSbry(isTvar(it))
          gtype=r3dvar
          scale=1.0_r8
          CALL netcdf_inq_varid (ng, model, ncname, Vname(1,ifield),    &
     &                           ncid, varid)
          IF (exit_flag.ne.NoError) RETURN

          status=nf_fread3d_bry (ng, model, ncname, ncid,               &
     &                           Vname(1,ifield), varid,                &
     &                           rec, gtype,                            &
     &                           LBij, UBij, 1, N(ng), Nbrec(ng),       &
     &                           scale, Fmin, Fmax,                     &
     &                           s_t_obc(:,:,:,:,Lwrk,it))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,ifield)), rec, TRIM(ncname)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
#  endif
#  ifdef ADJUST_STFLUX
!
!  Read in surface tracers flux.
!
      gtype=r3dvar
      scale=1.0_r8
      DO it=1,NT(ng)
        IF (Lstflux(it,ng)) THEN
          CALL netcdf_inq_varid (ng, model, ncname,                     &
     &                           Vname(1,idTsur(it)), ncid, varid)
          IF (exit_flag.ne.NoError) RETURN

          status=nf_fread3d(ng, model, ncname, ncid,                    &
     &                      Vname(1,idTsur(it)), varid,                 &
     &                      rec, gtype, Vsize,                          &
     &                      LBi, UBi, LBj, UBj, 1, Nfrec(ng),           &
     &                      scale, Fmin, Fmax,                          &
#   ifdef MASKING
     &                      rmask,                                      &
#   endif
     &                      s_tflux(:,:,:,Lwrk,it))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,idTsur(it))), rec,         &
     &                          TRIM(ncname)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
#  endif
# endif
!
!  If multiple files, close current file.
!
      IF (ndef.gt.0) THEN
        CALL netcdf_close (ng, model, ncid, ncname, .FALSE.)
      END IF
!
 10   FORMAT (' READ_STATE - unable to open NetCDF file: ',a)
 20   FORMAT (' READ_STATE - error while reading variable: ',a,2x,      &
     &        'at time record = ',i3,/,14x,'in NetCDF file: ',a)

      RETURN
      END SUBROUTINE read_state

#endif
      END MODULE rescale_HSS_mod
