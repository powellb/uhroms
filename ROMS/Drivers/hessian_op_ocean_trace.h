      MODULE ocean_control_mod
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2012 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  ROMS/TOMS Optimal Pertubations (Singular Vectors) Driver:           !
!                                                                      !
!  This driver computes the singular vectors of the propagator R(0,t)  !
!  which measure the  fastest  growing of all possible  perturbations  !
!  over a given time interval.  They  are  usually  used to study the  !
!  dynamics,  sensitivity,  and stability of the ocean circulation to  !
!  naturally occurring pertubations in the prediction system.          !
!                                                                      !
!  These  routines  control the  initialization,  time-stepping,  and  !
!  finalization of  ROMS/TOMS  model following ESMF conventions:       !
!                                                                      !
!     ROMS_initialize                                                  !
!     ROMS_run                                                         !
!     ROMS_finalize                                                    !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Moore, A.M. et al., 2004: A comprehensive ocean prediction and    !
!      analysis system based on the tangent linear and adjoint of a    !
!      regional ocean model, Ocean Modelling, 7, 227-258.              !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: ROMS_initialize
      PUBLIC  :: ROMS_run
      PUBLIC  :: ROMS_finalize

      CONTAINS

      SUBROUTINE ROMS_initialize (first, mpiCOMM)
!
!=======================================================================
!                                                                      !
!  This routine allocates and initializes ROMS/TOMS state variables    !
!  and internal and external parameters.                               !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
      USE mod_storage
      USE mod_fourdvar
      USE mod_netcdf

#ifdef MCT_LIB
!
# ifdef AIR_OCEAN
      USE ocean_coupler_mod, ONLY : initialize_ocn2atm_coupling
# endif
# ifdef WAVES_OCEAN
      USE ocean_coupler_mod, ONLY : initialize_ocn2wav_coupling
# endif
#endif
!
!  Imported variable declarations.
!
      logical, intent(inout) :: first

      integer, intent(in), optional :: mpiCOMM
!
!  Local variable declarations.
!
      logical :: allocate_vars = .TRUE.

#ifdef DISTRIBUTE
      integer :: MyError, MySize
#endif
      integer :: chunk_size, ng, thread
#ifdef _OPENMP
      integer :: my_threadnum
#endif

#ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Set distribute-memory (MPI) world communictor.
!-----------------------------------------------------------------------
!
      IF (PRESENT(mpiCOMM)) THEN
        OCN_COMM_WORLD=mpiCOMM
      ELSE
        OCN_COMM_WORLD=MPI_COMM_WORLD
      END IF
      CALL mpi_comm_rank (OCN_COMM_WORLD, MyRank, MyError)
      CALL mpi_comm_size (OCN_COMM_WORLD, MySize, MyError)
#endif
!
!-----------------------------------------------------------------------
!  On first pass, initialize model parameters a variables for all
!  nested/composed grids.  Notice that the logical switch "first"
!  is used to allow multiple calls to this routine during ensemble
!  configurations.
!-----------------------------------------------------------------------
!
      IF (first) THEN
        first=.FALSE.
!
!  Initialize parallel control switches. These scalars switches are
!  independent from standard input parameters.
!
        CALL initialize_parallel
!
!  Read in model tunable parameters from standard input. Allocate and
!  initialize variables in several modules after the number of nested
!  grids and dimension parameters are known.
!
        CALL inp_par (iTLM)
        IF (exit_flag.ne.NoError) RETURN
!
!  Set domain decomposition tile partition range.  This range is
!  computed only once since the "first_tile" and "last_tile" values
!  are private for each parallel thread/node.
!
!$OMP PARALLEL
#if defined _OPENMP
      MyThread=my_threadnum()
#elif defined DISTRIBUTE
      MyThread=MyRank
#else
      MyThread=0
#endif
      DO ng=1,Ngrids
        chunk_size=(NtileX(ng)*NtileE(ng)+numthreads-1)/numthreads
        first_tile(ng)=MyThread*chunk_size
        last_tile (ng)=first_tile(ng)+chunk_size-1
      END DO
!$OMP END PARALLEL
!
!  Initialize internal wall clocks. Notice that the timings does not
!  includes processing standard input because several parameters are
!  needed to allocate clock variables.
!
        IF (Master) THEN
          WRITE (stdout,10)
 10       FORMAT (/,' Process Information:',/)
        END IF
!
        DO ng=1,Ngrids
!$OMP PARALLEL
          DO thread=THREAD_RANGE
            CALL wclock_on (ng, iTLM, 0)
          END DO
!$OMP END PARALLEL
        END DO
!
!  Allocate and initialize modules variables.
!
!$OMP PARALLEL
        CALL mod_arrays (allocate_vars)
!$OMP END PARALLEL
!
!  Allocate and initialize 4D-Var arrays.
!
        CALL initialize_fourdvar

      END IF

#if defined MCT_LIB && (defined AIR_OCEAN || defined WAVES_OCEAN)
!
!-----------------------------------------------------------------------
!  Initialize coupling streams between model(s).
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
# ifdef AIR_OCEAN
        CALL initialize_ocn2atm_coupling (ng, MyRank)
# endif
# ifdef WAVES_OCEAN
        CALL initialize_ocn2wav_coupling (ng, MyRank)
# endif
      END DO
#endif
!
!-----------------------------------------------------------------------
!  Initialize tangent linear for all grids first in order to compute
!  the size of the state vector, Nstate.  This size is computed in
!  routine "wpoints".
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
#if defined BULK_FLUXES && defined NL_BULK_FLUXES
        BLK(ng)%name=FWD(ng)%name
#endif
!$OMP PARALLEL
        CALL tl_initial (ng)
!$OMP END PARALLEL
        IF (exit_flag.ne.NoError) RETURN
      END DO
!
!-----------------------------------------------------------------------
!  Read in Lanczos algorithm coefficients (cg_beta, cg_delta) from
!  file LCZ(ng)%name NetCDF (IS4DVAR adjoint file), as computed in the
!  I4D-Var Lanczos data assimilation algorithm for the first outer
!  loop.  They are needed here, in routine "tl_inner2state", to compute
!  the tangent linear model initial conditions as the weighted sum
!  of the Lanczos vectors. The weighting coefficient are computed
!  by solving a tri-diagonal system that uses "cg_beta" and "cg_gamma".
!-----------------------------------------------------------------------
!
      SourceFile='obs_sen_ocean.h, ROMS_initialize'

      DO ng=1,Ngrids
        CALL netcdf_get_fvar (ng, iADM, LCZ(ng)%name, 'cg_beta',        &
     &                        cg_beta)
        IF (exit_flag.ne. NoError) RETURN
        CALL netcdf_get_fvar (ng, iADM, LCZ(ng)%name, 'cg_delta',       &
     &                        cg_delta)
        IF (exit_flag.ne. NoError) RETURN
      END DO
!
!  Allocate arrays associated with Generalized Stability Theory (GST)
!  analysis.
!
      CALL allocate_storage
!
!  Initialize various IO flags.
!
      Nrun=0
      DO ng=1,Ngrids
        LdefADJ(ng)=.TRUE.
        LwrtADJ(ng)=.TRUE.
        LdefTLM(ng)=.TRUE.
        LwrtTLM(ng)=.TRUE.
        LwrtPER(ng)=.FALSE.
        LcycleTLM(ng)=.FALSE.
        LcycleADJ(ng)=.FALSE.
        nADJ(ng)=ntimes(ng)
        nTLM(ng)=ntimes(ng)
      END DO
!
!  Initialize ARPACK parameters.
!
      Lrvec=.TRUE.                ! Compute Ritz vectors
      bmat='I'                    ! standard eigenvalue problem
      which='LM'                  ! compute NEV largest eigenvalues
      howmany='A'                 ! compute NEV Ritz vectors
      DO ng=1,Ngrids
        ido(ng)=0                 ! reverse communication flag
        info(ng)=0                ! random initial residual vector
        iparam(1,ng)=1            ! exact shifts
        iparam(3,ng)=MaxIterGST   ! maximum number of Arnoldi iterations
        iparam(4,ng)=1            ! block size in the recurrence
        iparam(7,ng)=1            ! type of eigenproblem being solved
      END DO
!
!  ARPACK debugging parameters.
!
      logfil=stdout               ! output logical unit
      ndigit=-3                   ! number of decimal digits
      msaupd=1                    ! iterations, timings, Ritz
      msaup2=1                    ! norms, Ritz values
      msaitr=0
      mseigt=0
      msapps=0
      msgets=0
      mseupd=0
!
!  Determine size of the eigenproblem (Nsize) and size of work space
!  array SworkL (LworkL).
!
      DO ng=1,Ngrids
        Nconv(ng)=0
        Nsize(ng)=Ninner
      END DO

#ifdef CHECKPOINTING
!
!  If restart, read in checkpointing data GST restart NetCDF file.
!  Otherwise, create checkpointing restart NetCDF file.
!
      DO ng=1,Ngrids
        IF (LrstGST) THEN
          CALL get_gst (ng, iTLM)
          ido(ng)=-2
          laup2(1)=.FALSE.        ! cnorm
          laup2(2)=.FALSE.        ! getv0
          laup2(3)=.FALSE.        ! initv
          laup2(4)=.FALSE.        ! update
          laup2(5)=.TRUE.         ! ushift
        ELSE
          CALL def_gst (ng, iTLM)
        END IF
        IF (exit_flag.ne.NoError) RETURN
      END DO
#endif

      RETURN
      END SUBROUTINE ROMS_initialize

      SUBROUTINE ROMS_run (RunInterval)
!
!=======================================================================
!                                                                      !
!  This routine computes the singular vectors of R(0,t) by a single    !
!  integration  of  a  perturbation "u"  forward  in time  with the    !
!  tangent linear model over [0,t],  multiplication  of the  result    !
!  by "X",  followed by an integration of the  result  backwards in    !
!  time with the  adjoint model over [t,0].  This is  equivalmet to    !
!  the matrix-vector operation:                                        !
!                                                                      !
!       transpose[R(t,0)] X R(0,t) u                                   !
!                                                                      !
!  The above operator is symmetric and the  ARPACK library is used     !
!  to select eigenvectors and eigenvalues:                             !
!                                                                      !
!  Lehoucq, R.B., D.C. Sorensen, and C. Yang, 1997:  ARPACK user's     !
!    guide:  solution  of  large  scale  eigenvalue  problems with     !
!    implicit restarted Arnoldi Methods, Rice University, 140p.        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
      USE mod_stepping
      USE mod_storage
!
      USE propagator_mod
      USE white_noise_mod
#ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcastf, mp_bcasti
#endif
      USE packing_mod, ONLY : r_norm2
!
!  Imported variable declarations
!
      real(r8), intent(in) :: RunInterval            ! seconds
!
!  Local variable declarations.
!

      integer :: i, iter, ng, Rscheme

      real(r8) :: Rmin, Rmax, cff

      real(r8), dimension(Ninner) :: state, ad_state

      character (len=55) :: string

#ifdef LCZ_FINAL
      DO ng=1,Ngrids
        LdefHSS(ng)=.TRUE.
        CALL def_hessian (ng)
        IF (exit_flag.ne.NoError) RETURN
      END DO
#endif
!
      DO ng=1,Ngrids
        cff=0.0_r8
        DO iter=1,Ninner
          IF (Master) THEN
!AMM        Rscheme=1
!            Rscheme=0
!            CALL white_noise1d(ng,iTLM,Rscheme,1,Ninner,1,Ninner,       &
!     &                         Rmin,Rmax,state)
!            DO i=1,Ninner
!              IF (state(i).lt.0.0_r8) THEN
!                state(i)=-1.0_r8
!              ELSE
!                state(i)=1.0_r8
!              END IF
!            END DO
             DO i=1,Ninner
              state(i)=0.0_r8
             END DO
             state(iter)=1.0_r8
             print *,'TRACE: iter=',iter,' state=',state
          END IF
#ifdef DISTRIBUTE
          CALL mp_bcastf (ng, iTLM, state)
#endif
          CALL propagator (RunInterval, state, ad_state)
          IF (Master) THEN
!            cff=0.0_r8
!            DO i=1,Ninner
!              cff=cff+state(i)*ad_state(i)
!            END DO
            cff=cff+ad_state(iter)
            print *,'Trace estimate: iter=',iter,cff
          END IF
        END DO
      END DO

      RETURN
      END SUBROUTINE ROMS_run

      SUBROUTINE ROMS_finalize
!
!=======================================================================
!                                                                      !
!  This routine terminates ROMS/TOMS nonlinear and adjoint models      !
!  execution.                                                          !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
!
!  Local variable declarations.
!
      integer :: Fcount, ng, thread
!
!-----------------------------------------------------------------------
!  If blowing-up, save latest model state into RESTART NetCDF file.
!-----------------------------------------------------------------------
!
!  If cycling restart records, write solution into the next record.
!
      IF (exit_flag.eq.1) THEN
        DO ng=1,Ngrids
          IF (LwrtRST(ng)) THEN
            IF (Master) WRITE (stdout,10)
 10         FORMAT (/,' Blowing-up: Saving latest model state into ',   &
     &                ' RESTART file',/)
            Fcount=RST(ng)%Fcount
            IF (LcycleRST(ng).and.(RST(ng)%Nrec(Fcount).ge.2)) THEN
              RST(ng)%Rindex=2
              LcycleRST(ng)=.FALSE.
            END IF
            blowup=exit_flag
            exit_flag=NoError
            CALL wrt_rst (ng)
          END IF
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Stop model and time profiling clocks.  Close output NetCDF files.
!-----------------------------------------------------------------------
!
!  Stop time clocks.
!
      IF (Master) THEN
        WRITE (stdout,20)
 20     FORMAT (/,' Elapsed CPU time (seconds):',/)
      END IF

      DO ng=1,Ngrids
!$OMP PARALLEL
        DO thread=THREAD_RANGE
          CALL wclock_off (ng, iTLM, 0)
        END DO
!$OMP END PARALLEL
      END DO
!
!  Close IO files.
!
      CALL close_out

      RETURN
      END SUBROUTINE ROMS_finalize

      SUBROUTINE IRAM_error (info, string)
!
!=======================================================================
!                                                                      !
!  This routine decodes internal error messages from the Implicit      !
!  Restarted Arnoldi Method (IRAM) for the computation of optimal      !
!  perturbation Ritz eigenfunctions.                                   !
!                                                                      !
!=======================================================================
!
!
!  imported variable declarations.
!
      integer, intent(in) :: info

      character (len=*), intent(out) :: string
!
!-----------------------------------------------------------------------
!  Decode error message from IRAM.
!-----------------------------------------------------------------------
!
      IF (info.eq.0)  THEN
        string='Normal exit                                            '
      ELSE IF (info.eq.1) THEN
        string='Maximum number of iterations taken                     '
      ELSE IF (info.eq.3) THEN
        string='No shifts could be applied during an IRAM cycle        '
      ELSE IF (info.eq.-1) THEN
        string='Nstate must be positive                                '
      ELSE IF (info.eq.-2) THEN
        string='NEV must be positive                                   '
      ELSE IF (info.eq.-3) THEN
        string='NCV must be greater NEV and less than or equal Nstate  '
      ELSE IF (info.eq.-4) THEN
        string='Maximum number of iterations must be greater than zero '
      ELSE IF (info.eq.-5) THEN
        string='WHICH must be one of LM, SM, LA, SA or BE              '
      ELSE IF (info.eq.-6) THEN
        string='BMAT must be one of I or G                             '
      ELSE IF (info.eq.-7) THEN
        string='Length of private work array SworkL is not sufficient  '
      ELSE IF (info.eq.-8) THEN
        string='Error in DSTEQR in the eigenvalue calculation          '
      ELSE IF (info.eq.-9) THEN
        string='Starting vector is zero                                '
      ELSE IF (info.eq.-10) THEN
        string='IPARAM(7) must be 1, 2, 3, 4, 5                        '
      ELSE IF (info.eq.-11) THEN
        string='IPARAM(7) = 1 and BMAT = G are incompatable            '
      ELSE IF (info.eq.-12) THEN
        string='IPARAM(1) must be equal to 0 or 1                      '
      ELSE IF (info.eq.-13) THEN
        string='NEV and WHICH = BE are incompatable                    '
      ELSE IF (info.eq.-14) THEN
        string='Did not find any eigenvalues to sufficient accuaracy   '
      ELSE IF (info.eq.-15) THEN
        string='HOWMANY must be one of A or S if RVEC = .TRUE.         '
      ELSE IF (info.eq.-16) THEN
        string='HOWMANY = S not yet implemented                        '
      ELSE IF (info.eq.-17) THEN
        string='Different count of converge Ritz values in DSEUPD      '
      ELSE IF (info.eq.-9999) THEN
        string='Could not build and Arnoldi factorization              '
      END IF

      RETURN
      END SUBROUTINE IRAM_error

      END MODULE ocean_control_mod
