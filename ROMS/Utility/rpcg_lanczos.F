#include "cppdefs.h"

#if defined WEAK_CONSTRAINT && defined RPCG
      SUBROUTINE rpcg_lanczos (ng, model, outLoop, innLoop, NinnLoop,   &
     &                    Lcgini)
!
!svn $Id: rpcg_lanczos.F 527 2011-01-06 19:42:35Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2011 The ROMS/TOMS Group       Selime Gurol      !
!    Licensed under a MIT/X style license            Andrew M. Moore   !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Weak Constraint 4-Dimensional Variational (4DVar) Restricted        !
!      Pre-conditioned Conjugate Gradient (RPCG) Algorithm             !

# if defined TL_W4DVAR || defined W4DVAR || defined W4DVAR_SENSITIVITY
!                                                                      !
!  The indirect representer method solves the system:                  !
!                                                                      !
!              (R_n + Cobs) * Beta_n = h_n                             !
!                                                                      !
!              h_n = Xo - H * X_n                                      !
!                                                                      !
!  where R_n is the representer matrix, Cobs is the observation-error  !
!  covariance,  Beta_n  are the representer coefficients,  h_n is the  !
!  misfit between observations (Xo) and model (H*X_n),  and  H is the  !
!  linearized observation operator. Here, _n denotes iteration.        !
!                                                                      !
!  This system does not need to be solved explicitly by inverting the  !
!  symmetric stabilized representer matrix, P_n:                       !
!                                                                      !
!              P_n = R_n + Cobs                                        !
!                                                                      !
!  but by computing the action of P_n on any vector PSI, such that     !
!                                                                      !
!              P_n * PSI = R_n * PSI + Cobs * PSI                      !
!                                                                      !
!  The representer matrix is not explicitly computed but evaluated by  !
!  one integration backward of the adjoint model  and one integration  !
!  forward of the tangent linear model for any forcing vector PSI.     !
!                                                                      !
!  Notice that "ObsScale" vector is used for screenning observations.  !
!  This scale is one (zero) for good (bad) observations.               !
!                                                                      !
!  Currently, parallelization of this algorithm is not needed because  !
!  each parallel node has a full copy of the assimilation vectors.     !
!                                                                      !
!  This code solves Ax=b by minimizing the cost function               !
!  0.5*x*A*x-x*b, assuming an initial guess of x=x0. In this case the  !
!  gradient is Ax-b and the stabilized representer matrix is A.        !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Chua, B. S. and A. F. Bennett,  2001:  An inverse ocean modeling  !
!      sytem, Ocean Modelling, 3, 137-165.                             !

# elif defined TL_W4DPSAS || \
       defined W4DPSAS    || defined W4DPSAS_SENSITIVITY
!                                                                      !
!  Solve the system (following Courtier, 1997):                        !
!                                                                      !
!              (H M_n B (M_n)' H' + Cobs) * w_n = d_n                  !
!                                                                      !
!              d_n = yo - H * Xb_n                                     !
!                                                                      !
!  where  M_n is the tangent linear model matrix and  _n denotes a     !
!  sequence of outer-loop estimates, Cobs is the observation-error     !
!  covariance,  B is the background error  covariance,  d_n is the     !
!  misfit between observations (yo) and model (H * Xb_n), and H is     !
!  the linearized observation operator.                                !
!                                                                      !
!  The analysis increment is:                                          !
!                                                                      !
!             dx_n = B M' H' w_n                                       !
!                                                                      !
!  so that Xa = Xb + dx_n.                                             !
!                                                                      !
!  The system does not need to be  solved explicitly  by inverting     !
!  the symmetric matrix, P_n:                                          !
!                                                                      !
!              P_n = H M_n B (M_n)' H' + Cobs                          !
!                                                                      !
!  but by computing the action of P_n on any vector PSI, such that     !
!                                                                      !
!              P_n * PSI =  H M_n B (M_n)' H' * PSI + Cobs * PSI       !
!                                                                      !
!  The (H M_n B (M_n)' H') matrix is not  explicitly computed  but     !
!  evaluated by  one integration backward of the adjoint model and     !
!  one  integration  forward of the  tangent linear model  for any     !
!  forcing vector PSI.                                                 !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Courtier, P., 1997: Dual formulation of four-dimensional          !
!      variational assimilation, Quart. J. Roy. Meteor. Soc.,          !
!      123, 2449-2461.                                                 !
# endif
!                                                                      !
!  The restricted preconditioned conjugate gradient (RPCG) algorithm   !
!  is used to compute the solution. Specfically, the following system  !
!  is solved:                                                          !
!                                                                      !
!              (Cobs^-1 H M_n B (M_n)' H'+I)*w=Cobs^-1 d_n             !
!                                                                      !
!  All inner-products are defined with respect to a norm involving     !
!  H M_n B (M_n)' H', and the convergence of J in model space is       !
!  guaranteed to converge at the same rate as I4DVAR in primal space.  !
!                                                                      !
!  RPCG Reference:                                                     !
!                                                                      !
!    Gratton, S. and J. Tshimanga, 2009: An observation-space          !
!      formulation of variational assimilation using a restricted      !
!      preconditioned conjugate gradient algorithm. QJRMS, 135,        !
!      1573-1585.
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_iounits
      USE mod_scalars

# ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcastf, mp_bcasti, mp_bcastl
# endif
      implicit none
!
!  Imported variable declarations
!
      logical, intent(in) :: Lcgini
      integer, intent(in) :: ng, model, outLoop, innLoop, NinnLoop
!
!  Local variable declarations.
!
      logical :: Ltrans, Laug
      integer :: i, ic, j, iobs, ivec, Lscale, info

      real(r8) :: zbet, eps, preducv, preducy
      real(r8) :: Jopt, Jf, Jmod, Jdata, Jb, Jobs, Jact, cff

      real(r8), dimension(NinnLoop) :: zu, zgam, zfact, ztemp3
      real(r8), dimension(NinnLoop) :: ztemp1, ztemp2, zu1, zu2
      real(r8), dimension(NinnLoop) :: ztemp4
      real(r8), dimension(Ndatum(ng)+1) :: px, pgrad, zw, zt
      real(r8), dimension(Ndatum(ng)+1) :: zhv, zht, zd
      real(r8), dimension(Ninner,3) :: zwork
      real(r8), dimension(2*(NinnLoop-1)) :: work
      real(r8), dimension(Ninner,Ninner) :: zgv

      character (len=13) :: string
!
!=======================================================================
!  Weak constraint 4DVar conjugate gradient, Lanczos-based algorithm.
!=======================================================================
!
!  This conjugate gradient algorithm is not run in parallel since the
!  weak constraint is done in observation space. Mostly all the
!  variables are 1D arrays. Therefore, in parallel applications (only
!  distributed-memory is possible) the master node does the computations
!  and then broadcasts results to remaining nodes in the communicator.
!
!  This version of congrad solves A(x+x0)=b for x, by minimizing
!  J=0.5*x'Ax-x'(b+Ax0), where x0 is a first-guess estimate of the
!  representer coefficients from the previous outer-loop.
!  For the first outer-loop, x0=0. In the code x0=cg_pxsave and
!  x=px.
!
      Ltrans=.FALSE.
!
      MASTER_THREAD : IF (Master) THEN
!
!  Initialize internal parameters.  This is needed here for output
!  reasons.
!
        Jopt=0.0_r8
        Jf=0.0_r8
        Jdata=0.0_r8
        Jmod=0.0_r8
        Jb=Jb0(outLoop-1)
        Jobs=0.0_r8
        Jact=0.0_r8
        eps=0.0_r8
        preducv=0.0_r8
        preducy=0.0_r8
!
!  Multiply TLmodVal by ObsScale to effectively remove any rejected
!  observations from the analysis.
!
        IF (innLoop.gt.0) THEN
          DO iobs=1,Ndatum(ng)
            TLmodVal(iobs)=ObsScale(iobs)*TLmodVal(iobs)
          END DO
        END IF
!
!-----------------------------------------------------------------------
!  Initialization step, innloop = 0.
!-----------------------------------------------------------------------
!
        MINIMIZE : IF (innLoop.eq.0) THEN

# if defined W4DPSAS             || defined TL_W4DPSAS || \
     defined W4DPSAS_SENSITIVITY
!
          DO iobs=1,Ndatum(ng)
            BCKmodVal(iobs)=NLmodVal(iobs)
          END DO
# endif
!
!  For the first outer-loop, x0=0.
!
          IF (outLoop.eq.1) THEN
            DO iobs=1,Ndatum(ng)+1
              px(iobs)=0.0_r8
            END DO
            DO iobs=1,Ndatum(ng)
              Hbk(iobs,outLoop)=0.0_r8
            END DO
          ELSE
!! NOTE: CODE WON'T WORK FOR W4DVAR AT THE MOMENT SINCE THIS IS
!! THE WRONG TLmodVal HERE!!!
            DO iobs=1,Ndatum(ng)
              Hbk(iobs,outLoop)=-TLmodVal(iobs)
            END DO
            IF (Lprecond) THEN
              DO iobs=1,Ndatum(ng)+1
                FOURDVAR(ng)%cg_Dpxsum(iobs)=FOURDVAR(ng)%cg_pxsum(iobs)
              END DO
              Lscale=1                 ! Spectral LMP
              Laug=.FALSE.
              CALL rprecond (ng, Lscale, Laug, outLoop, NinnLoop-1,     &
     &                       FOURDVAR(ng)%cg_Dpxsum)
              DO iobs=1,Ndatum(ng)+1
                FOURDVAR(ng)%cg_Dpxsum(iobs)=                           &
     &         -FOURDVAR(ng)%cg_Dpxsum(iobs)+FOURDVAR(ng)%cg_pxsum(iobs)
              END DO
            END IF
          END IF
!
!  If this is the first pass of the inner loop, set up the vectors and
!  store the first guess. The initial starting guess is assumed to be
!  zero in which case the residual is just: (obs-model).
!
            DO iobs=1,Ndatum(ng)
!
!  Selime's code:  r=b.
!
# if defined W4DPSAS             || defined TL_W4DPSAS || \
     defined W4DPSAS_SENSITIVITY
              pgrad(iobs)=ObsScale(iobs)*                               &
     &                    (ObsVal(iobs)-BCKmodVal(iobs))
# else
              pgrad(iobs)=ObsScale(iobs)*                               &
     &                    (ObsVal(iobs)-TLmodVal(iobs))
# endif
              IF (ObsErr(iobs).NE.0.0_r8) THEN
                 pgrad(iobs)=pgrad(iobs)/ObsErr(iobs)
              END IF
              vgrad0(iobs)=pgrad(iobs)
              vgrad0s(iobs)=vgrad0(iobs)
            END DO
!
!AUGMENTED
!  Augment the residual vector.
!
            pgrad(Ndatum(ng)+1)=1.0_r8
            vgrad0(Ndatum(ng)+1)=pgrad(Ndatum(ng)+1)
            vgrad0s(Ndatum(ng)+1)=vgrad0(Ndatum(ng)+1)
!
!  If preconditioning, convert pgrad from v-space to y-space.
!
!Selime's code: z1=G*r.
!
            IF (Lprecond.and.(outLoop.gt.1)) THEN
              Lscale=1                 ! Spectral LMP
              Laug=.TRUE.
              CALL rprecond (ng, Lscale, Laug, outLoop, NinnLoop-1,     &
     &                       pgrad)
            END IF

!
!   Selime's code: zgrad0=z1.
!
!   Selime's code: vgrad0=r.
!
            DO iobs=1,Ndatum(ng)+1
              zgrad0(iobs,outLoop)=pgrad(iobs)
            END DO

!
!   Selime's code: ADmodVal=z1.
!
            DO iobs=1,Ndatum(ng)
              ADmodVal(iobs)=zgrad0(iobs,outLoop)
            END DO
!
          preducv=1.0_r8
          preducy=1.0_r8
!
!-----------------------------------------------------------------------
!  Minimization step, innLoop > 0.
!-----------------------------------------------------------------------
!
        ELSE
!
!  COMPLETE THE CALCULATION OF PREVIOUS LANCZOS VECTOR.
!
!  Selime's code: t1=GDG'*z1=TLmodVal.
!
!
!  Complete the computation of the Lanczos vector from previous inner-loop
!
          IF (innLoop.EQ.1) THEN
            cg_Gnorm_v(outLoop)=0.0_r8
            DO iobs=1,Ndatum(ng)
              cg_Gnorm_v(outLoop)=cg_Gnorm_v(outLoop)+                  &
     &                            TLmodVal(iobs)*vgrad0(iobs)
            END DO
!
! AUGMENTED. Here Hbk=H(x0-xb). Hbk=0 and Jb0=0 when outer=1.
! Jb0 is computed using sum_i (G'*w)_i, i.e. the sum of the outputs 
! of the adjoint model at the end of AD_LOOP2.
!
            DO iobs=1,Ndatum(ng)
              cg_Gnorm_v(outLoop)=cg_Gnorm_v(outLoop)+                  &
     &            Hbk(iobs,outLoop)*zgrad0(Ndatum(ng)+1,outLoop)*       &
     &            vgrad0(iobs)+  &
     &            Hbk(iobs,outLoop)*vgrad0(Ndatum(ng)+1)*               &
     &            zgrad0(iobs,outLoop)
            END DO
            cg_Gnorm_v(outLoop)=cg_Gnorm_v(outLoop)+                    &
     &        Jb0(outLoop-1)*vgrad0(Ndatum(ng)+1)*                      &
     &        zgrad0(Ndatum(ng)+1,outLoop)
!
            cg_Gnorm_v(outLoop)=SQRT(cg_Gnorm_v(outLoop))
!
!   Selime's code: beta_first=sqrt(r't1)=cg_Gnorm_v.
!                  cg_QG(1)=beta_first
!
            cg_QG(1,outLoop)=cg_Gnorm_v(outLoop)
!
! Normalize TLmodVal here and save in TLmodVal_S.
!
!AUGMENTED.
!
            DO iobs=1,Ndatum(ng)+1
!
!   Selime's code: v1=r/beta_first=zcglwk(1).
!
              zcglwk(iobs,1,outLoop)=vgrad0(iobs)/                      &
     &                                          cg_Gnorm_v(outLoop)
!
!   Selime's code: z1=z1/beta_first=vcglwk(1).
!
              vcglwk(iobs,1,outLoop)=zgrad0(iobs,outLoop)/              &
     &                                          cg_Gnorm_v(outLoop)
            END DO
!
!AUGMENTED
!
            DO iobs=1,Ndatum(ng)
!
!   Selime's code: t1=t1/beta_first=TLmodVal
!
              TLmodVal_S(iobs,innLoop,outLoop)=TLmodVal(iobs)
              TLmodVal(iobs)=TLmodVal(iobs)/cg_Gnorm_v(outLoop)
            END DO
!
            preducv=1.0_r8
            preducy=1.0_r8
!
! Compute the actual penalty function terms.
!
            Jobs=0.0_r8
!
! NO AUGMENTED TERM HERE.
!
            Jact=Jb
            DO iobs=1,Ndatum(ng)
              IF (ObsErr(iobs).NE.0.0_r8) THEN
                Jact=Jact+ObsErr(iobs)*vgrad0s(iobs)*vgrad0s(iobs)
              END IF
            END DO
            Jobs=Jact-Jb
!
          ELSE
!
!  First grab the non-normalized Lanczos vector.
!
!
!AUGMENTED.
!
            DO iobs=1,Ndatum(ng)+1
!
!  Selime's code: w=pgrad, t1=TLmodVal.
!
              pgrad(iobs)=zcglwk(iobs,innLoop,outLoop)
            END DO
!
!  Selime's code:  beta=pgrad*HBH'*pgrad=cg_beta
!
            cg_beta(innLoop,outLoop)=0.0_r8
            DO iobs=1,Ndatum(ng)
              cg_beta(innLoop,outLoop)=cg_beta(innLoop,outLoop)+        &
     &                                 pgrad(iobs)*TLmodVal(iobs)
            END DO
!
!AUGMENTED.
!
            DO iobs=1,Ndatum(ng)
              cg_beta(innLoop,outLoop)=cg_beta(innLoop,outLoop)+        &
     &             Hbk(iobs,outLoop)*                                   &
     &             vcglwk(Ndatum(ng)+1,innLoop,outLoop)*                &
     &             pgrad(iobs)+                                         &
     &             Hbk(iobs,outLoop)*pgrad(Ndatum(ng)+1)*               &
     &             vcglwk(iobs,innLoop,outLoop)
            END DO
            cg_beta(innLoop,outLoop)=cg_beta(innLoop,outLoop)+          &
     &                        Jb0(outLoop-1)*pgrad(Ndatum(ng)+1)*       &
     &                        vcglwk(Ndatum(ng)+1,innLoop,outLoop)
!
            cg_beta(innLoop,outLoop)=SQRT(cg_beta(innLoop,outLoop))
!
! Normalize TLmodVal here and save in TLmodVal_S.
!
!AUGMENTED
            DO iobs=1,Ndatum(ng)+1
!
!   Selime's code: v1=w/beta=zcglwk
!
              zcglwk(iobs,innLoop,outLoop)=pgrad(iobs)/                 &
     &                                       cg_beta(innLoop,outLoop)
!   Selime's code: z1=w/beta=vcglwk
              vcglwk(iobs,innLoop,outLoop)=vcglwk(iobs,innLoop,outLoop)/&
     &                                       cg_beta(innLoop,outLoop)
            END DO
!
!   Selime's code: t1=t1/beta=TLmodVal
!
!AUGMENTED
!
            DO iobs=1,Ndatum(ng)
              TLmodVal_S(iobs,innLoop,outLoop)=TLmodVal(iobs)
              TLmodVal(iobs)=TLmodVal(iobs)/cg_beta(innLoop,outLoop)
            END DO
!
!   Selime's code:  cg_QG=zgrad0*HBH'*zcglwk
!
            cg_QG(innLoop,outLoop)=0.0_r8
            DO iobs=1,Ndatum(ng)
              cg_QG(innLoop,outLoop)=cg_QG(innLoop,outLoop)+            &
     &                                 TLmodVal(iobs)*                  &
     &                                 vgrad0(iobs)
            END DO
!
!AUGMENTED.
!
            DO iobs=1,Ndatum(ng)
              cg_QG(innLoop,outLoop)=cg_QG(innLoop,outLoop)+            &
     &                Hbk(iobs,outLoop)*                                &
     &                vcglwk(Ndatum(ng)+1,innLoop,outLoop)*             &
     &                vgrad0(iobs)                         +            &
     &                Hbk(iobs,outLoop)*vcglwk(iobs,innLoop,outLoop)*   &
     &                vgrad0(Ndatum(ng)+1)
            END DO
            cg_QG(innLoop,outLoop)=cg_QG(innLoop,outLoop)+              &
     &           Jb0(outLoop-1)*vcglwk(Ndatum(ng)+1,innLoop,outLoop)*   &
     &           vgrad0(Ndatum(ng)+1)
!TEST
!TEST            IF (innLoop.gt.1) cg_QG(innLoop,outLoop)=0.0_r8
!TEST
!
!  Calculate the new solution based upon the new, orthonormalized
!  Lanczos vector. First, the tridiagonal system is solved by
!  decomposition and forward/back substitution.
!
            zbet=cg_delta(1,outLoop)
            zu(1)=cg_QG(1,outLoop)/zbet
            DO ivec=2,innLoop-1
              zgam(ivec)=cg_beta(ivec,outLoop)/zbet
              zbet=cg_delta(ivec,outLoop)-cg_beta(ivec,outLoop)*        &
     &             zgam(ivec)
              zu(ivec)=(cg_QG(ivec,outLoop)-cg_beta(ivec,outLoop)*      &
     &                  zu(ivec-1))/zbet
            END DO
            zwork(innLoop-1,3)=zu(innLoop-1)

            DO ivec=innLoop-2,1,-1
              zu(ivec)=zu(ivec)-zgam(ivec+1)*zu(ivec+1)
              zwork(ivec,3)=zu(ivec)
            END DO
!
!AUGMENTED.
!
            DO iobs=1,Ndatum(ng)+1
              zw(iobs)=zgrad0(iobs,outLoop)+                            &
     &                 cg_beta(innLoop,outLoop)*                        &
     &                 vcglwk(iobs,innLoop,outLoop)*                    &
     &                 zwork(innLoop-1,3)
            END DO

!
!AUGMENTED.
!
            DO iobs=1,Ndatum(ng)+1
              px(iobs)=0.0_r8
              DO ivec=1,innLoop-1
                px(iobs)=px(iobs)+zcglwk(iobs,ivec,outLoop)*            &
     &                            zwork(ivec,3)
                zw(iobs)=zw(iobs)+zcglwk(iobs,ivec,outLoop)*            &
     &                                         cg_QG(ivec,outLoop)
              END DO
            END DO
!
!  If preconditioning, convert px from y-space to v-space.
!  We will always keep px in v-space.
!
            IF (Lprecond.and.(outLoop.gt.1)) THEN
              Lscale=1                 ! Spectral LMP
              Laug=.TRUE.
              CALL rprecond (ng, Lscale, Laug, outLoop, NinnLoop-1, px)
            END IF
!
!  Compute the reduction in the gradient Ax-b (in y-space).
!
!AMM            preducy=0.0_r8
!AMM            DO iobs=1,Ndatum(ng)
!AMM              preducy=preducy+zw(iobs)*zw(iobs)
!AMM            END DO
!AMM            preducy=SQRT(preducy)/cg_Gnorm_y(outLoop)
!
!  Compute the reduction in the gradient Ax-b (in v-space).
!
            IF (Lprecond.and.(outLoop.gt.1)) THEN
              Lscale=1                 ! Spectral LMP
              Laug=.TRUE.
              CALL rprecond (ng, Lscale, Laug, outLoop, NinnLoop-1, zw)
            END IF
!
            preducv=0.0_r8
!
!AUGMENTED.
!
            DO iobs=1,Ndatum(ng)+1
              preducv=preducv+zw(iobs)*zw(iobs)
            END DO
            preducv=SQRT(preducv)/cg_Gnorm_v(outLoop)
!
!  Estimate the residual (in y-space) by:
!
!    ||Ax-b|| = cg_beta(k+1,outLoop)*zwork(innLoop,3)
!
          eps=ABS(cg_beta(innLoop,outLoop)*zwork(innLoop-1,3))
!
!  Estimate the various penalty function components based on
!  Bennett (2002), page 22 and pages 42-44:
!
!     Jopt  = value of the penalty function when true px
!             has been identified.
!     Jf    = initial data misfit for the first-guess
!     Jdata = data misfit for the state estimate
!     Jmod  = the model penalty function (the sum of the
!             dynamical, initial and boundary penalties)
!     Jb    = the ACTUAL model penalty function.
!     Jobs  = the ACTUAL data penalty function.
!     Jact  = the ACTUAL total penalty function.
!
!  NOTE:  Unlike IS4DVAR where we compute the actual cost function
!         for each iteration that decreases with increasing
!         number of iterations, the various penalty terms
!         represent the optimal values (i.e. the values when
!         optimal state estimate has been identified). Therefore
!         the various penalty terms will only be correct when the
!         optimal state has been found. Therefore, as the
!         W4DVAR proceeds, Jopt, Jdata and Jmod will asymptote
!         to their final values and WILL NOT necessary decrease
!         with increasing iteration number.
!
            IF (outLoop.eq.1) THEN
              DO iobs=1,Ndatum(ng)
                IF (ObsErr(iobs).ne.0.0_r8) THEN
!                  Jopt=Jopt+px(iobs)*vgrad0(iobs)/SQRT(ObsErr(iobs))
                  Jopt=Jopt+px(iobs)*TLmodVal_S(iobs,1,outLoop)
!                  Jf=Jf+vgrad0(iobs)*vgrad0(iobs)/ObsErr(iobs)
                  Jf=Jf+vgrad0(iobs)*TLmodVal_S(iobs,1,outLoop)
                END IF
!AMM: For Jdata we need GDG'*px but we don't have this except at the end
!AMM                Jdata=Jdata+px(iobs)*px(iobs)
              END DO
!AMM              Jmod=Jopt-Jdata
            ELSE
              DO iobs=1,Ndatum(ng)
!AMM                Jopt=Jopt-(px(iobs)+cg_pxsave(iobs))*                   &
!AMM     &                             cg_innov(iobs)
                IF (ObsErr(iobs).ne.0.0_r8) THEN
                  Jf=Jf+cg_innov(iobs)*cg_innov(iobs)
                END IF
!AMM                Jdata=Jdata+                                            &
!AMM     &                (px(iobs)+cg_pxsave(iobs))*                       &
!AMM     &                (px(iobs)+cg_pxsave(iobs))
              END DO
              Jmod=Jopt-Jdata
            END IF
!
!  Compute the actual penalty function terms.
!
            DO ivec=1,innLoop-1
              IF (ivec.eq.1) THEN
                zfact(ivec)=1.0_r8/cg_Gnorm_v(outLoop)
              ELSE
                zfact(ivec)=1.0_r8/cg_beta(ivec,outLoop)
              ENDIF
            END DO
!
!AUGMENTED.
!
            DO iobs=1,Ndatum(ng)+1
              zd(iobs)=vgrad0s(iobs)
            END DO
!
!AUGMENTED.
!
            DO iobs=1,Ndatum(ng)+1
              zhv(iobs)=zd(iobs)
            END DO
            DO ivec=1,innLoop-1
              ztemp1(ivec)=0.0_r8
              ztemp3(ivec)=0.0_r8
              ztemp4(ivec)=0.0_r8
              DO iobs=1,Ndatum(ng)
                ztemp1(ivec)=ztemp1(ivec)+                              &
     &                       TLmodVal_S(iobs,ivec,outLoop)*zhv(iobs)*   &
     &                       zfact(ivec)
                ztemp3(ivec)=ztemp3(ivec)+                              &
     &                       vcglwk(iobs,ivec,outLoop)*Hbk(iobs,outLoop)
                ztemp4(ivec)=ztemp4(ivec)+                              &
     &                       TLmodVal_S(iobs,ivec,outLoop)*             &
     &                       zgrad0(iobs,outLoop)*zfact(ivec)
              END DO
!
!AUGMENTED
!  Note the factor zfact is not needed in the next loop since it has
!  already been applied to vcglwk above.
!
              DO iobs=1,Ndatum(ng)
                ztemp1(ivec)=ztemp1(ivec)+zhv(iobs)*                    &
     &          Hbk(iobs,outLoop)*vcglwk(Ndatum(ng)+1,ivec,outLoop)+    &
     &          Hbk(iobs,outLoop)*vcglwk(iobs,ivec,outLoop)*            &
     &          zhv(Ndatum(ng)+1) 
                ztemp4(ivec)=ztemp4(ivec)+zgrad0(iobs,outLoop)*         &
     &          Hbk(iobs,outLoop)*vcglwk(Ndatum(ng)+1,ivec,outLoop)+    &
     &          Hbk(iobs,outLoop)*vcglwk(iobs,ivec,outLoop)*            &
     &          zgrad0(Ndatum(ng)+1,outLoop) 
              END DO
              ztemp1(ivec)=ztemp1(ivec)+zhv(Ndatum(ng)+1)*              &
     &              Jb0(outLoop-1)*vcglwk(Ndatum(ng)+1,ivec,outLoop)
              ztemp3(ivec)=ztemp3(ivec)+                                &
     &              Jb0(outLoop-1)*vcglwk(Ndatum(ng)+1,ivec,outLoop)
              ztemp4(ivec)=ztemp4(ivec)+zgrad0(Ndatum(ng)+1,outLoop)*   &
     &              Jb0(outLoop-1)*vcglwk(Ndatum(ng)+1,ivec,outLoop)
!TEST
              ztemp4(ivec)=0.0_r8
!TEST
            END DO
            zbet=cg_delta(1,outLoop)
!TEST            zu1(1)=ztemp1(1)/zbet
!TEST            zu1(1)=ztemp4(1)/zbet
            zu1(1)=cg_QG(1,outLoop)/zbet
!TEST
            DO ivec=2,innLoop-1
              zgam(ivec)=cg_beta(ivec,outLoop)/zbet
              zbet=cg_delta(ivec,outLoop)-                              &
     &             cg_beta(ivec,outLoop)*zgam(ivec)
!TEST              zu1(ivec)=(ztemp1(ivec)-                                  &
              zu1(ivec)=(ztemp4(ivec)-                                  &
     &                   cg_beta(ivec,outLoop)*zu1(ivec-1))/zbet
            END DO
            ztemp2(innLoop-1)=zu1(innLoop-1)

            DO ivec=innLoop-2,1,-1
              zu1(ivec)=zu1(ivec)-zgam(ivec+1)*zu1(ivec+1)
              ztemp2(ivec)=zu1(ivec)
            END DO
!
            Jb=Jb0(outLoop-1)
            Jact=Jb
            DO iobs=1,Ndatum(ng)
              IF (ObsErr(iobs).NE.0.0_r8) THEN
                Jact=Jact+ObsErr(iobs)*vgrad0s(iobs)*vgrad0s(iobs)
              END IF
            END DO
            DO ivec=1,innLoop-1
              Jb=Jb+ztemp2(ivec)*ztemp2(ivec)
!TEST              Jact=Jact-ztemp1(ivec)*ztemp2(ivec)+                      &
!TEST     &                 +2.0_r8*ztemp2(ivec)*ztemp3(ivec)
              Jact=Jact-ztemp1(ivec)*ztemp2(ivec)
            END DO
            DO iobs=1,Ndatum(ng)
              Jb=Jb-2.0_r8*px(iobs)*Hbk(iobs,outLoop)
!TEST
              Jact=Jact+2.0_r8*px(iobs)*Hbk(iobs,outLoop)
!TEST
            END DO
            Jb=Jb-2.0_r8*px(Ndatum(ng)+1)*Jb0(outLoop-1)
!TEST
            Jact=Jact+2.0_r8*px(Ndatum(ng)+1)*Jb0(outLoop-1)
!TEST
            Jobs=Jact-Jb

          END IF
!
!  NEXT CALCULATION
!
!  After the initialization, for every other inner loop, calculate a
!  new Lanczos vector, store it in the matrix, and update search.
!
!AUGMENTED.
!
          DO iobs=1,Ndatum(ng)
!
!  Selime's code: w=t1.
!
            pgrad(iobs)=TLmodVal(iobs)+                                 &
     &            Hbk(iobs,outLoop)*vcglwk(Ndatum(ng)+1,innLoop,outLoop)
          END DO
          DO iobs=1,Ndatum(ng)
!
!  Convert gradient contribution from x-space to v-space.
!
            IF (ObsErr(iobs).NE.0.0_r8) THEN
!
!  Selime's code: w=t1/R. No augmented term here since Raug^-1=[R^-1 0].
!
              pgrad(iobs)=pgrad(iobs)/ObsErr(iobs)
            END IF
          END DO
          pgrad(Ndatum(ng)+1)=0.0_r8
!
!  Selime's code: w=t1/R+z1
!
          DO iobs=1,Ndatum(ng)
            pgrad(iobs)=pgrad(iobs)+vcglwk(iobs,innLoop,outLoop)
          END DO
          pgrad(Ndatum(ng)+1)=pgrad(Ndatum(ng)+1)+                       &
     &                             vcglwk(Ndatum(ng)+1,innLoop,outLoop)
!
          IF (innLoop.gt.1) THEN
!
!  Selime's code: w=w-v0*beta
!
!AUGMENTED.
            DO iobs=1,Ndatum(ng)+1
              pgrad(iobs)=pgrad(iobs)-cg_beta(innLoop,outLoop)*         &
     &                                 zcglwk(iobs,innLoop-1,outLoop)
            END DO
          END IF
!
          cg_delta(innLoop,outLoop)=0.0_r8
!
!  Selime's code: alpha=w'*t1=cg_delta
!
!
!AUGMENTED.
          DO iobs=1,Ndatum(ng)
            cg_delta(innLoop,outLoop)=cg_delta(innLoop,outLoop)+        &
     &                                TLmodVal(iobs)*pgrad(iobs)+       &
     &       Hbk(iobs,outLoop)*vcglwk(Ndatum(ng)+1,innLoop,outLoop)*    &
     &              pgrad(iobs)+                                        &
     &              Hbk(iobs,outLoop)*vcglwk(iobs,innLoop,outLoop)*     &
     &              pgrad(Ndatum(ng)+1)
          END DO
          cg_delta(innLoop,outLoop)=cg_delta(innLoop,outLoop)+          &
     &            Jb0(outLoop-1)*vcglwk(Ndatum(ng)+1,innLoop,outLoop)*  &
     &              pgrad(Ndatum(ng)+1)
!
!  Exit, if not a positive definite algorithm.
!
          IF (cg_delta(innLoop,outLoop).le.0.0_r8) THEN
            WRITE (stdout,20) cg_delta(innLoop,outLoop), outLoop,       &
     &                        innLoop
            exit_flag=8
            GO TO 10
          END IF
!
!  Selime's code: w=w-alpha*v1
!
!
!AUGMENTED.
!
          DO iobs=1,Ndatum(ng)+1
            pgrad(iobs)=pgrad(iobs)-cg_delta(innLoop,outLoop)*          &
     &                  zcglwk(iobs,innLoop,outLoop)
          END DO
!
!  Orthonormalize against previous directions.
!
!  Identify the appropriate Lanczos vector normalization coefficient.
!
          DO ivec=1,innLoop
            IF (ivec.eq.1) THEN
              zfact(ivec)=1.0_r8/cg_Gnorm_v(outLoop)
            ELSE
              zfact(ivec)=1.0_r8/cg_beta(ivec,outLoop)
            ENDIF
          END DO
!
!  Selime's code: cg_dla=pgrad*HBH'*v1=pgrad*TLmodVal_S/zfact
!
!
!AUGMENTED
!
          DO ivec=innLoop,1,-1
            cg_dla(ivec,outLoop)=0.0_r8
            DO iobs=1,Ndatum(ng)
              cg_dla(ivec,outLoop)=cg_dla(ivec,outLoop)+pgrad(iobs)*    &
     &                             TLmodVal_S(iobs,ivec,outLoop)*       &
     &                             zfact(ivec)+pgrad(iobs)*             &
     &             Hbk(iobs,outLoop)*vcglwk(Ndatum(ng)+1,ivec,outLoop)+ &
     &             Hbk(iobs,outLoop)*vcglwk(iobs,ivec,outLoop)*         &
     &             pgrad(Ndatum(ng)+1)
            END DO
            cg_dla(ivec,outLoop)=cg_dla(ivec,outLoop)+                  &
     &                Jb0(outLoop-1)*vcglwk(Ndatum(ng)+1,ivec,outLoop)* &
     &                pgrad(Ndatum(ng)+1)
            DO iobs=1,Ndatum(ng)+1
              pgrad(iobs)=pgrad(iobs)-                                  &
     &                    cg_dla(ivec,outLoop)*                         &
     &                    vcglwk(iobs,ivec,outLoop)
            END DO
          END DO
!
!  Save the non-normalized Lanczos vector. zcglwk is used as temporary 
!  storage. The calculation will be completed at the start of the next
!  inner-loop.
!
!
!  Selime's code: w=zcglwk.
!
!AUGMENTED
!
          DO iobs=1,Ndatum(ng)+1
            zcglwk(iobs,innLoop+1,outLoop)=pgrad(iobs)
          END DO
!
!  If preconditioning, convert pgrad from v-space to y-space.
!
!
!  Selime's code: z1=Gw.
!
          IF (Lprecond.and.(outLoop.gt.1)) THEN
            Lscale=1                 ! Spectral LMP
            Laug=.TRUE.
            CALL rprecond (ng, Lscale, Laug, outLoop, NinnLoop-1, pgrad)
          END IF

!
!  Selime's code: z1=vcglwk.
!
!AUGMENTED.
!
          DO iobs=1,Ndatum(ng)+1
            vcglwk(iobs,innLoop+1,outLoop)=pgrad(iobs)
          END DO
!
!  Put the new trial solution into the adjoint vector for the next
!  loop or put the final solution into the adjoint vector on the
!  final inner-loop.
!
          IF (innLoop.eq.NinnLoop) THEN
!
!  Note: px is already in v-space so there is no need to convert
!        if preconditioning; cg_pxsave is also in v-space.
!
            DO iobs=1,Ndatum(ng)
              ADmodVal(iobs)=px(iobs)
            END DO
            DO iobs=1,Ndatum(ng)+1
              FOURDVAR(ng)%cg_pxsum(iobs)=FOURDVAR(ng)%cg_pxsum(iobs)+  &
     &                                    FOURDVAR(ng)%cg_pxsave(iobs)
              FOURDVAR(ng)%cg_pxsave(iobs)=px(iobs)
!AMM              FOURDVAR(ng)%cg_pxsum(iobs)=FOURDVAR(ng)%cg_pxsum(iobs)+  &
!AMM     &                                    px(iobs)
            END DO
          ELSE
!
!  Selime's code: z1=ADmodVal.
!
            DO iobs=1,Ndatum(ng)
              ADmodVal(iobs)=vcglwk(iobs,innLoop+1,outLoop)
            END DO
!
          END IF
!
        END IF MINIMIZE
!
!-----------------------------------------------------------------------
!  Report minimization progress.
!-----------------------------------------------------------------------
!
        IF (innLoop.gt.0) THEN
          WRITE (stdout,30) Ndatum(ng),                                 &
     &                      outLoop, innLoop-1, eps,                    &
     &                      outLoop, innLoop-1, preducy,                &
     &                      outLoop, innLoop-1, preducv,                &
     &                      outLoop, innLoop-1, Jf,                     &
     &                      outLoop, innLoop-1, Jdata,                  &
     &                      outLoop, innLoop-1, Jmod,                   &
     &                      outLoop, innLoop-1, Jopt,                   &
     &                      outLoop, innLoop-1, Jb,                     &
     &                      outLoop, innLoop-1, Jobs,                   &
     &                      outLoop, innLoop-1, Jact,                   &
     &                      outLoop, innLoop-1
        END IF
        DO ivec=1,innLoop
          WRITE (stdout,40) ivec, cg_delta(ivec,outLoop),               &
     &                      cg_beta(ivec,outLoop), zwork(ivec,3)
        END DO
        WRITE (stdout,50) innLoop+1, cg_beta(innLoop+1,outLoop)
!
!  If preconditioning, report the Ritz eigenvalues used.
!
        IF ((LhessianEV.or.Lprecond).and.(outLoop.gt.1)) THEN
          IF (NritzEV.gt.0) THEN
            WRITE (stdout,60) outLoop, innLoop, nConvRitz
            ic=0
            DO ivec=1,NinnLoop-1
              IF (ivec.gt.(NinnLoop-1-nConvRitz)) THEN
                string='         '
                ic=ic+1
                WRITE (stdout,70) ivec, cg_Ritz(ivec,outLoop-1),        &
     &                            cg_RitzErr(ivec,outLoop-1),           &
     &                            TRIM(ADJUSTL(string)),                &
     &                            'Used=', ic
              ELSE
                string='             '
                WRITE (stdout,80) ivec, cg_Ritz(ivec,outLoop-1),        &
     &                            cg_RitzErr(ivec,outLoop-1),           &
     &                            TRIM(ADJUSTL(string))
              END IF
            END DO
          ELSE
            WRITE (stdout,90) outLoop, innLoop, RitzMaxErr
            ic=0
            DO ivec=1,NinnLoop
              IF (cg_RitzErr(ivec,outLoop-1).le.RitzMaxErr) THEN
                string='converged'
                ic=ic+1
                WRITE (stdout,70) ivec, cg_Ritz(ivec,outLoop-1),        &
     &                            cg_RitzErr(ivec,outLoop-1),           &
     &                            TRIM(ADJUSTL(string)),                &
     &                            'Good=', ic
              ELSE
                string='not converged'
                WRITE (stdout,80) ivec, cg_Ritz(ivec,outLoop-1),        &
     &                            cg_RitzErr(ivec,outLoop-1),           &
     &                            TRIM(ADJUSTL(string))
              END IF
            END DO
          END IF
        END IF
!
!-----------------------------------------------------------------------
!  If last inner loop, innLoop = NinnLoop.
!-----------------------------------------------------------------------
!
!  Compute the eigenvalues and eigenvectors of the tridiagonal matrix.
!
        IF (innLoop.eq.NinnLoop) THEN
          IF (LhessianEV.or.Lprecond) THEN
            DO ivec=1,innLoop-1
              cg_Ritz(ivec,outLoop)=cg_delta(ivec,outLoop)
            END DO
            DO ivec=1,innLoop-2
              zwork(ivec,1)=cg_beta(ivec+1,outLoop)
            END DO
!
!  Use the LAPACK routine DSTEQR to compute the eigenvectors and
!  eigenvalues of the tridiagonal matrix. If applicable, the
!  eigenpairs are computed by master thread only.
!
            CALL DSTEQR ('I', innLoop-1, cg_Ritz(1,outLoop), zwork(1,1),&
     &                   zgv, Ninner, work, info)
            IF (info.ne.0) THEN
              PRINT *,' CONGRAD - Error in DSTEQR: info = ',info
              exit_flag=8
              GO TO 10
            END IF
!
            DO j=1,Ninner
              DO i=1,Ninner
                cg_zv(i,j,outLoop)=zgv(i,j)
              END DO
            END DO
!
!  Estimate the Ritz value error bounds.
!
            DO ivec=1,innLoop-1
              cg_RitzErr(ivec,outLoop)=ABS(cg_beta(innLoop,outLoop)*  &
     &                                 cg_zv(innLoop-1,ivec,outLoop))
            END DO
!
!  Check for exploding or negative Ritz values.
!
            DO ivec=1,innLoop-1
              IF (cg_RitzErr(ivec,outLoop).lt.0.0_r8) THEN
                PRINT *,'CONGRAD - negative Ritz value found.'
                exit_flag=8
                GO TO 10
              END IF
            END DO
!
!  Change the error bounds for the acceptable eigenvectors.
!
            RitzMaxErr=HevecErr
            DO ivec=1,NinnLoop-1
              cg_RitzErr(ivec,outLoop)=cg_RitzErr(ivec,outLoop)/        &
     &                                 cg_Ritz(NinnLoop-1,outLoop)
            END DO
!
!  Report Eigenvalues and their relative accuracy.
!
            WRITE (stdout,100) outLoop, InnLoop, RitzMaxErr
            ic=0
            DO ivec=1,NinnLoop-1
              IF (cg_RitzErr(ivec,outLoop).le.RitzMaxErr) THEN
                string='converged'
                ic=ic+1
                WRITE (stdout,70) ivec, cg_Ritz(ivec,outLoop),          &
     &                            cg_RitzErr(ivec,outLoop),             &
     &                            TRIM(ADJUSTL(string)),                &
     &                            'Good=', ic
              ELSE
                string='not converged'
                WRITE (stdout,80) ivec, cg_Ritz(ivec,outLoop),          &
     &                            cg_RitzErr(ivec,outLoop),             &
     &                            TRIM(ADJUSTL(string))
              END IF
            END DO
!
!  Calculate the converged eigenvectors (vcglev) of the [R_n + Cobs].
!
            CALL RPevecs (ng, outLoop, NinnLoop-1)
          END IF
        END IF

      END IF MASTER_THREAD
!
!-----------------------------------------------------------------------
!  Come here if not a possite definite algorithm.
!-----------------------------------------------------------------------
!
 10   CONTINUE
# ifdef DISTRIBUTE
      CALL mp_bcasti (ng, model, exit_flag)
# endif
      IF (exit_flag.ne.NoError) RETURN

# ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Broadcast new solution to other nodes.
!-----------------------------------------------------------------------
!
      CALL mp_bcastf (ng, model, ADmodVal)
#  if defined SENSITIVITY_4DVAR || \
      defined TL_W4DPSAS        || defined TL_W4DVAR || \
      defined W4DPSAS           || defined W4DVAR
!BSP       CALL mp_bcastf (ng, model, TLmodVal_S(:,:,outLoop))
!BSP       CALL mp_bcastf (ng, model, Hbk(:,outLoop))
!BSP       CALL mp_bcastf (ng, model, Jb0(outLoop-1))
#  endif
!BSP       CALL mp_bcasti (ng, model, info)
!BSP       CALL mp_bcastf (ng, model, cg_beta(:,outLoop))
!BSP       CALL mp_bcastf (ng, model, cg_QG(:,outLoop))
!BSP       CALL mp_bcastf (ng, model, cg_delta(:,outLoop))
!BSP       CALL mp_bcastf (ng, model, cg_dla(:,outLoop))
!BSP       CALL mp_bcastf (ng, model, cg_Gnorm_y(:))
!BSP       CALL mp_bcastf (ng, model, cg_Gnorm_v(:))
!BSP       CALL mp_bcastf (ng, model, FOURDVAR(ng)%cg_pxsave(:))
!AMM       CALL mp_bcastf (ng, model, cg_pxsave(:))
!BSP       CALL mp_bcastf (ng, model, cg_innov(:))
!BSP       CALL mp_bcastf (ng, model, zgrad0(:,outLoop))
!BSP       CALL mp_bcastf (ng, model, zcglwk(:,:,outLoop))
!BSP       CALL mp_bcastf (ng, model, vcglwk(:,:,outLoop))
      IF ((LhessianEV.or.Lprecond).and.(innLoop.eq.NinnLoop)) THEN
        CALL mp_bcastf (ng, model, cg_Ritz(:,outLoop))
        CALL mp_bcastf (ng, model, cg_RitzErr(:,outLoop))
        CALL mp_bcastf (ng, model, cg_zv(:,:,outLoop))
        CALL mp_bcastf (ng, model, vcglev(:,:,outLoop))
      END IF
# endif
!
!-----------------------------------------------------------------------
!  Write out conjugate gradient vectors into 4DVAR NetCDF file.
!-----------------------------------------------------------------------
!
# ifdef BGQC
!
!  Save the cost function values for background quality control decisions
!
      Jact_S(innLoop)=Jact
!
# endif
      IF (innLoop.gt.0) THEN
         CALL cg_write_rpcg (ng, model, innLoop, outLoop,               &
     &                  Jf, Jdata, Jmod, Jopt, Jb, Jobs, Jact,          &
     &                  preducv, preducy)
      END IF
!
 20   FORMAT (/,' CONGRAD - Fatal error, not possitive definite',       &
     &        ' algorithm:',/,                                          &
     &        /,11x,'cg_delta = ',1p,e15.8,0p,3x,'(',i3.3,', ',i3.3,')')
 30   FORMAT (/,' CONGRAD - Conjugate Gradient Information:',/,         &
     &        /,11x,'Ndatum     = ',i7.7,/,/,                           &
     &        1x,'(',i3.3,',',i3.3,'): ',                               &
     &        'Residual estimate,                      eps = ',         &
     &        1p,e14.7,/,1x,'(',i3.3,',',i3.3,'): ',                    &
     &        'Reduction in gradient norm,  Greduc y-space = ',         &
     &        1p,e14.7,/,1x,'(',i3.3,',',i3.3,'): ',                    &
     &        'Reduction in gradient norm,  Greduc v-space = ',         &
     &        1p,e14.7,/,1x,'(',i3.3,',',i3.3,'): ',                    &
     &        'First guess initial data misfit,         Jf = ',         &
     &        1p,e14.7,/,1x,'(',i3.3,',',i3.3,'): ',                    &
     &        'State estimate data misfit,           Jdata = ',         &
     &        1p,e14.7,/,1x,'(',i3.3,',',i3.3,'): ',                    &
     &        'Model penalty function,                Jmod = ',         &
     &        1p,e14.7,/,1x,'(',i3.3,',',i3.3,'): ',                    &
     &        'Optimal penalty function,              Jopt = ',         &
     &        1p,e14.7,/,/,1x,'(',i3.3,',',i3.3,'): ',                  &
     &        'Actual Model penalty function,           Jb = ',         &
     &        1p,e14.7,/,1x,'(',i3.3,',',i3.3,'): ',                    &
     &        'Actual data penalty function,          Jobs = ',         &
     &        1p,e14.7,/,/,1x,'(',i3.3,',',i3.3,'): ',                  &
     &        'Actual total penalty function,         Jact = ',         &
     &        1p,e14.7,/,/,1x,'(',i3.3,',',i3.3,'): ',                  &
     &        'Lanczos vectors - cg_delta, cg_beta, zwork:',/)
 40   FORMAT (6x,i3.3,4x,3(1p,e15.8,5x))
 50   FORMAT (6x,i3.3,24x,1p,e15.8)
 60   FORMAT (/,1x,'(',i3.3,',',i3.3,'): ',                             &
     &        'Ritz eigenvalues used in preconditioning, ',             &
     &        'nConvRitz = ',i3.3,/)
 70   FORMAT (6x,i3.3,2x,1p,e14.7,2x,1p,e14.7,2x,a,5x,'('a,i3.3,')')
 80   FORMAT (6x,i3.3,2x,1p,e14.7,2x,1p,e14.7,2x,a)
 90   FORMAT (/,1x,'(',i3.3,',',i3.3,'): ',                             &
     &        'Ritz eigenvalues used in preconditioning, ',             &
     &        'RitzMaxErr = ',1p,e12.5,/)
100   FORMAT (/,1x,'(',i3.3,',',i3.3,'): ',                             &
     &        'New Ritz eigenvalues and their accuracy,  ',             &
     &        'RitzMaxErr = ',1p,e12.5,/)

      RETURN
      END SUBROUTINE rpcg_lanczos

      SUBROUTINE RPevecs (ng, outLoop, NinnLoop)
!
!=======================================================================
!                                                                      !
!  This routine computes the converged eigenvectors of (R_n+Cobs).     !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_iounits
!
      implicit none
!
!  Imported variable declarations
!
      integer, intent(in) :: ng, outLoop, NinnLoop
!
!  Local variable declarations.
!
      integer :: iobs, ivec, jn, jk, ij, jkk, ingood

      real(r8) :: dla, zfact
      real(r8), dimension(Ndatum(ng)+1) :: zt
!
!-----------------------------------------------------------------------
!  Compute and store converged eigenvectors of (R_n+Cobs).
!-----------------------------------------------------------------------
!
!  Count and collect the converged eigenvalues.
!
      ingood=0
      DO ivec=NinnLoop,1,-1
        ingood=ingood+1
        Ritz(ingood)=cg_Ritz(ivec,outLoop)
      END DO
!
!  Calculate the converged eigenvectors of (R_n+Cobs).
!
      ij=0
      DO jn=NinnLoop,1,-1
        ij=ij+1
!AUGMENTED
        DO iobs=1,Ndatum(ng)
          vcglev(iobs,ij,outLoop)=0.0_r8
        END DO
        DO jk=1,NinnLoop
          DO iobs=1,Ndatum(ng)
            vcglev(iobs,ij,outLoop)=vcglev(iobs,ij,outLoop)+            &
     &                              (vcglwk(iobs,jk,outLoop)-           &
     &                               FOURDVAR(ng)%cg_pxsum(iobs)*       &
     &                             vcglwk(Ndatum(ng)+1,jk,outLoop))*    &
     &                              cg_zv(jk,jn,outLoop)
          END DO
        END DO

!
! Orthonormalize.
!
        DO jk=1,ij-1
          DO iobs=1,Ndatum(ng)
            zt(iobs)=0.0_r8
          END DO
          DO iobs=1,Ndatum(ng)
            DO jkk=1,NinnLoop
              IF (jkk.eq.1) THEN
                zfact=1.0_r8/cg_Gnorm_v(outLoop)
              ELSE
                zfact=1.0_r8/cg_beta(jkk,outLoop)
              ENDIF
              zt(iobs)=zt(iobs)+(TLmodVal_S(iobs,jkk,outLoop)*zfact+    &
     &                       Hbk(iobs,outLoop)*                         &
     &                       vcglwk(Ndatum(ng)+1,jkk,outLoop))*          &
     &                       cg_zv(jkk,jk,outLoop)
            END DO
          END DO
          dla=0.0_r8
          DO iobs=1,Ndatum(ng)
            dla=dla+zt(iobs)*vcglev(iobs,ij,outLoop)
          END DO
          DO iobs=1,Ndatum(ng)
            vcglev(iobs,ij,outLoop)=vcglev(iobs,ij,outLoop)-            &
     &                              dla*zt(iobs)
          END DO
        END DO
        dla=0.0_r8
!AUGMENTED
        DO iobs=1,Ndatum(ng)
          dla=dla+vcglev(iobs,ij,outLoop)*vcglev(iobs,ij,outLoop)
        END DO
!AUGMENTED+1
        DO iobs=1,Ndatum(ng)
          vcglev(iobs,ij,outLoop)=vcglev(iobs,ij,outLoop)/SQRT(dla)
        END DO
      END DO

      RETURN
      END SUBROUTINE RPevecs

      SUBROUTINE rprecond (ng, Lscale, Laug, outLoop, NinnLoop, py)
!
!=======================================================================
!                                                                      !
!  This routine is the preconditioner in observation space.            !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_scalars
      USE mod_parallel
      USE mod_fourdvar
!
      implicit none
!
!  Imported variable declarations
!
      logical, intent(in) :: Laug
      integer, intent(in) :: ng, Lscale, outLoop, NinnLoop
!AUGMENTED
      real(r8), intent(inout) :: py(Ndatum(ng)+1)
!
!  Local variable declarations.
!
      integer :: is, ie, inc, iss, i, jk
      integer :: nol, nols, nole, ninc
      integer :: ingood
      integer :: iobs, nvec, ivec

      real(r8) :: dla, dlar, fac, facritz, zfact
      real(r8), dimension(NinnLoop) :: zvect
      real(r8), dimension(Ndatum(ng)+1) :: zt
!
!  Set the do-loop indices for the sequential preconditioner
!  loop.
!
      nols=1
      nole=outLoop-1
      ninc=1
!
!  Apply the preconditioners derived from all previous outer-loops
!  sequentially.
!
      DO nol=nols,nole,ninc
!
!  Determine the number of Ritz vectors to use.
!  For Lritz=.TRUE., choose HvecErr to be larger.
!
        ingood=0
        DO i=1,NinnLoop
          IF (cg_RitzErr(i,nol).le.RitzMaxErr) THEN
            ingood=ingood+1
          END IF
        END DO
        IF (NritzEV.gt.0) THEN
          nConvRitz=NritzEV
          ingood=nConvRitz
        ELSE
          nConvRitz=ingood
        END IF
!
        IF (Lscale.gt.0) THEN
          is=1
          ie=ingood
          inc=1
        ELSE
          is=ingood
          ie=1
          inc=-1
        END IF
!
!  Lscale determines the form of the preconditioner:
!
!     1= LMP
!    -1= Inverse LMP
!
        DO nvec=is,ie,inc
          DO iobs=1,Ndatum(ng)
            zt(iobs)=0.0_r8
          END DO
          DO iobs=1,Ndatum(ng)
            DO jk=1,NinnLoop
              IF (jk.eq.1) THEN
                zfact=1.0_r8/cg_Gnorm_v(nol)
              ELSE
                zfact=1.0_r8/cg_beta(jk,nol)
              ENDIF
              zt(iobs)=zt(iobs)+(TLmodVal_S(iobs,jk,nol)*zfact+         &
     &                       Hbk(iobs,nol)*vcglwk(Ndatum(ng)+1,jk,nol))*&
     &                       cg_zv(jk,nvec,nol)
            END DO
          END DO
!
          dla=0.0_r8
          DO iobs=1,Ndatum(ng)
            dla=dla+zt(iobs)*py(iobs)
          END DO
          IF (Lritz) THEN
            dlar=0.0_r8
            DO iobs=1,Ndatum(ng)
              dlar=dlar+(TLmodVal_S(iobs,NinnLoop+1,nol)/               &
     &            cg_beta(NinnLoop+1,nol)+                              &
     &            Hbk(iobs,nol)*vcglwk(Ndatum(ng)+1,NinnLoop+1,nol))*   &
     &            py(iobs)
            END DO
          END IF
!
!  NOTE: Lscale=1 and Lscale=-1 cases are not complete yet.
!
          IF (Lscale.eq.-1) THEN
            fac=(cg_Ritz(NinnLoop+1-nvec,nol)-1.0_r8)*dla
          ELSE
            fac=(1.0_r8/cg_Ritz(NinnLoop+1-nvec,nol)-1.0_r8)*dla
          END IF
!
          DO iobs=1,Ndatum(ng)
             py(iobs)=py(iobs)+fac*vcglev(iobs,nvec,nol)
          END DO
!
          IF (Lritz) THEN
!NOTE: We still need the code for Lscale=-1.
           IF (Lscale.eq.1) THEN
            facritz=cg_beta(NinnLoop+1,nol)*                            &
     &              cg_zv(NinnLoop,NinnLoop+1-nvec,nol)/                &
     &              cg_Ritz(NinnLoop+1-nvec,nol)
            DO iobs=1,Ndatum(ng)
               py(iobs)=py(iobs)+facritz*facritz*vcglev(iobs,nvec,nol)  &
     &                     -facritz*dlar*(vcglwk(iobs,NinnLoop+1,nol)-  &
     &                    FOURDVAR(ng)%cg_pxsum(iobs)*                  &
     &                    vcglwk(Ndatum(ng)+1,NinnLoop+1,nol))
            END DO
            DO iobs=1,Ndatum(ng)
               py(iobs)=py(iobs)-facritz*dla*vcglev(iobs,nvec,nol)
            END DO
           END IF
          END IF
!
        END DO
      END DO
!
      IF (Laug) THEN
        DO iobs=1,Ndatum(ng)
           py(iobs)=py(iobs)+                                           &
     &              FOURDVAR(ng)%cg_Dpxsum(iobs)*py(Ndatum(ng)+1)
        END DO
      END IF

      END SUBROUTINE rprecond

      SUBROUTINE cg_write_rpcg (ng, model, innLoop, outLoop,            &
     &                     Jf, Jdata, Jmod, Jopt, Jb, Jobs, Jact,       &
     &                     preducv, preducy)
!
!=======================================================================
!                                                                      !
!  This routine writes conjugate gradient vectors into 4DVAR NetCDF    !
!  for restart purposes.                                               !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      Use mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations
!
      integer, intent(in) :: ng, model, innLoop, outLoop

      real(r8), intent(in) :: Jf, Jdata, Jmod, Jopt, preducv, preducy
      real(r8), intent(in) :: Jb, Jobs, Jact
!
!  Local variable declarations.
!
      integer :: status
!
      SourceFile='congrad.F, cg_write_rpcg'
!
!-----------------------------------------------------------------------
!  Write out conjugate gradient variables.
!-----------------------------------------------------------------------
!
!  Write out outer and inner iteration.
!
      CALL netcdf_put_ivar (ng, model, DAV(ng)%name, 'outer',           &
     &                      outer, (/0/), (/0/),                        &
     &                      ncid = DAV(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_ivar (ng, model, DAV(ng)%name, 'inner',           &
     &                      inner, (/0/), (/0/),                        &
     &                      ncid = DAV(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out number of converged Ritz eigenvalues.
!
      IF (innLoop.eq.Ninner) THEN
        CALL netcdf_put_ivar (ng, model, DAV(ng)%name, 'nConvRitz',     &
     &                        nConvRitz, (/outLoop/), (/1/),            &
     &                        ncid = DAV(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out converged Ritz eigenvalues.
!
      IF (innLoop.eq.Ninner) THEN
        CALL netcdf_put_fvar (ng, model, DAV(ng)%name, 'Ritz',          &
     &                        Ritz, (/1,outLoop/), (/nConvRitz,1/),     &
     &                        ncid = DAV(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out conjugate gradient norm.
!
      IF (innLoop.gt.0) THEN
        CALL netcdf_put_fvar (ng, model, DAV(ng)%name, 'cg_beta',       &
     &                        cg_beta, (/1,1/), (/Ninner+1,Nouter/),    &
     &                        ncid = DAV(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out Lanczos algorithm coefficients.
!
      IF (innLoop.gt.0) THEN
        CALL netcdf_put_fvar (ng, model, DAV(ng)%name, 'cg_delta',      &
     &                        cg_delta, (/1,1/), (/Ninner,Nouter/),     &
     &                        ncid = DAV(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write normalization coefficients for Lanczos vectos.
!
      IF (innLoop.gt.0) THEN
        CALL netcdf_put_fvar (ng, model, DAV(ng)%name, 'cg_dla',        &
     &                        cg_dla, (/1,1/), (/Ninner,Nouter/),       &
     &                        ncid = DAV(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Initial gradient normalization factor.
!
      IF (innLoop.eq.1) THEN
        CALL netcdf_put_fvar (ng, model, DAV(ng)%name, 'cg_Gnorm_y',    &
     &                        cg_Gnorm_y, (/1/), (/Nouter/),            &
     &                        ncid = DAV(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN

        CALL netcdf_put_fvar (ng, model, DAV(ng)%name, 'cg_Gnorm_v',    &
     &                        cg_Gnorm_v, (/1/), (/Nouter/),            &
     &                        ncid = DAV(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Lanczos vector normalization factor.
!
      IF (innLoop.gt.0) THEN
        CALL netcdf_put_fvar (ng, model, DAV(ng)%name, 'cg_QG',         &
     &                        cg_QG, (/1,1/), (/Ninner,Nouter/),        &
     &                        ncid = DAV(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Reduction in the gradient norm.
!
      IF (innLoop.gt.0) THEN
        CALL netcdf_put_fvar (ng, model, DAV(ng)%name, 'cg_Greduc_y',   &
     &                        preducy, (/innLoop,outLoop/), (/1,1/),    &
     &                        ncid = DAV(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN

        CALL netcdf_put_fvar (ng, model, DAV(ng)%name, 'cg_Greduc_v',   &
     &                        preducv, (/innLoop,outLoop/), (/1,1/),    &
     &                        ncid = DAV(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Eigenvalues of Lanczos recurrence relationship.
!
      IF (innLoop.gt.0) THEN
        CALL netcdf_put_fvar (ng, model, DAV(ng)%name, 'cg_Ritz',       &
     &                        cg_Ritz, (/1,1/), (/Ninner,Nouter/),      &
     &                        ncid = DAV(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Eigenvalues relative error.
!
      IF (innLoop.gt.0) THEN
        CALL netcdf_put_fvar (ng, model, DAV(ng)%name, 'cg_RitzErr',    &
     &                        cg_RitzErr, (/1,1/), (/Ninner,Nouter/),   &
     &                        ncid = DAV(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Eigenvectors of Lanczos recurrence relationship.
!
      IF (innLoop.gt.0) THEN
        CALL netcdf_put_fvar (ng, model, DAV(ng)%name, 'cg_zv',         &
     &                        cg_zv(:,:,outLoop),                       &
     &                        (/1,1,outLoop/), (/Ninner,Ninner,1/),     &
     &                        ncid = DAV(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  First guess initial data misfit.
!
      IF (innLoop.gt.0) THEN
        CALL netcdf_put_fvar (ng, model, DAV(ng)%name, 'Jf',            &
     &                        Jf, (/innLoop+1,outLoop/), (/1,1/),       &
     &                        ncid = DAV(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
!
!  State estimate data misfit.
!
        CALL netcdf_put_fvar (ng, model, DAV(ng)%name, 'Jdata',         &
     &                        Jdata, (/innLoop+1,outLoop/), (/1,1/),    &
     &                        ncid = DAV(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
!
!  Model penalty function.
!
        CALL netcdf_put_fvar (ng, model, DAV(ng)%name, 'Jmod',          &
     &                        Jmod, (/innLoop+1,outLoop/), (/1,1/),     &
     &                        ncid = DAV(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
!
!  Optimal penalty function.
!
        CALL netcdf_put_fvar (ng, model, DAV(ng)%name, 'Jopt',          &
     &                        Jopt, (/innLoop+1,outLoop/), (/1,1/),     &
     &                        ncid = DAV(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
!
!  Actual model penalty function.
!
        CALL netcdf_put_fvar (ng, model, DAV(ng)%name, 'Jb',            &
     &                        Jb, (/innLoop+1,outLoop/), (/1,1/),       &
     &                        ncid = DAV(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
!
!  Actual data penalty function.
!
        CALL netcdf_put_fvar (ng, model, DAV(ng)%name, 'Jobs',          &
     &                        Jobs, (/innLoop+1,outLoop/), (/1,1/),     &
     &                        ncid = DAV(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
!
!  Actual total penalty function.
!
        CALL netcdf_put_fvar (ng, model, DAV(ng)%name, 'Jact',          &
     &                        Jact, (/innLoop+1,outLoop/), (/1,1/),     &
     &                        ncid = DAV(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Initial gradient for minimization.
!
      IF (innLoop.eq.1) THEN
        CALL netcdf_put_fvar (ng, model, DAV(ng)%name, 'zgrad0',        &
     &                        zgrad0(:,outLoop),                        &
     &                        (/1,outLoop/), (/Ndatum(ng)+1,1/),        &
     &                        ncid = DAV(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
        CALL netcdf_put_fvar (ng, model, DAV(ng)%name, 'vgrad0',        &
     &                        vgrad0,                                   &
     &                        (/1/), (/Ndatum(ng)+1/),                  &
     &                        ncid = DAV(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
        CALL netcdf_put_fvar (ng, model, DAV(ng)%name, 'Hbk',           &
     &                        Hbk(:,outLoop),                           &
     &                        (/1,outLoop/), (/Ndatum(ng),1/),          &
     &                        ncid = DAV(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
        CALL netcdf_put_fvar (ng, model, DAV(ng)%name, 'Jb0',           &
     &                        Jb0(0:), (/1/), (/Nouter+1/),             &
     &                        ncid = DAV(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Lanczos vectors.
!
      IF (innLoop.gt.0) THEN
        CALL netcdf_put_fvar (ng, model, DAV(ng)%name, 'zcglwk',        &
     &                        zcglwk(:,innLoop,outLoop),                &
     &                      (/1,innLoop,outLoop/), (/Ndatum(ng)+1,1,1/),&
     &                      ncid = DAV(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
        CALL netcdf_put_fvar (ng, model, DAV(ng)%name, 'vcglwk',        &
     &                        vcglwk(:,innLoop,outLoop),                &
     &                      (/1,innLoop,outLoop/), (/Ndatum(ng)+1,1,1/),&
     &                      ncid = DAV(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Saved TLmodVal_S.
!
      CALL netcdf_put_fvar (ng, model, DAV(ng)%name, 'TLmodVal_S',      &
     &                      TLmodVal_S(:,innLoop,outLoop),              &
     &                      (/1,innLoop,outLoop/), (/Ndatum(ng),1,1/),  &
     &                      ncid = DAV(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Synchronize model/observation NetCDF file to disk.
!-----------------------------------------------------------------------
!
      CALL netcdf_sync (ng, model, DAV(ng)%name, DAV(ng)%ncid)

      RETURN
      END SUBROUTINE cg_write_rpcg

#else
      SUBROUTINE rpcg_lanczos
      RETURN
      END SUBROUTINE rpcg_lanczos
#endif
