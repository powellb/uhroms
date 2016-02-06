      SUBROUTINE read_BioPar (model, inp, out, Lwrite)
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2011 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine reads in microbial ecosystem model input               !
!  parameters. They are specified in input script "microbe.in".        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_biology
      USE mod_ncparam
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations
!
      logical, intent(in) :: Lwrite
      integer, intent(in) :: model, inp, out
!
!  Local variable declarations.
!
      integer :: Npts, Nval, Nmax, i, itrc, ng, status

      integer :: decode_line, load_i, load_l, load_r

      logical, dimension(NBT,Ngrids) :: Ltrc

      real(r8), dimension(NBT,Ngrids) :: Rbio

      real(r8), dimension(100) :: Rval

      character (len=40 ) :: KeyWord
      character (len=256) :: line
      character (len=256), dimension(100) :: Cval
!
!-----------------------------------------------------------------------
!  Read in microbial model parameters.
!-----------------------------------------------------------------------
!
      IF (.not.allocated(BioIni)) allocate ( BioIni(MT,Ngrids) )
      DO WHILE (.TRUE.)
        READ (inp,'(a)',ERR=10,END=20) line
        status=decode_line(line, KeyWord, Nval, Cval, Rval)
        IF (status.gt.0) THEN
          IF (TRIM(KeyWord).eq.'Lbiology') THEN
            Npts=load_l(Nval, Cval, Ngrids, Lbiology)
          ELSE IF (TRIM(KeyWord).eq.'BioIter') THEN
            Npts=load_i(Nval, Rval, Ngrids, BioIter)
          ELSE IF (TRIM(KeyWord).eq.'BioIni(iEntero)') THEN
            Npts=load_r(Nval, Rval, Ngrids, BioIni(iEntero,1))
          ELSE IF (TRIM(KeyWord).eq.'BioIni(iVulA)') THEN
            Npts=load_r(Nval, Rval, Ngrids, BioIni(iVulA,1))
          ELSE IF (TRIM(KeyWord).eq.'BioIni(iVulB)') THEN
            Npts=load_r(Nval, Rval, Ngrids, BioIni(iVulB,1))
          ELSE IF (TRIM(KeyWord).eq.'PARfracUV') THEN
            Npts=load_r(Nval, Rval, Ngrids, PARfracUV)
          ELSE IF (TRIM(KeyWord).eq.'PARfracBlue') THEN
            Npts=load_r(Nval, Rval, Ngrids, PARfracBlue)
          ELSE IF (TRIM(KeyWord).eq.'AttSWUV') THEN
            Npts=load_r(Nval, Rval, Ngrids, AttSWUV)
          ELSE IF (TRIM(KeyWord).eq.'AttSWBlue') THEN
            Npts=load_r(Nval, Rval, Ngrids, AttSWBlue)
          ELSE IF (TRIM(KeyWord).eq.'Entero_uvd') THEN
            Npts=load_r(Nval, Rval, Ngrids, Ent_uvd)
          ELSE IF (TRIM(KeyWord).eq.'VulnificusA_uvd') THEN
            Npts=load_r(Nval, Rval, Ngrids, VulA_uvd)
          ELSE IF (TRIM(KeyWord).eq.'VulnificusB_uvd') THEN
            Npts=load_r(Nval, Rval, Ngrids, VulB_uvd)
          ELSE IF (TRIM(KeyWord).eq.'Entero_blug') THEN
            Npts=load_r(Nval, Rval, Ngrids, Ent_blug)
          ELSE IF (TRIM(KeyWord).eq.'VulnificusA_blug') THEN
            Npts=load_r(Nval, Rval, Ngrids, VulA_blug)
          ELSE IF (TRIM(KeyWord).eq.'VulnificusB_blug') THEN
            Npts=load_r(Nval, Rval, Ngrids, VulB_blug)
          ELSE IF (TRIM(KeyWord).eq.'wEntero') THEN
            Npts=load_r(Nval, Rval, Ngrids, wEntero)
          ELSE IF (TRIM(KeyWord).eq.'wVulnificusA') THEN
            Npts=load_r(Nval, Rval, Ngrids, wVulA)
          ELSE IF (TRIM(KeyWord).eq.'wVulnificusB') THEN
            Npts=load_r(Nval, Rval, Ngrids, wVulB)
          ELSE IF (TRIM(KeyWord).eq.'zVulnificusA') THEN
            Npts=load_r(Nval, Rval, Ngrids, zVulA)
          ELSE IF (TRIM(KeyWord).eq.'zVulnificusB') THEN
            Npts=load_r(Nval, Rval, Ngrids, zVulB)
          ELSE IF (TRIM(KeyWord).eq.'ccEntero') THEN
            Npts=load_r(Nval, Rval, Ngrids, ccEntero)
          ELSE IF (TRIM(KeyWord).eq.'crEntero') THEN
            Npts=load_r(Nval, Rval, Ngrids, crEntero)
          ELSE IF (TRIM(KeyWord).eq.'ccVulnificusA') THEN
            Npts=load_r(Nval, Rval, Ngrids, ccVulA)
          ELSE IF (TRIM(KeyWord).eq.'crVulnificusA') THEN
            Npts=load_r(Nval, Rval, Ngrids, crVulA)
          ELSE IF (TRIM(KeyWord).eq.'ccVulnificusB') THEN
            Npts=load_r(Nval, Rval, Ngrids, ccVulB)
          ELSE IF (TRIM(KeyWord).eq.'crVulnificusB') THEN
            Npts=load_r(Nval, Rval, Ngrids, crVulB)
          ELSE IF (TRIM(KeyWord).eq.'nVulnificusA_win') THEN
            Npts=load_i(Nval, Rval, Ngrids, nVulA_win)
            Npts=0
            Nmax=0
            DO ng=1,Ngrids
              nVulA_win(ng)=MAX(1,nVulA_win(ng))
              Npts=MAX(Npts,nVulA_win(ng))
              Nmax=MAX(Nmax,N(ng))
            END DO
            IF (.not.allocated(zVulA_win)) THEN
              allocate ( zVulA_win(Npts, Nmax, Ngrids) )
            END IF
            zVulA_win=0.0_r8
          ELSE IF (TRIM(KeyWord).eq.'nVulnificusB_win') THEN
            Npts=load_i(Nval, Rval, Ngrids, nVulB_win)
            Npts=0
            Nmax=0
            DO ng=1,Ngrids
              nVulB_win(ng)=MAX(1,nVulB_win(ng))
              Npts=MAX(Npts,nVulB_win(ng))
              Nmax=MAX(Nmax,N(ng))
            END DO
            IF (.not.allocated(zVulB_win)) THEN
              allocate ( zVulB_win(Npts, Nmax, Ngrids) )
            END IF
            zVulB_win=0.0_r8
          ELSE IF (TRIM(KeyWord).eq.'nVulnificusA_lag') THEN
            Npts=load_i(Nval, Rval, Ngrids, nVulA_lag)
            Npts=0
            Nmax=0
            DO ng=1,Ngrids
              nVulA_lag(ng)=MAX(1,nVulA_lag(ng))
              Npts=MAX(Npts,nVulA_lag(ng))
              Nmax=MAX(Nmax,N(ng))
            END DO
            IF (.not.allocated(zVulA_avg)) THEN
              allocate ( zVulA_avg(Npts, Nmax, Ngrids) )
            END IF
            zVulA_avg=0.0_r8
            IF (.not.allocated(zVulA_std)) THEN
              allocate ( zVulA_std(Npts, Nmax, Ngrids) )
            END IF
            zVulA_std=0.0_r8
            IF (.not.allocated(VulA_pop)) THEN
              allocate ( VulA_pop(Npts, Nmax, Ngrids) )
            END IF
            VulA_pop=0.0_r8
          ELSE IF (TRIM(KeyWord).eq.'nVulnificusB_lag') THEN
            Npts=load_i(Nval, Rval, Ngrids, nVulB_lag)
            Npts=0
            Nmax=0
            DO ng=1,Ngrids
              nVulB_lag(ng)=MAX(1,nVulB_lag(ng))
              Npts=MAX(Npts,nVulB_lag(ng))
              Nmax=MAX(Nmax,N(ng))
            END DO
            IF (.not.allocated(zVulB_avg)) THEN
              allocate ( zVulB_avg(Npts, Nmax, Ngrids) )
            END IF
            zVulB_avg=0.0_r8
            IF (.not.allocated(zVulB_std)) THEN
              allocate ( zVulB_std(Npts, Nmax, Ngrids) )
            END IF
            zVulB_std=0.0_r8
            IF (.not.allocated(VulB_pop)) THEN
              allocate ( VulB_pop(Npts, Nmax, Ngrids) )
            END IF
            VulB_pop=0.0_r8
          ELSE IF (TRIM(KeyWord).eq.'nVulAWeights') THEN
            Npts=load_i(Nval, Rval, Ngrids, nVulAWeights)
            DO ng=1,Ngrids
              nVulAWeights(ng)=MIN(40,nVulAWeights(ng))
            END DO
          ELSE IF (TRIM(KeyWord).eq.'vulnificusA_weights') THEN
            Npts=load_r(Nval, Rval, MAX(1,nVulAWeights), vulAwght)
          ELSE IF (TRIM(KeyWord).eq.'vulnificusA_temp') THEN
            Npts=load_r(Nval, Rval, MAX(1,nVulAWeights), vulAtemp)
          ELSE IF (TRIM(KeyWord).eq.'vulnificusA_salt') THEN
            Npts=load_r(Nval, Rval, MAX(1,nVulAWeights), vulAsalt)
          ELSE IF (TRIM(KeyWord).eq.'nVulBWeights') THEN
            Npts=load_i(Nval, Rval, Ngrids, nVulBWeights)
            DO ng=1,Ngrids
              nVulBWeights(ng)=MIN(40,nVulBWeights(ng))
            END DO
          ELSE IF (TRIM(KeyWord).eq.'vulnificusB_weights') THEN
            Npts=load_r(Nval, Rval, MAX(1,nVulBWeights), vulBwght)
          ELSE IF (TRIM(KeyWord).eq.'vulnificusB_temp') THEN
            Npts=load_r(Nval, Rval, MAX(1,nVulBWeights), vulBtemp)
          ELSE IF (TRIM(KeyWord).eq.'vulnificusB_salt') THEN
            Npts=load_r(Nval, Rval, MAX(1,nVulBWeights), vulBsalt)
          ELSE IF (TRIM(KeyWord).eq.'TNU2') THEN
            Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idbio(itrc)
                nl_tnu2(i,ng)=Rbio(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'TNU4') THEN
            Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idbio(itrc)
                nl_tnu4(i,ng)=Rbio(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'ad_TNU2') THEN
            Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idbio(itrc)
                ad_tnu2(i,ng)=Rbio(itrc,ng)
                tl_tnu2(i,ng)=Rbio(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'ad_TNU4') THEN
            Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idbio(itrc)
                ad_tnu4(i,ng)=Rbio(itrc,ng)
                ad_tnu4(i,ng)=Rbio(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'AKT_BAK') THEN
            Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idbio(itrc)
                Akt_bak(i,ng)=Rbio(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'ad_AKT_fac') THEN
            Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idbio(itrc)
                ad_Akt_fac(i,ng)=Rbio(itrc,ng)
                tl_Akt_fac(i,ng)=Rbio(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'TNUDG') THEN
            Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idbio(itrc)
                Tnudg(i,ng)=Rbio(itrc,ng)
              END DO
            END DO
#ifdef TS_PSOURCE
          ELSE IF (TRIM(KeyWord).eq.'LtracerSrc') THEN
            Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idbio(itrc)
                LtracerSrc(i,ng)=Ltrc(itrc,ng)
              END DO
            END DO
#endif
          ELSE IF (TRIM(KeyWord).eq.'Hout(idTvar)') THEN
            Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idTvar(idbio(itrc))
                IF (i.eq.0) THEN
                  IF (Master) WRITE (out,30)                            &
     &                              'idTvar(idbio(', itrc, '))'
                  exit_flag=5
                  RETURN
                END IF
                Hout(i,ng)=Ltrc(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Hout(idTsur)') THEN
            Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idTsur(idbio(itrc))
                IF (i.eq.0) THEN
                  IF (Master) WRITE (out,30)                            &
     &                              'idTsur(idbio(', itrc, '))'
                  exit_flag=5
                  RETURN
                END IF
                Hout(i,ng)=Ltrc(itrc,ng)
              END DO
            END DO
#ifdef AVERAGES
          ELSE IF (TRIM(KeyWord).eq.'Aout(idTvar)') THEN
            Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idTvar(idbio(itrc))
                Aout(i,ng)=Ltrc(itrc,ng)
              END DO
            END DO
#endif
#ifdef DIAGNOSTICS_TS
          ELSE IF (TRIM(KeyWord).eq.'Dout(iTrate)') THEN
            Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
            DO ng=1,Ngrids
              DO i=1,NBT
                itrc=idbio(i)
                Dout(idDtrc(itrc,iTrate),ng)=Ltrc(i,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Dout(iThadv)') THEN
            Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
            DO ng=1,Ngrids
              DO i=1,NBT
                itrc=idbio(i)
                Dout(idDtrc(itrc,iThadv),ng)=Ltrc(i,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Dout(iTxadv)') THEN
            Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
            DO ng=1,Ngrids
              DO i=1,NBT
                itrc=idbio(i)
                Dout(idDtrc(itrc,iTxadv),ng)=Ltrc(i,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Dout(iTyadv)') THEN
            Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
            DO ng=1,Ngrids
              DO i=1,NBT
                itrc=idbio(i)
                Dout(idDtrc(itrc,iTyadv),ng)=Ltrc(i,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Dout(iTvadv)') THEN
            Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
            DO ng=1,Ngrids
              DO i=1,NBT
                itrc=idbio(i)
                Dout(idDtrc(itrc,iTvadv),ng)=Ltrc(i,ng)
              END DO
            END DO
# if defined TS_DIF2 || defined TS_DIF4
          ELSE IF (TRIM(KeyWord).eq.'Dout(iThdif)') THEN
            Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
            DO ng=1,Ngrids
              DO i=1,NBT
                itrc=idbio(i)
                Dout(idDtrc(itrc,iThdif),ng)=Ltrc(i,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Dout(iTxdif)') THEN
            Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
            DO ng=1,Ngrids
              DO i=1,NBT
                itrc=idbio(i)
                Dout(idDtrc(itrc,iTxdif),ng)=Ltrc(i,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Dout(iTydif)') THEN
            Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
            DO ng=1,Ngrids
              DO i=1,NBT
                itrc=idbio(i)
                Dout(idDtrc(itrc,iTydif),ng)=Ltrc(i,ng)
              END DO
            END DO
#  if defined MIX_GEO_TS || defined MIX_ISO_TS
          ELSE IF (TRIM(KeyWord).eq.'Dout(iTsdif)') THEN
            Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
            DO ng=1,Ngrids
              DO i=1,NBT
                itrc=idbio(i)
                Dout(idDtrc(itrc,iTsdif),ng)=Ltrc(i,ng)
              END DO
            END DO
#  endif
# endif
          ELSE IF (TRIM(KeyWord).eq.'Dout(iTvdif)') THEN
            Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
            DO ng=1,Ngrids
              DO i=1,NBT
                itrc=idbio(i)
                Dout(idDtrc(itrc,iTvdif),ng)=Ltrc(i,ng)
              END DO
            END DO
#endif
          END IF
        END IF
      END DO
      DO ng=1,Ngrids
        zVulA_avg(:,:,ng)=0.0_r8
        zVulA_std(:,:,ng)=0.0_r8
        zVulB_avg(:,:,ng)=0.0_r8
        zVulB_std(:,:,ng)=0.0_r8
      END DO
      zVulA_std=1
      zVulB_std=1
  10  IF (Master) WRITE (out,40) line
      exit_flag=4
      RETURN
  20  CONTINUE
!
!-----------------------------------------------------------------------
!  Report input parameters.
!-----------------------------------------------------------------------
!
      IF (Lwrite) THEN
        DO ng=1,Ngrids
          IF (Lbiology(ng)) THEN
            WRITE (out,50) ng
            WRITE (out,60) BioIter(ng), 'BioIter',                      &
     &            'Number of iterations for nonlinear convergence.'
            WRITE (out,70) BioIni(iEntero,ng), 'BioIni(iEntero)',       &
     &            'Enterococcus initial concentration (nmol/m3).'
            WRITE (out,70) BioIni(iVulA,ng), 'BioIni(iVulA)',           &
     &            'Vulnificus A initial concentration (nmol/m3).'
            WRITE (out,70) BioIni(iVulB,ng), 'BioIni(iVulB)',           &
     &            'Vulnificus B initial concentration (nmol/m3).'
            WRITE (out,80) PARfracUV(ng), 'PARfracUV',                  &
     &            'Fraction of shortwave radiation that is',            &
     &            'UV active (nondimensional).'
            WRITE (out,80) PARfracBlue(ng), 'PARfracBlue',              &
     &            'Fraction of shortwave radiation that is',            &
     &            'Blue light active (nondimensional).'
            WRITE (out,70) AttSWUV(ng), 'AttSWUV',                      &
     &            'UV Light attenuation of seawater (m-1).'
            WRITE (out,70) AttSWBlue(ng), 'AttSWBlue',                  &
     &            'Blue Light attenuation of seawater (m-1).'
            WRITE (out,70) Ent_blug(ng), 'Entero_blug',                 &
     &            'Enterococcus growth due to Blue Light (day-1).'
            WRITE (out,70) VulA_blug(ng), 'VulA_blug',                  &
     &            'Vulnificus A growth due to Blue Light (day-1).'
            WRITE (out,70) VulB_blug(ng), 'VulB_blug',                  &
     &            'Vulnificus B growth due to Blue Light (day-1).'
            WRITE (out,70) Ent_uvd(ng), 'Entero_uvd',                   &
     &            'Enterococcus decay due to UV Absorption (day-1).'
            WRITE (out,70) VulA_uvd(ng), 'VulA_uvd',                    &
     &            'Vulnificus A decay due to UV Absorption (day-1).'
            WRITE (out,70) VulB_uvd(ng), 'VulB_uvd',                    &
     &            'Vulnificus B decay due to UV Absorption (day-1).'
            WRITE (out,70) wEntero(ng), 'wEntero',                      &
     &            'Enterococcus sinking rate (m/day).'
            WRITE (out,70) wVulA(ng), 'wVulnificusA',                   &
     &            'vibrio vulnificus A sinking rate (m/day).'
            WRITE (out,70) wVulB(ng), 'wVulnificusB',                   &
     &            'vibrio vulnificus B sinking rate (m/day).'
            WRITE (out,70) zVulA(ng), 'zVulnificusA',                   &
     &            'vibrio vulnificus A mortality rate (nmol/day).'
            WRITE (out,70) zVulB(ng), 'zVulnificusB',                   &
     &            'vibrio vulnificus B mortality rate (nmol/day).'
            WRITE (out,70) ccVulA(ng), 'ccVulnificusA',                 &
     &            'vibrio vulnificus A carrying capacity (nmol).'
            WRITE (out,70) crVulA(ng), 'crVulnificusA',                 &
     &            'vibrio vulnificus A carrying capacity growth ratio.'
            WRITE (out,70) ccVulB(ng), 'ccVulnificusB',                 &
     &            'vibrio vulnificus B carrying capacity (nmol).'
            WRITE (out,70) crVulB(ng), 'crVulnificusB',                 &
     &            'vibrio vulnificus B carrying capacity growth ratio.'
            WRITE (out,60) NvulA_win(ng), 'nVulnificusA_win',           &
     &         ' vibrio vulnificus A growth average window (steps).'
            WRITE (out,60) NvulB_win(ng), 'nVulnificusB_win',           &
     &         ' vibrio vulnificus B growth average window (steps).'
            WRITE (out,60) NvulA_lag(ng), 'nVulnificusA_lag',           &
     &         ' vibrio vulnificus A mortality lag from growth (steps).'
            WRITE (out,60) NvulB_lag(ng), 'nVulnificusB_lag',           &
     &         ' vibrio vulnificus B mortality lag from growth (steps).'
            IF (NvulAWeights(ng).GT.0) THEN
              WRITE  (out, 60) NvulAWeights(ng), 'NvulAWeights',        &
     &            'number of weights to use for microbes.'
              WRITE (out, *) 'VIBRIO VULNIFICUS A WEIGHTS'
              write (out, 250) 'TEMP',  'SALT', 'WEIGHT'
              DO i=1,NvulAWeights(ng)
                WRITE (out,251) vulAtemp(i), vulAsalt(i), vulAwght(i), i
              END DO
            END IF
            IF (NvulBWeights(ng).GT.0) THEN
              WRITE  (out, 60) NvulBWeights(ng), 'NvulBWeights',        &
     &            'number of weights to use for microbes.'
              WRITE (out, *) 'VIBRIO VULNIFICUS B WEIGHTS'
              write (out, 250) 'TEMP',  'SALT', 'WEIGHT'
              DO i=1,NvulBWeights(ng)
                WRITE (out,252) vulBtemp(i), vulBsalt(i), vulBwght(i), i
              END DO
            END IF
#ifdef TS_DIF2
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE (out,90) nl_tnu2(i,ng), 'nl_tnu2', i,               &
     &              'NLM Horizontal, harmonic mixing coefficient',      &
     &              '(m2/s) for tracer ', i, TRIM(Vname(1,idTvar(i)))
# ifdef ADJOINT
              WRITE (out,90) ad_tnu2(i,ng), 'ad_tnu2', i,               &
     &              'ADM Horizontal, harmonic mixing coefficient',      &
     &              '(m2/s) for tracer ', i, TRIM(Vname(1,idTvar(i)))
# endif
# if defined TANGENT || defined TL_IOMS
              WRITE (out,90) tl_tnu2(i,ng), 'tl_tnu2', i,               &
     &              'TLM Horizontal, harmonic mixing coefficient',      &
     &              '(m2/s) for tracer ', i, TRIM(Vname(1,idTvar(i)))
# endif
            END DO
#endif
#ifdef TS_DIF4
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE (out,90) nl_tnu4(i,ng), 'nl_tnu4', i,               &
     &              'NLM Horizontal, biharmonic mixing coefficient',    &
     &              '(m4/s) for tracer ', i, TRIM(Vname(1,idTvar(i)))
# ifdef ADJOINT
              WRITE (out,90) ad_tnu4(i,ng), 'ad_tnu4', i,               &
     &              'ADM Horizontal, biharmonic mixing coefficient',    &
     &              '(m4/s) for tracer ', i, TRIM(Vname(1,idTvar(i)))
# endif
# if defined TANGENT || defined TL_IOMS
              WRITE (out,90) tl_tnu4(i,ng), 'tl_tnu4', i,              &
     &              'TLM Horizontal, biharmonic mixing coefficient',    &
     &              '(m4/s) for tracer ', i, TRIM(Vname(1,idTvar(i)))
# endif
            END DO
#endif
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE(out,90) Akt_bak(i,ng), 'Akt_bak', i,                &
     &             'Background vertical mixing coefficient (m2/s)',     &
     &             'for tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
#ifdef FORWARD_MIXING
            DO itrc=1,NBT
              i=idbio(itrc)
# ifdef ADJOINT
              WRITE (out,90) ad_Akt_fac(i,ng), 'ad_Akt_fac', i,         &
     &              'ADM basic state vertical mixing scale factor',     &
     &              'for tracer ', i, TRIM(Vname(1,idTvar(i)))
# endif
# if defined TANGENT || defined TL_IOMS
              WRITE (out,90) tl_Akt_fac(i,ng), 'tl_Akt_fac', i,         &
     &              'TLM basic state vertical mixing scale factor',     &
     &              'for tracer ', i, TRIM(Vname(1,idTvar(i)))
# endif
            END DO
#endif
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE (out,90) Tnudg(i,ng), 'Tnudg', i,                   &
     &              'Nudging/relaxation time scale (days)',             &
     &              'for tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
#ifdef TS_PSOURCE
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE (out,100) LtracerSrc(i,ng), 'LtracerSrc',           &
     &              i, 'Processing point sources/Sink on tracer ', i,   &
     &              TRIM(Vname(1,idTvar(i)))
            END DO
#endif
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (Hout(idTvar(i),ng)) WRITE (out,110)                   &
     &            Hout(idTvar(i),ng), 'Hout(idTvar)',                   &
     &            'Write out tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (Hout(idTsur(i),ng)) WRITE (out,110)                   &
     &            Hout(idTsur(i),ng), 'Hout(idTsur)',                   &
     &            'Write out tracer flux ', i, TRIM(Vname(1,idTvar(i)))
            END DO
#ifdef AVERAGES
            WRITE (out,'(1x)')
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (Aout(idTvar(i),ng)) WRITE (out,110)                   &
     &            Aout(idTvar(i),ng), 'Aout(idTvar)',                   &
     &            'Write out averaged tracer ', i,                      &
     &            TRIM(Vname(1,idTvar(i)))
            END DO
#endif
#ifdef DIAGNOSTICS_TS
            WRITE (out,'(1x)')
            DO i=1,NBT
              itrc=idbio(i)
              IF (Dout(idDtrc(itrc,iTrate),ng))                         &
     &          WRITE (out,110) .TRUE., 'Dout(iTrate)',                 &
     &            'Write out rate of change of tracer ', itrc,          &
     &            TRIM(Vname(1,idTvar(itrc)))
            END DO
            DO i=1,NBT
              itrc=idbio(i)
              IF (Dout(idDtrc(itrc,iThadv),ng))                         &
     &          WRITE (out,110) .TRUE., 'Dout(iThadv)',                 &
     &            'Write out horizontal advection, tracer ', itrc,      &
     &            TRIM(Vname(1,idTvar(itrc)))
            END DO
            DO i=1,NBT
              itrc=idbio(i)
              IF (Dout(idDtrc(itrc,iTxadv),ng))                         &
     &          WRITE (out,110) .TRUE., 'Dout(iTxadv)',                 &
     &            'Write out horizontal X-advection, tracer ', itrc,    &
     &            TRIM(Vname(1,idTvar(itrc)))
            END DO
            DO i=1,NBT
              itrc=idbio(i)
              IF (Dout(idDtrc(itrc,iTyadv),ng))                         &
     &          WRITE (out,110) .TRUE., 'Dout(iTyadv)',                 &
     &            'Write out horizontal Y-advection, tracer ', itrc,    &
     &            TRIM(Vname(1,idTvar(itrc)))
            END DO
            DO i=1,NBT
              itrc=idbio(i)
              IF (Dout(idDtrc(itrc,iTvadv),ng))                         &
     &          WRITE (out,110) .TRUE., 'Dout(iTvadv)',                 &
     &            'Write out vertical advection, tracer ', itrc,        &
     &            TRIM(Vname(1,idTvar(itrc)))
            END DO
# if defined TS_DIF2 || defined TS_DIF4
            DO i=1,NBT
              itrc=idbio(i)
              IF (Dout(idDtrc(itrc,iThdif),ng))                         &
     &          WRITE (out,110) .TRUE., 'Dout(iThdif)',                 &
     &            'Write out horizontal diffusion, tracer ', itrc,      &
     &            TRIM(Vname(1,idTvar(itrc)))
            END DO
            DO i=1,NBT
              itrc=idbio(i)
              IF (Dout(idDtrc(i,iTxdif),ng))                            &
     &          WRITE (out,110) .TRUE., 'Dout(iTxdif)',                 &
     &            'Write out horizontal X-diffusion, tracer ', itrc,    &
     &            TRIM(Vname(1,idTvar(itrc)))
            END DO
            DO i=1,NBT
              itrc=idbio(i)
              IF (Dout(idDtrc(itrc,iTydif),ng))                         &
     &          WRITE (out,110) .TRUE., 'Dout(iTydif)',                 &
     &            'Write out horizontal Y-diffusion, tracer ', itrc,    &
     &            TRIM(Vname(1,idTvar(itrc)))
            END DO
#  if defined MIX_GEO_TS || defined MIX_ISO_TS
            DO i=1,NBT
              itrc=idbio(i)
              IF (Dout(idDtrc(itrc,iTsdif),ng))                         &
     &          WRITE (out,110) .TRUE., 'Dout(iTsdif)',                 &
     &            'Write out horizontal S-diffusion, tracer ', itrc,    &
     &            TRIM(Vname(1,idTvar(itrc)))
            END DO
#  endif
# endif
            DO i=1,NBT
              itrc=idbio(i)
              IF (Dout(idDtrc(itrc,iTvdif),ng))                         &
     &          WRITE (out,110) .TRUE., 'Dout(iTvdif)',                 &
     &            'Write out vertical diffusion, tracer ', itrc,        &
     &            TRIM(Vname(1,idTvar(itrc)))
            END DO
#endif
          END IF
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Rescale biological tracer parameters.
!-----------------------------------------------------------------------
!
!  Take the square root of the biharmonic coefficients so it can
!  be applied to each harmonic operator.
!
      DO ng=1,Ngrids
        DO itrc=1,NBT
          i=idbio(itrc)
          nl_tnu4(i,ng)=SQRT(ABS(nl_tnu4(i,ng)))
#ifdef ADJOINT
          ad_tnu4(i,ng)=SQRT(ABS(ad_tnu4(i,ng)))
#endif
#if defined TANGENT || defined TL_IOMS
          tl_tnu4(i,ng)=SQRT(ABS(tl_tnu4(i,ng)))
#endif
!
!  Compute inverse nudging coefficients (1/s) used in various tasks.
!
          IF (Tnudg(i,ng).gt.0.0_r8) THEN
            Tnudg(i,ng)=1.0_r8/(Tnudg(i,ng)*86400.0_r8)
          ELSE
            Tnudg(i,ng)=0.0_r8
          END IF
        END DO
      END DO

  30  FORMAT (/,' read_BioPar - variable info not yet loaded, ',        &
     &        a,i2.2,a)
  40  FORMAT (/,' read_BioPar - Error while processing line: ',/,a)
  50  FORMAT (/,/,' Microbial Model Parameters, Grid: ',i2.2,                &
     &        /,  ' ===============================',/)
  60  FORMAT (1x,i10,2x,a,t30,a)
  70  FORMAT (1p,e11.4,2x,a,t30,a)
  80  FORMAT (1p,e11.4,2x,a,t30,a,/,t32,a)
  90  FORMAT (1p,e11.4,2x,a,'(',i2.2,')',t30,a,/,t32,a,i2.2,':',1x,a)
 100  FORMAT (10x,l1,2x,a,'(',i2.2,')',t30,a,i2.2,':',1x,a)
 110  FORMAT (10x,l1,2x,a,t30,a,i2.2,':',1x,a)
 250  FORMAT (1p,a11,' ',a11,' ',a11)
 251  FORMAT (1p,e11.4,' ',e11.4,' ',e11.4,2x,'weight A ',i2.2,'.')
 252  FORMAT (1p,e11.4,' ',e11.4,' ',e11.4,2x,'weight B ',i2.2,'.')

      RETURN
      END SUBROUTINE read_BioPar
