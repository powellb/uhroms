/*
** svn $Id$
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2011 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
************************************************************************
**                                                                    **
**  Writes microbial input parameters into                            **
**  output NetCDF files. It is included in routine "wrt_info.F".      **
**                                                                    **
************************************************************************
*/

!
!  Write out microbial biological model parameters.
!
      CALL netcdf_put_ivar (ng, model, ncname, 'BioIter',               &
     &                      BioIter(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'AttSWUV',               &
     &                      AttSWUV(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'AttSWBlue',             &
     &                      AttSWBlue(ng), (/0/), (/0/),                &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'Entero_uvd',            &
     &                      Ent_uvd(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'VulnificusA_uvd',       &
     &                      VulA_uvd(ng), (/0/), (/0/),                 &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'VulnificusB_uvd',       &
     &                      VulB_uvd(ng), (/0/), (/0/),                 &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'Entero_blug',           &
     &                      Ent_blug(ng), (/0/), (/0/),                 &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'VulnificusA_blug',      &
     &                      VulA_blug(ng), (/0/), (/0/),                &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'VulnificusB_blug',      &
     &                      VulB_blug(ng), (/0/), (/0/),                &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'wEntero',               &
     &                      wEntero(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'wVulnificusA',          &
     &                      wVulA(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'wVulnificusB',          &
     &                      wVulB(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'zVulnificusA',          &
     &                      zVulA(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'zVulnificusB',          &
     &                      zVulB(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_ivar (ng, model, ncname, 'nVulnificusA_win',      &
     &                      nVulA_win(ng), (/0/), (/0/),                &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_ivar (ng, model, ncname, 'nVulnificusB_win',      &
     &                      nVulB_win(ng), (/0/), (/0/),                &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_ivar (ng, model, ncname, 'nVulnificusA_lag',      &
     &                      nVulA_lag(ng), (/0/), (/0/),                &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_ivar (ng, model, ncname, 'nVulnificusB_lag',      &
     &                      nVulB_lag(ng), (/0/), (/0/),                &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN


!
!  Vulnificus Weights
!
      IF (nVulAWeights(ng).GT.0) THEN
        CALL netcdf_put_fvar (ng, model, ncname, 'vulnificusA_weights', &
     &                        vulAwght, (/1/), (/nVulAWeights/),        &
     &                        ncid = ncid)
        CALL netcdf_put_fvar (ng, model, ncname, 'vulnificusA_temp',    &
     &                        vulAtemp, (/1/), (/nVulAWeights/),        &
     &                        ncid = ncid)
        CALL netcdf_put_fvar (ng, model, ncname, 'vulnificusA_salt',    &
     &                        vulAsalt, (/1/), (/nVulAWeights/),        &
     &                        ncid = ncid)
      END IF
      IF (nVulBWeights(ng).GT.0) THEN
        CALL netcdf_put_fvar (ng, model, ncname, 'vulnificusB_weights', &
     &                        vulBwght, (/1/), (/nVulBWeights/),        &
     &                        ncid = ncid)
        CALL netcdf_put_fvar (ng, model, ncname, 'vulnificusB_temp',    &
     &                        vulBtemp, (/1/), (/nVulBWeights/),        &
     &                        ncid = ncid)
        CALL netcdf_put_fvar (ng, model, ncname, 'vulnificusB_salt',    &
     &                        vulBsalt, (/1/), (/nVulBWeights/),        &
     &                        ncid = ncid)
      END IF
