/*
** svn $Id$
*************************************************** Brian Powell ***
** Copyright (c) 2002-2011 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
************************************************************************
**                                                                    **
**  Defines microbial model parameters.                               **
**                                                                    **
************************************************************************
*/

!
!  Define Enterococcus biological model parameters.
!
      Vinfo( 1)='BioIter'
      Vinfo( 2)='number of iterations to achieve convergence'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='AttSWUV'
      Vinfo( 2)='UV light attenuation due to sea water'
      Vinfo( 3)='meter-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='AttSWBlue'
      Vinfo( 2)='Blue light attenuation due to sea water'
      Vinfo( 3)='meter-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='Ent_GrowthBlue'
      Vinfo( 2)='Growth of Enterococcus due to Blue-Light repair'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='Ent_DecayUV'
      Vinfo( 2)='Mortality of Enterococcus due to UV absorption'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='wEntero'
      Vinfo( 2)='enterococcus sinking rate'
      Vinfo( 3)='m day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='wVulnificusA'
      Vinfo( 2)='vibrio Vulnificus A sinking rate'
      Vinfo( 3)='m day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='wVulnificusB'
      Vinfo( 2)='vibrio Vulnificus B sinking rate'
      Vinfo( 3)='m day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='zVulnificusA'
      Vinfo( 2)='vibrio Vulnificus A mortality rate'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='zVulnificusB'
      Vinfo( 2)='vibrio Vulnificus B mortality rate'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='nVulnificusA_lag'
      Vinfo( 2)='vibrio Vulnificus A mortality rate time-lag'
      Vinfo( 3)='time-steps'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='nVulnificusB_lag'
      Vinfo( 2)='vibrio Vulnificus B mortality rate time-lag'
      Vinfo( 3)='time-steps'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
!
!  Vulnificus weight parameters.
!
        status=def_dim(ng, model, ncid, ncname, 'nVulnificusAFit',       &
     &                 25, wghtdimA)
        IF (exit_flag.ne.NoError) RETURN

        Vinfo( 1)='vulnificusA_weights'
        Vinfo( 2)='weights for growth fit function'
        status=def_var(ng, model, ncid, varid, NF_TYPE,                 &
     &                 1, (/wghtdimA/), Aval, Vinfo, ncname,             &
     &               SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN

        Vinfo( 1)='vulnificusA_temp'
        Vinfo( 2)='temperature for growth fit function'
        Vinfo( 3)='Celcius'
        status=def_var(ng, model, ncid, varid, NF_TYPE,                 &
     &                 1, (/wghtdimA/), Aval, Vinfo, ncname,             &
     &               SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN

        Vinfo( 1)='vulnificusA_salt'
        Vinfo( 2)='salt for growth fit function'
        Vinfo( 3)='ppt'
        status=def_var(ng, model, ncid, varid, NF_TYPE,                 &
     &                 1, (/wghtdimA/), Aval, Vinfo, ncname,             &
     &               SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
      
        status=def_dim(ng, model, ncid, ncname, 'nVulnificusBFit',       &
     &                 25, wghtdimB)
        IF (exit_flag.ne.NoError) RETURN

        Vinfo( 1)='vulnificusB_weights'
        Vinfo( 2)='weights for growth fit function'
        status=def_var(ng, model, ncid, varid, NF_TYPE,                 &
     &                 1, (/wghtdimB/), Aval, Vinfo, ncname,             &
     &               SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN

        Vinfo( 1)='vulnificusB_temp'
        Vinfo( 2)='temperature for growth fit function'
        Vinfo( 3)='Celcius'
        status=def_var(ng, model, ncid, varid, NF_TYPE,                 &
     &                 1, (/wghtdimB/), Aval, Vinfo, ncname,             &
     &               SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN

        Vinfo( 1)='vulnificusB_salt'
        Vinfo( 2)='salt for growth fit function'
        Vinfo( 3)='ppt'
        status=def_var(ng, model, ncid, varid, NF_TYPE,                 &
     &                 1, (/wghtdimB/), Aval, Vinfo, ncname,             &
     &               SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
