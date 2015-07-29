/*
** svn $Id$
************************************************* Brian Powell, 2014 ***
** Copyright (c) 2002-2015 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
************************************************************************
**                                                                    **
**  Assigns metadata indices for the isotope model variables          **
**  that are used in input and output NetCDF files.                   **
**  The metadata information is read from "varinfo.dat".              **
**                                                                    **
**  This file is included in file "mod_ncparam.F", routine            **
**  "initialize_ncparm".                                              **
**                                                                    **
************************************************************************
*/

/*
**  Model state biological tracers.
*/

            CASE ('idTvar(i16O)')
              idTvar(i16O)=varid
            CASE ('idTvar(i18O)')
              idTvar(i18O)=varid

#if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
    defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR   || \
    defined SO_SEMI

/*
**  Adjoint sensitivity state biological tracers.
*/

            CASE ('idTads(i16O)')
              idTads(i16O)=varid
            CASE ('idTads(i18O)')
              idTads(i18O)=varid
#endif

/*
**  Biological tracers open boundary conditions.
*/

            CASE ('idTbry(iwest,i16O)')
              idTbry(iwest,i16O)=varid
            CASE ('idTbry(ieast,i16O)')
              idTbry(ieast,i16O)=varid
            CASE ('idTbry(isouth,i16O)')
              idTbry(isouth,i16O)=varid
            CASE ('idTbry(inorth,i16O)')
              idTbry(inorth,i16O)=varid

            CASE ('idTbry(iwest,i18O)')
              idTbry(iwest,i18O)=varid
            CASE ('idTbry(ieast,i18O)')
              idTbry(ieast,i18O)=varid
            CASE ('idTbry(isouth,i18O)')
              idTbry(isouth,i18O)=varid
            CASE ('idTbry(inorth,i18O)')
              idTbry(inorth,i18O)=varid

#ifdef TS_PSOURCE

/*
**  Biological tracers point Source/Sinks (river runoff).
*/


            CASE ('idRtrc(i16O)')
              idRtrc(i16O)=varid
            CASE ('idRtrc(i18O)')
              idRtrc(i18O)=varid
            CASE ('idRtrc(iVulB)')
              idRtrc(iVulB)=varid
#endif
