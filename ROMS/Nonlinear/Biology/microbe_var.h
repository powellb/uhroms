/*
** svn $Id$
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2011 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
************************************************************************
**
**  Enterococcus                                                      **
**                                                                    **
************************************************************************
*/

/*
**  Model state biological tracers.
*/

              CASE ('idTvar(iEntero)')
                idTvar(iEntero)=varid
              CASE ('idTvar(iVulA)')
                idTvar(iVulA)=varid
              CASE ('idTvar(iVulB)')
                idTvar(iVulB)=varid

#if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
    defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR   || \
    defined SO_SEMI

/*
**  Adjoint sensitivity state biological tracers.
*/

              CASE ('idTads(iEntero)')
                idTads(iEntero)=varid
              CASE ('idTads(iVulA)')
                idTads(iVulA)=varid
              CASE ('idTads(iVulB)')
                idTads(iVulB)=varid
#endif

/*
**  Biological tracers open boundary conditions.
*/

              CASE ('idTbry(iwest,iEntero)')
                idTbry(iwest,iEntero)=varid
              CASE ('idTbry(ieast,iEntero)')
                idTbry(ieast,iEntero)=varid
              CASE ('idTbry(isouth,iEntero)')
                idTbry(isouth,iEntero)=varid
              CASE ('idTbry(inorth,iEntero)')
                idTbry(inorth,iEntero)=varid

              CASE ('idTbry(iwest,iVulA)')
                idTbry(iwest,iVulA)=varid
              CASE ('idTbry(ieast,iVulA)')
                idTbry(ieast,iVulA)=varid
              CASE ('idTbry(isouth,iVulA)')
                idTbry(isouth,iVulA)=varid
              CASE ('idTbry(inorth,iVulA)')
                idTbry(inorth,iVulA)=varid

              CASE ('idTbry(iwest,iVulB)')
                idTbry(iwest,iVulB)=varid
              CASE ('idTbry(ieast,iVulB)')
                idTbry(ieast,iVulB)=varid
              CASE ('idTbry(isouth,iVulB)')
                idTbry(isouth,iVulB)=varid
              CASE ('idTbry(inorth,iVulB)')
                idTbry(inorth,iVulB)=varid

#ifdef TS_PSOURCE

/*
**  Biological tracers point Source/Sinks (river runoff).
*/


              CASE ('idRtrc(iEntero)')
                idRtrc(iEntero)=varid
              CASE ('idRtrc(iVulA)')
                idRtrc(iVulA)=varid
              CASE ('idRtrc(iVulB)')
                idRtrc(iVulB)=varid
#endif
