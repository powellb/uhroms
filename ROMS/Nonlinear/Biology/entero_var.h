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

#if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
    defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR   || \
    defined SO_SEMI

/*
**  Adjoint sensitivity state biological tracers.
*/

              CASE ('idTads(iEntero)')
                idTads(iNO3_)=varid
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

#ifdef TS_PSOURCE

/*
**  Biological tracers point Source/Sinks (river runoff).
*/


              CASE ('idRtrc(iEntero)')
                idRtrc(iEntero)=varid
#endif
