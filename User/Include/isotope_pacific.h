/*
*******************************************************************************
** Copyright (c) 2002-2013 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Isotope Pacific Application
**
*/

#undef  AFT_EIGENMODES          /* Adjoint Finite Time Eigenmodes */
#undef  CORRELATION             /* Background-error Correlation Check */
#undef  FORCING_SV              /* Forcing Singular Vectors */
#undef  FT_EIGENMODES           /* Finite Time Eigenmodes */
#undef  IS4DVAR                 /* Incremental, strong constraint 4DVAR */
#define NLM_DRIVER              /* Nonlinear Basic State trajectory */
#undef  OPT_PERTURBATION        /* Optimal perturbations */
#undef  PICARD_TEST             /* Picard Iterations Test */
#undef  R_SYMMETRY              /* Representer Matrix Symmetry Test */
#undef  SANITY_CHECK            /* Sanity Check */
#undef  SO_SEMI                 /* Stochastic Optimals: Semi-norm */
#undef  TLM_CHECK               /* Tangent Linear Model Check */
#undef  W4DPSAS                 /* Weak constraint 4D-PSAS */
#undef  W4DVAR                  /* Weak constraint 4DVAR */
#undef  VERIFICATION            /* NL Observation Verification Driver */
#undef  NORMALIZATION           /* Background error Covariance Normalization */
#undef  AD_SENSITIVITY          /* Adjoint Sensitivity Driver */

/* Switch for debugging */
#undef DIAGNOSTIC

/*
**-----------------------------------------------------------------------------
**  Nonlinear basic state settings.
**-----------------------------------------------------------------------------
*/
#define SOLVE3D
#define UV_ADV
#define UV_COR
#define UV_QDRAG
#define UV_VIS2
#define MIX_S_UV
#define DJ_GRADPS
#define TS_U3HADVECTION
#define TS_C4VADVECTION
#define UV_U3HADVECTION
#define UV_C4VADVECTION
#define TS_DIF2
#define MIX_GEO_TS
#define NONLIN_EOS
#define VISC_GRID
#define DIFF_GRID
#define SALINITY
#define AVERAGES
#define MASKING
#define SPHERICAL
#define CURVGRID

#define SPONGE

#define PERFECT_RESTART
#define MY25_MIXING
#if defined MY25_MIXING
#  define KANTHA_CLAYSON
#  define N2S2_HORAVG
#endif

#define ANA_BSFLUX
#define ANA_BTFLUX
#define ANA_BPFLUX

#define BULK_FLUXES

#if defined BULK_FLUXES
#  define LONGWAVE_OUT
#  define SOLAR_SOURCE
#  define QCORRECTION
#  undef SRELAXATION
#  define SCORRECTION
#  undef EMINUSP
#endif

#define AVERAGES

#define ISOTOPE
#define ISOTOPE_FLUXES

#ifdef DIAGNOSTIC
#  define DIAGNOSTICS_BIO
#  define DIAGNOSTICS_UV
#  define DIAGNOSTICS_TS
#endif

#if defined AD_SENSITIVITY
#  define FORWARD_READ
#  define FORWARD_MIXING
#endif

#if defined TLM_DRIVER
#  define FORWARD_READ
#endif

#if defined VERIFICATION
#  undef ANA_INITIAL
#endif

#if defined HESSIAN
#  define CONVOLVE
#endif

#if defined W4DPSAS || defined W4DVAR
#  undef  AVERAGES
#  define FULL_GRID
#  undef ANA_INITIAL
#  undef VCONVOLUTION
#  define CONVOLVE
#  undef RTRACE
#  undef FORWARD_MIXING
#  undef FORWARD_RHS
#  define FORWARD_WRITE
#  define FORWARD_READ
#  define OUT_DOUBLE
#endif

