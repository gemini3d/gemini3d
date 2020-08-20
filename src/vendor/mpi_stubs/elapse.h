/*
 *
 *  This file is part of MUMPS 5.3.3, released
 *  on Mon Jun 15 09:57:25 UTC 2020
 *
 *
 *  Copyright 1991-2020 CERFACS, CNRS, ENS Lyon, INP Toulouse, Inria,
 *  Mumps Technologies, University of Bordeaux.
 *
 *  This version of MUMPS is provided to you free of charge. It is
 *  released under the CeCILL-C license:
 *  https://cecill.info/licences/Licence_CeCILL-C_V1-en.html
 *
 */

#ifndef MUMPS_CALL
#if defined(_WIN32)
/* Modify/choose between next 2 lines depending
 *  * on your Windows calling conventions */
/* #define MUMPS_CALL __stdcall */
#define MUMPS_CALL
#else
#define MUMPS_CALL
#endif
#endif

#if (defined(_WIN32) && ! defined(__MINGW32__)) || defined(UPPER)
#define mumps_elapse MUMPS_ELAPSE
#elif defined(Add__)
#define mumps_elapse mumps_elapse__
#elif defined(Add_)
#define mumps_elapse mumps_elapse_
#endif

void MUMPS_CALL mumps_elapse(double *val);
