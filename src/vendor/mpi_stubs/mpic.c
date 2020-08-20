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
#include "mpi.h"
LIBSEQ_INT LIBSEQ_CALL MPI_Init(LIBSEQ_INT *pargc, char ***pargv)
{
  return 0;
}

LIBSEQ_INT LIBSEQ_CALL MPI_Comm_rank( MPI_Comm comm, LIBSEQ_INT *rank)
{
  *rank=0;
  return 0;
}
LIBSEQ_INT LIBSEQ_CALL MPI_Finalize(void)
{
   return 0;
}

/* Internal: for MPI_IS_IN_PLACE tests from Fortran */

void LIBSEQ_CALL MUMPS_CHECKADDREQUAL(char *a, char*b, LIBSEQ_INT *i)
{
  if (a - b == 0)
   {
     *i=1;
   }
 else
   {
     *i=0;
   }
}

void LIBSEQ_CALL MUMPS_CHECKADDREQUAL_(char *a, char*b, LIBSEQ_INT *i)
 {
   MUMPS_CHECKADDREQUAL(a,b,i);
 }
void LIBSEQ_CALL mumps_checkaddrequal_(char *a, char*b, LIBSEQ_INT *i)
 {
   MUMPS_CHECKADDREQUAL(a,b,i);
 }
void LIBSEQ_CALL mumps_checkaddrequal__(char *a, char*b, LIBSEQ_INT *i)
 {
   MUMPS_CHECKADDREQUAL(a,b,i);
 }
