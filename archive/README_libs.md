# Gemini required libraries

Note that using optimized BLAS and LAPACK routines may vastly improve code performance.
However, the instructions below assume you are downloading and compiling these libraries manually.

## BLAS

The source code can be
[downloaded](http://www.netlib.org/blas/)

It can be built by navigating into the BLAS directory and running:

    make clean
    make

results will be in `libblas.a` (same directory as source)

## LAPACK

The LAPACK source code can be downloaded from: http://www.netlib.org/lapack/

Make sure that the BLAS location in LAPACK's make.inc points to the BLAS library compiled above (or whatever optimized library you wish to use) then:

    make lib

Results will be in source directory as 'libtmglib.a’ and 'liblapack.a’.

The remainder require a functioning OpenMPI installation - ubuntu and most other linux distributions has repository packages for this…

## LAPACK95

LAPACK95 can be downloaded from:  http://www.netlib.org/lapack95/

To build LAPACK95:

```sh
cd SRC
make clean
make double
```

Make sure that the make.inc file has the correct path for your LAPACK libraries.
lapack95.a results.

### BLACS

The source code can be downloaded from:  http://www.netlib.org/blacs/

Make sure that the top-level BLACs directory in the Bmake.inc points to the right place then:

```sh
make mpi what=clean
make mpi
```

Also make sure the variables for MPI correspond to the right directory and filename

### SCALAPACK

Download source code from:  http://www.netlib.org/scalapack/

Go into the base directory and then:

```sh
cd SRC
make clean
make
```

Check SLmake.inc to make sure that the location of your LAPACK libraries are correct.

### MUMPS

The source code must be requested via
[webpage](http://mumps.enseeiht.fr)
```sh
make clean
make
```

Make sure the makefile points to the correct path for the scalapack libraries.
