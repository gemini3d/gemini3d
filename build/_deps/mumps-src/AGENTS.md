
## MUMPS < 5.9.1 Flang workaround

With MUMPS < 5.9.1, some MUMPS sources use `omp_lib` but don't guard it with `!$ sentinel`.

A few of the MUMPS source files have unconditional `use omp_lib` statements:

```
*fac_b.F
*fac_compact_factors_m.F
*fac_par_m.F
*sol_omp_m.F
mumps_mpitoomp_m.F
```

Gfortran ignores these unused `use omp_lib` statements, but Flang rightly does not ignore them and requires "omp_lib.mod" to be present, even if OpenMP isn't being used for parallelism.
It's not appropriate to unconditionally link OpenMP libraries, because that defeats MUMPS_openmp=no.
Thus we created "src/dummy_openmp.f90" that contains a dummy module "omp_lib" that satisfies the Flang syntax requirement without actually linking OpenMP libraries.

We used OBJECT library as INTERFACE libraries don't build the "omp_lib.mod" file.
