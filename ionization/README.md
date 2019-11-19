## Glow ionization notes

Zxden is excitation of 12 species vs. altittude.
Zxden comes from Glow.
In Glow Zxden is a 2-D array, 12 x Nalt.

Nalt is zx1 in gemini.
In gemini, we strongly prefer to have zxden(lx1,lx2,lx3,12)
We will permute zxden immediately in Gemini and pass it through MPI.

A key factor to note is that zxden is `real32` like Glow, regardless of Gemini precision.
I kept zxden as real32 to save memory and in general to motivate generalizing Gemini to be real-precision polymorphic where appropriate.
`h5fortran` is real32/real64 polymorphic throughout.

## Fortran Glow procedure call order

(in reverse from Glow to Gemini)

1. procedure glow_run is a shim to Glow in a Fortran submodule in glow_run.F90
2. procedure ionrate_glow98 in ionization.f90 maps 2D glow dummy `zxden(12,Nalt)` to Gemini 4D actual `zxden(1:size(nn,1), 1:size(nn,2), 1:size(nn,3), 12)` with a nested loop. Note that `lx2==size(nn,2)` and `lx3==size(nn,3)`

  ```fortran
  call glow_run(..., zxden=zxden(:, ix2, ix3,:))
  ```
3. procedure fluid_adv in multifluid.f90 just passes through `zxden(:,:,:,:)`
4. program gemini in gemini.f90 does `allocate(zxden(lx1, lx2, lx3, 12))`, passing zxden from fluid_adv() to output_aur(). Remember, both MPI workers and MPI root do this, there is an `if` statement in output_aur() that selects gather_send() or gather_recv() respectively.
5. procedure output_aur() in a submodule of io.f90 in aurora.f90 consists simply of an `if` statement--root calls `output_aur_root()` and workers call `output_aur_workers()`
6. the MPI root goes to output_aur_root() in aurora_hdf5.f90. This module collects all the workers' zxden using generic procedure `gather_recv() => gather_recv32_2D_23()` in mpirecv32.f90 and writes to disk in HDF5 format. gather_recv32_2D_23(val, out, valall) takes in buffer array val w/o ghost cells and returns valall with ghost cells.
7. the MPI workers go to output_aur_workers() in aurora.f90. There, generic procedure `gather_send() => gather_send32_2D_23() in mpisend32.f90.