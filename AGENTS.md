* the `bind(C)` interface is a vital capability that couples to external programs written in C and C++.
* This project does not use `logical(C_BOOL)` at `bind(C)` interfaces due to possibility of bugs with distinct length integers (1 byte) vs. logical (4 bytes) in Fortran. Instead, we use `integer(C_INT)` and convert to and from logical in Fortran.
* It's perfectly fine to use Fortran `logical` inside this program, and even inside procedures that are `bind(C)`, as long as the interface to C uses `integer(C_INT)` and the conversion is handled correctly in Fortran.
* to convert from `integer(C_INT)` to `logical`, we check if the integer is not equal to zero. This is a common convention in C where zero represents false and any non-zero value represents true. For example:

```fortran
logical :: flag_fortran
integer(C_INT) :: flag_c
flag_fortran = flag_c /= 0
```
* to convert from `logical` to `integer(C_INT)`, we can use a simple conditional expression that returns 1 for true and 0 for false. This ensures that the values are correctly interpreted when passed to C functions. For example:

```fortran
logical :: flag_fortran
integer(C_INT) :: flag_c
flag_c = merge(1, 0, flag_fortran)
```
* examples of these interfaces are in src/libgemini_c.f90 `get_config_vars_C()` subroutine, where we convert the `flagneuBG` from `integer(C_INT)` to `logical` for use in Fortran, and then convert it back to `integer(C_INT)` before returning to C.


For reference, there were only a couple places that used `logical(C_BOOL)` in the original code.

* src/libgemini.f90:`type(c_params)` has members of type `integer(C_INT)` that represent boolean flags: `fortran_nml, fortran_cli, debug, dryrun`.
* src/libgemini_c.f90:`get_config_vars_C()` has a parameter `integer(C_INT) :: flagneuBG` used to receive the value from C.
