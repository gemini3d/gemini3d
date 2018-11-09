
add_library(msis vendor/msis00/msis00_gfortran.f)
target_compile_options(msis PRIVATE -std=legacy -w -fno-implicit-none)

# --- for setting up an equilibrium simulation --

add_executable(msis_setup setup/MSIS00/call_msis_gfortran.f90)
target_link_libraries(msis_setup PUBLIC msis)

