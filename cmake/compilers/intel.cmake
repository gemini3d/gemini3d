add_compile_options(
$<IF:$<BOOL:${WIN32}>,/QxHost,-xHost>
"$<$<COMPILE_LANGUAGE:Fortran>:$<IF:$<BOOL:${WIN32}>,/warn:declarations,-implicitnone>>"
$<$<COMPILE_LANGUAGE:Fortran>:-traceback>
$<$<COMPILE_LANG_AND_ID:Fortran,Intel>:$<IF:$<BOOL:${WIN32}>,/Qopenmp,-qopenmp>>
$<$<COMPILE_LANG_AND_ID:C,IntelLLVM>:-fiopenmp>
$<$<COMPILE_LANG_AND_ID:CXX,IntelLLVM>:-fiopenmp>
$<$<COMPILE_LANG_AND_ID:Fortran,IntelLLVM>:-fiopenmp>
$<$<COMPILE_LANGUAGE:Fortran>:-heap-arrays>
"$<$<COMPILE_LANGUAGE:Fortran>:$<IF:$<BOOL:${WIN32}>,/Qdiag-disable:5268;/Qdiag-disable:7712,-diag-disable 5268;-diag-disable 7712>>"
)


# -fiopenmp:
# undefined reference to `omp_get_max_threads'
# undefined reference to `__kmpc_global_thread_num' and more similar

# heap-arrays: avoid stack overflow

add_link_options(-qopenmp)
# undefined reference to `__kmpc_begin'


# --- IMPORTANT: bounds checking
add_compile_options("$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:-check>")
# -check is an alias for -check all
# --- IMPORTANT

add_compile_options("$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:-debug;-fpe0;-fp-stack-check>")
# -debug is an alias for -debug all

# Fortran 2018 standard too many false warnings
