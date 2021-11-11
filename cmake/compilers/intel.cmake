add_compile_options(
$<IF:$<BOOL:${WIN32}>,/QxHost,-xHost>
"$<$<COMPILE_LANGUAGE:Fortran>:$<IF:$<BOOL:${WIN32}>,/warn:declarations,-implicitnone>>"
$<$<COMPILE_LANGUAGE:Fortran>:-traceback>
$<$<COMPILE_LANG_AND_ID:Fortran,Intel>:$<IF:$<BOOL:${WIN32}>,/Qopenmp,-qopenmp>>
$<$<COMPILE_LANG_AND_ID:C,IntelLLVM>:$<IF:$<BOOL:${WIN32}>,/Qiopenmp,-fiopenmp>>
$<$<COMPILE_LANG_AND_ID:CXX,IntelLLVM>:$<IF:$<BOOL:${WIN32}>,/Qiopenmp,-fiopenmp>>
$<$<COMPILE_LANG_AND_ID:Fortran,IntelLLVM>:$<IF:$<BOOL:${WIN32}>,/Qiopenmp,-fiopenmp>>
$<$<COMPILE_LANGUAGE:Fortran>:-heap-arrays>
$<$<COMPILE_LANGUAGE:Fortran>:$<IF:$<BOOL:${WIN32}>,/Qdiag-disable:5268$<COMMA>7712$<COMMA>10182,-diag-disable=5268$<COMMA>7712$<COMMA>10182>>
$<$<COMPILE_LANG_AND_ID:Fortran,IntelLLVM>:$<IF:$<BOOL:${WIN32}>,/Qdiag-disable:5415,-diag-disable=5415>>
)
# remark #5415: Feature not yet implemented: Some 'check' options temporarily disabled.
# warning #10182: disabling optimization; runtime debug checks enabled


# -fiopenmp:
# undefined reference to `omp_get_max_threads'
# undefined reference to `__kmpc_global_thread_num' and more similar

# heap-arrays: avoid stack overflow

add_link_options(-qopenmp)
# undefined reference to `__kmpc_begin'


# --- IMPORTANT: bounds checking
# add_compile_options("$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug,RelWithDebInfo>>:-check>")
# -check is an alias for -check all. However, MUMPS trips on -check, so we have to use a less stringent check.
add_compile_options("$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug,RelWithDebInfo>>:-CB>")
# -CB is an alias for -check bounds.
# --- IMPORTANT

add_compile_options("$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug,RelWithDebInfo>>:-debug>")
# -debug is an alias for -debug all
# -fpe0 causes MUMPS failures (internal to MUMPS)

# Fortran 2018 standard too many false warnings
