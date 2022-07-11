# Ref: https://gitlab.kitware.com/cmake/cmake/-/blob/master/Modules/ProcessorCount.cmake
# https://gitlab.kitware.com/search?group_id=415&project_id=541&repository_ref=master&scope=blobs&search=GetNumberOfPhysicalCPU&snippets=false

add_library(cpu_count OBJECT cpu.f90 cpu.cpp)
