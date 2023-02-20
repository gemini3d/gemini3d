set_property(DIRECTORY PROPERTY LABELS diffusion)

add_library(PDEparabolic OBJECT PDEparabolic.f90)
target_link_libraries(PDEparabolic PRIVATE const gbsv)

add_library(diffusion OBJECT diffusion.f90
$<TARGET_OBJECTS:grid>
)
target_link_libraries(diffusion PRIVATE const meshobj PDEparabolic)
