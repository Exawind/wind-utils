
add_executable(calc_ndtw2d NDTW2D.cpp calc_ndtw2d.cpp)
target_link_libraries(calc_ndtw2d ${Trilinos_LIBRARIES} ${YAML_LIBRARY})
