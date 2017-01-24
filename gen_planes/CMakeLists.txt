
add_executable(gen_zplanes ZPlanes.cpp gen_planes.cpp)
target_link_libraries(gen_zplanes ${Trilinos_LIBRARIES} ${YAML_LIBRARY})
