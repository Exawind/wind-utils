
add_executable(boxturb
  BoxTurb.cpp
  BoxTurbIO.cpp
  boxturb_main.cpp)

target_link_libraries(
  boxturb
  PUBLIC
  struct_grid
  ${Trilinos_TPL_LIBRARIES}
  ${Trilinos_LIBRARIES}
  yaml-cpp
  ${MPI_LIBRARIES})
target_include_directories(boxturb SYSTEM PUBLIC ${YAML_CPP_INCLUDE_DIR})

if(MPI_COMPILE_FLAGS)
  set_target_properties(boxturb PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif(MPI_COMPILE_FLAGS)

if(MPI_LINK_FLAGS)
  set_target_properties(boxturb PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif(MPI_LINK_FLAGS)
