
add_executable(boxturb
  BoxTurb.cpp
  BoxTurbIO.cpp

  boxturb_main.cpp)

target_link_libraries(
  boxturb
  struct_grid
  ${Trilinos_TPL_LIBRARIES}
  ${Trilinos_LIBRARIES}
  ${YAML_CPP_LIBRARIES}
  ${MPI_LIBRARIES})

add_executable(teststk
  test.cpp)
target_link_libraries(
  teststk
  ${Trilinos_TPL_LIBRARIES}
  ${Trilinos_LIBRARIES}
  ${YAML_CPP_LIBRARIES}
  ${MPI_LIBRARIES})

if(MPI_COMPILE_FLAGS)
  set_target_properties(boxturb PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif(MPI_COMPILE_FLAGS)

if(MPI_LINK_FLAGS)
  set_target_properties(boxturb PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif(MPI_LINK_FLAGS)
