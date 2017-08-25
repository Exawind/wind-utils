
add_library(nalumeshutils
  HexBlockMesh.cpp
  )
target_link_libraries(nalumeshutils utilscore
  ${Trilinos_LIBRARIES}
  ${YAML_LIBRARIES}
  ${MPI_LIBRARIES})

add_executable(abl_mesh
  abl_mesh.cpp)
target_link_libraries(abl_mesh
  nalumeshutils utilscore
  ${Trilinos_LIBRARIES}
  ${YAML_LIBRARIES}
  ${MPI_LIBRARIES})

if(MPI_COMPILE_FLAGS)
  set_target_properties(abl_mesh PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif(MPI_COMPILE_FLAGS)

if(MPI_LINK_FLAGS)
  set_target_properties(abl_mesh PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif(MPI_LINK_FLAGS)
