
add_library(utilscore
  CFDMesh.cpp)
target_link_libraries(utilscore ${Trilinos_LIBRARIES} ${YAML_LIBRARIES})

if(MPI_COMPILE_FLAGS)
  set_target_properties(utilscore PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif(MPI_COMPILE_FLAGS)

if(MPI_LINK_FLAGS)
  set_target_properties(utilscore PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif(MPI_LINK_FLAGS)
