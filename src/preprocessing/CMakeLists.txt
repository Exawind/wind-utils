
add_executable(nalu_preprocess
  PreProcessingTask.cpp
  NDTW2D.cpp
  SamplingPlanes.cpp
  ABLFields.cpp
  RotateMesh.cpp

  PreProcessDriver.cpp
  nalu_preprocess.cpp
  )

target_link_libraries(nalu_preprocess
  utilscore ${Trilinos_LIBRARIES} ${YAML_LIBRARIES} ${MPI_LIBRARIES})

if(MPI_COMPILE_FLAGS)
  set_target_properties(nalu_preprocess PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif(MPI_COMPILE_FLAGS)

if(MPI_LINK_FLAGS)
  set_target_properties(nalu_preprocess PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif(MPI_LINK_FLAGS)
