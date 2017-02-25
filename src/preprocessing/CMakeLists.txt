
add_executable(nalu_preprocess
  PreProcessingTask.cpp
  NDTW2D.cpp
  SamplingPlanes.cpp
  ABLFields.cpp

  PreProcessDriver.cpp
  nalu_preprocess.cpp
  )

target_link_libraries(nalu_preprocess
  utilscore ${Trilinos_LIBRARIES} ${YAML_LIBRARIES} ${MPI_LIBRARIES})
