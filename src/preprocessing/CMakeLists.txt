
add_executable(nalu_preprocess
  PreProcessingTask.cpp
  NDTW2D.cpp
  SamplingPlanes.cpp
  ABLFields.cpp
  ChannelFields.cpp
  HITFields.cpp
  RotateMesh.cpp
  BdyIOPlanes.cpp
  NestedRefinement.cpp
  TranslateMesh.cpp
  InflowHistory.cpp

  PreProcessDriver.cpp
  nalu_preprocess.cpp
  )

target_link_libraries(nalu_preprocess
  utilscore
  ${Trilinos_TPL_LIBRARIES}
  ${Trilinos_LIBRARIES}
  ${YAML_CPP_LIBRARIES}
  ${MPI_LIBRARIES})

target_compile_options(nalu_preprocess
  PRIVATE $<$<COMPILE_LANGUAGE:CXX>:${WIND_UTILS_WARN_FLAGS}>)

if(MPI_COMPILE_FLAGS)
  set_target_properties(nalu_preprocess PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif(MPI_COMPILE_FLAGS)

if(MPI_LINK_FLAGS)
  set_target_properties(nalu_preprocess PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif(MPI_LINK_FLAGS)

install(TARGETS nalu_preprocess
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

file(GLOB PREPROCESSING_HEADERS *.h)
install(FILES ${PREPROCESSING_HEADERS}
  DESTINATION include)
