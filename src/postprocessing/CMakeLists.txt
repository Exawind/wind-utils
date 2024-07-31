
add_executable(nalu_postprocess
  PostProcessingTask.cpp
  PostProcessDriver.cpp
  ABLStatistics.cpp

  nalu_postprocess.cpp)

target_link_libraries(nalu_postprocess
  PUBLIC
  utilscore
  ${Trilinos_TPL_LIBRARIES}
  ${Trilinos_LIBRARIES}
  yaml-cpp
  ${MPI_LIBRARIES})
target_include_directories(nalu_postprocess SYSTEM PUBLIC ${YAML_CPP_INCLUDE_DIR})

target_compile_options(nalu_postprocess
  PRIVATE $<$<COMPILE_LANGUAGE:CXX>:${WIND_UTILS_WARN_FLAGS}>)

if(MPI_COMPILE_FLAGS)
  set_target_properties(nalu_postprocess PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif(MPI_COMPILE_FLAGS)

if(MPI_LINK_FLAGS)
  set_target_properties(nalu_postprocess PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif(MPI_LINK_FLAGS)

install(TARGETS nalu_postprocess
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

file(GLOB POSTPROCESSING_HEADERS *.h)
install(FILES ${POSTPROCESSING_HEADERS}
  DESTINATION include)
