
add_executable(turbsim_netcdf
  TurbulenceFile.cpp
  WindSimFile.cpp

  turbsim_netcdf.cpp
  )

target_compile_options(turbsim_netcdf
  PRIVATE $<$<COMPILE_LANGUAGE:CXX>:${WIND_UTILS_WARN_FLAGS}>)

target_link_libraries(turbsim_netcdf
  PUBLIC
  utilscore
  ${Trilinos_TPL_LIBRARIES}
  ${Trilinos_LIBRARIES}
  ${yaml-cpp}
  ${MPI_LIBRARIES})
target_include_directories(turbsim_netcdf SYSTEM PUBLIC ${YAML_CPP_INCLUDE_DIR})

if(MPI_COMPILE_FLAGS)
  set_target_properties(turbsim_netcdf PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif(MPI_COMPILE_FLAGS)

if(MPI_LINK_FLAGS)
  set_target_properties(turbsim_netcdf PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif(MPI_LINK_FLAGS)

install(TARGETS turbsim_netcdf
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)
