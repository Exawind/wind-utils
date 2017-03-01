add_executable(wrftonalu
  module_constants.F90
  module_dm.F90
  module_exodus.F90
  module_ncderrcheck.F90
  module_str2int.F90
  module_utmdeg_converter.F90
  wrftonalu.F90
  )

target_link_libraries(wrftonalu
  ${NETCDF_F77_LIBRARIES})
