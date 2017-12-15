# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

# WRF to Nalu conversion utility requires NetCDF-Fortran library
set (NETCDF_F77 "YES")
find_package (NetCDF REQUIRED COMPONENTS F77)

include_directories(${NETCDF_INCLUDE_DIRS})

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
  ${Trilinos_TPL_LIBRARIES}
  ${Trilinos_LIBRARIES}
  ${NETCDF_F77_LIBRARIES}
  ${MPI_LIBRARIES})

install(TARGETS wrftonalu
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)
