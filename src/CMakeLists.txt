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

add_subdirectory(core)
add_subdirectory(preprocessing)
add_subdirectory(postprocessing)
add_subdirectory(mesh)

if (ENABLE_WRFTONALU)
  set (NETCDF_F77 "YES")
  find_package (NetCDF REQUIRED)

  include_directories(${NETCDF_INCLUDE_DIRS})
  add_subdirectory(wrftonalu)
endif (ENABLE_WRFTONALU)
