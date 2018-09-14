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

cmake_minimum_required(VERSION 3.1)

project(NaluWindUtils CXX C Fortran)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")
set(CMAKE_PREFIX_PATH ${Trilinos_DIR} ${CMAKE_PREFIX_PATH})

option(ENABLE_WRFTONALU "Enable WRF to Nalu Fortran conversion utility" off)
option(ENABLE_SPHINX_DOCS "Enable Sphinx-based user/developer manual" off)
option(ENABLE_DOXYGEN_DOCS "Enable API documentation using Doxygen" off)
option(ENABLE_SPHINX_API_DOCS "Enable linking Doxygen API documentation to Sphinx" off)

find_package(MPI REQUIRED)
find_package(Trilinos REQUIRED)
find_package(YAML-CPP REQUIRED)
find_package(OpenMP)

if (CMAKE_VERSION VERSION_LESS "3.1")
  set(CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS}")
else()
  set(CMAKE_CXX_STANDARD 11)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif()

include_directories(SYSTEM ${MPI_INCLUDE_PATH})
include_directories(SYSTEM ${Trilinos_INCLUDE_DIRS})
include_directories(SYSTEM ${Trilinos_TPL_INCLUDE_DIRS})
include_directories(SYSTEM ${YAML_CPP_INCLUDE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

if (OPENMP_FOUND)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
endif()

add_subdirectory(src)

if (ENABLE_SPHINX_DOCS OR ENABLE_DOXYGEN_DOCS)
  add_subdirectory(doc)
endif()

include(CTest)
add_subdirectory(tests)

# Install library dependency information
install(EXPORT ${PROJECT_NAME}Libraries
  DESTINATION lib/cmake/${PROJECT_NAME}
  FILE ${PROJECT_NAME}Libraries.cmake)

# Create NaluWindUtils config so that other codes can find NaluWindUtils
include(CMakePackageConfigHelpers)

set(INCLUDE_INSTALL_DIR include/)
set(LIB_INSTALL_DIR lib/)

configure_package_config_file(
  cmake/${PROJECT_NAME}.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}.cmake
  INSTALL_DESTINATION lib/cmake/${PROJECT_NAME}
  PATH_VARS INCLUDE_INSTALL_DIR LIB_INSTALL_DIR)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}.cmake
  DESTINATION lib/cmake/${PROJECT_NAME})
