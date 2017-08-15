
cmake_minimum_required(VERSION 2.8.11)

project(NALU_UTILS CXX C Fortran)

option(ENABLE_WRFTONALU "Enable WRF to Nalu conversion utility" off)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

find_package(MPI REQUIRED)

set(CMAKE_PREFIX_PATH ${Trilinos_DIR} ${CMAKE_PREFIX_PATH})
set(CMAKE_PREFIX_PATH ${Trilinos_PATH} ${CMAKE_PREFIX_PATH})
find_package(Trilinos REQUIRED)

find_library(YAML_LIBRARIES NAMES libyaml-cpp.a PATHS ${YAML_DIR}/lib)
find_path(YAML_INCLUDES yaml.h PATHS ${YAML_DIR}/include/yaml-cpp)
include_directories(${YAML_INCLUDES}/..)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
include_directories(${MPI_INCLUDE_PATH})
include_directories(${Trilinos_INCLUDE_DIRS})
include_directories(${Trilinos_TPL_INCLUDE_DIRS})
include_directories(${CMAKE_SOURCE_DIR})

add_subdirectory(core)
add_subdirectory(preprocessing)

if (ENABLE_WRFTONALU)
  set (NETCDF_F77 "YES")
  find_package (NetCDF REQUIRED)

  include_directories(${NETCDF_INCLUDE_DIRS})
  add_subdirectory(wrftonalu)
endif (ENABLE_WRFTONALU)
