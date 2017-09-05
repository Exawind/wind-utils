#!/bin/bash
# Example script to configure NaluWindUtils using Spack

# Make sure SPACK_ROOT (spack root directory) is an environment variable
if [ -z ${SPACK_ROOT+x} ]; then
    echo "Please set SPACK_ROOT variable. Exiting configure script.";
    exit 1
fi

# Set the Spack compiler flavor
COMPILER=gcc

# Spack executable
SPACK=${SPACK_ROOT}/bin/spack

# Point to appropriate spack installed libraries or local installations
trilinos_install_dir=$(${SPACK} location -i nalu-trilinos %${COMPILER})
yaml_install_dir=$(${SPACK} location -i yaml-cpp %${COMPILER})
netcdf_install_dir=$(${SPACK} location -i netcdf %${COMPILER})
netcdffortran_install_dir=$(${SPACK} location -i netcdf-fortran@4.4.3 %${COMPILER})

# Load necessary modules created by spack
module use ${SPACK_ROOT}/share/spack/modules/$(${SPACK} arch)
module load $(${SPACK} module find binutils %${COMPILER})
module load $(${SPACK} module find cmake %${COMPILER})
module load $(${SPACK} module find openmpi %${COMPILER})

# Clean before cmake configure
set +e
rm -rf CMakeFiles
rm -f CMakeCache.txt
set -e

EXTRA_ARGS=$@

cmake \
  -DTrilinos_DIR:PATH=$trilinos_install_dir \
  -DYAML_DIR:PATH=$yaml_install_dir \
  -DNETCDF_DIR:PATH=$netcdf_install_dir \
  -DNETCDF_F77_ROOT:PATH=$netcdffortran_install_dir \
  -DCMAKE_BUILD_TYPE=Debug \
$EXTRA_ARGS \
..
