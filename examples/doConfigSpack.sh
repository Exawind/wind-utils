# Example script to configure NaluWindUtils using Spack

# Set the Spack compiler flavor
COMPILER=gcc

# Point to appropriate spack installed libraries or local installations
trilinos_install_dir=$(spack location -i nalu-trilinos %${COMPILER})
netcdf_install_dir=$(spack location -i netcdf-fortran %${COMPILER})
yaml_install_dir=$(spack location -i yaml-cpp %${COMPILER})

EXTRA_ARGS=$@

cmake \
  -DTrilinos_DIR:PATH=$trilinos_install_dir \
  -DNETCDF_INCLUDE_DIR:PATH=$netcdf_install_dir/include \
  -DNETCDF_LIBRARY:PATH=$netcdf_install_dir/lib/libnetcdff.a \
  -DYAML_DIR:PATH=$yaml_install_dir \
  -DCMAKE_BUILD_TYPE=Debug \
$EXTRA_ARGS \
../src
