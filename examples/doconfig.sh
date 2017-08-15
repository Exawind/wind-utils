nalu_build_dir=${HOME}/nalu/
trilinos_install_dir=$nalu_build_dir/install/trilinos
netcdf_install_dir=$nalu_build_dir/install
yaml_install_dir=$nalu_build_dir/install

EXTRA_ARGS=$@

# Cleanup old cache before we configure
# Note:  This does not remove files produced by make.  Use "make clean" for this.
#find . -name "CMakeFiles" -exec rm -rf {} \;
#rm -f CMakeCache.txt

cmake \
  -DTrilinos_DIR:PATH=$trilinos_install_dir \
  -DYAML_DIR:PATH=$yaml_install_dir \
  -DCMAKE_BUILD_TYPE=RELEASE \
  -DENABLE_WRFTONALU:BOOL=OFF \
  -DNetCDF_INCLUDE_DIR:PATH=$netcdf_install_dir/include \
  -DNetCDF_LIBRARY:PATH=$netcdf_install_dir/lib \
$EXTRA_ARGS \
../src
