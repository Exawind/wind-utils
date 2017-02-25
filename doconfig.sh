nalu_build_dir=${HOME}/nalu/
trilinos_install_dir=$nalu_build_dir/install/trilinos
yaml_install_dir=$nalu_build_dir/install

EXTRA_ARGS=$@

cmake \
  -DTrilinos_DIR:PATH=$trilinos_install_dir \
  -DYAML_DIR:PATH=$yaml_install_dir \
  -DENABLE_INSTALL:BOOL=OFF \
  -DCMAKE_BUILD_TYPE=RELEASE \
$EXTRA_ARGS \
../src
