
function(preproc_test testname)
  add_test(
    ${testname}
    ${CMAKE_BINARY_DIR}/src/preprocessing/nalu_preprocess
    -i ${CMAKE_CURRENT_SOURCE_DIR}/${testname}.yaml)
  set_tests_properties(${testname} PROPERTIES LABELS "preprocess")
endfunction(preproc_test)

preproc_test(init_abl_fields)
preproc_test(generate_planes)
preproc_test(bdy_io_mesh)
preproc_test(percept_refine)
preproc_test(rotate_mesh)
preproc_test(move_mesh)
preproc_test(channel_fields)
