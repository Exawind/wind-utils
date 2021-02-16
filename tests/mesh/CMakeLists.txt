
function(mesh_test testname)
  add_test(
    ${testname}
    ${CMAKE_CURRENT_BINARY_DIR}/src/mesh/abl_mesh
      -i ${CMAKE_CURRENT_SOURCE_DIR}/${testname}.yaml)
  set_tests_properties(${testname} PROPERTIES LABELS "mesh")
endfunction(mesh_test)

mesh_test(abl_mesh_bbox)
mesh_test(abl_mesh_vertices)
mesh_test(channel_mesh)
mesh_test(user_specified_spacing)
mesh_test(plot3d_to_exo)
add_test(
  slice_mesh
  ${CMAKE_CURRENT_BINARY_DIR}/src/mesh/slice_mesh
  -i ${CMAKE_CURRENT_SOURCE_DIR}/slice_mesh.yaml)
