
enable_testing()
execute_process(
  COMMAND ${CMAKE_COMMAND} -E create_symlink
  ${CMAKE_CURRENT_SOURCE_DIR}/grids
  ${CMAKE_CURRENT_BINARY_DIR}/grids)

add_subdirectory(mesh)
add_subdirectory(preprocessing)
