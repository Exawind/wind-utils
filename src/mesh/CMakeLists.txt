
# add_library(nalumeshutils SHARED
#   # Mesh spacing functions
#   spacing/MeshSpacing.cpp
#   spacing/ConstantSpacing.cpp
#   spacing/GeometricStretching.cpp

#   HexBlockBase.cpp
#   HexBlockMesh.cpp
#   Plot3DMesh.cpp
#   )
# target_link_libraries(nalumeshutils utilscore
#   ${Trilinos_TPL_LIBRARIES}
#   ${Trilinos_LIBRARIES}
#   ${YAML_CPP_LIBRARIES}
#   ${MPI_LIBRARIES})

add_executable(abl_mesh
  # Mesh spacing functions
  spacing/MeshSpacing.cpp
  spacing/ConstantSpacing.cpp
  spacing/GeometricStretching.cpp

  HexBlockBase.cpp
  HexBlockMesh.cpp
  Plot3DMesh.cpp
  abl_mesh.cpp)

add_executable(slice_mesh
  Slice.cpp
  slice_mesh.cpp
  )

target_link_libraries(abl_mesh
  struct_grid
  utilscore
  ${Trilinos_TPL_LIBRARIES}
  ${Trilinos_LIBRARIES}
  ${YAML_CPP_LIBRARIES}
  ${MPI_LIBRARIES})

target_link_libraries(slice_mesh
  utilscore
  ${Trilinos_TPL_LIBRARIES}
  ${Trilinos_LIBRARIES}
  ${YAML_CPP_LIBRARIES}
  ${MPI_LIBRARIES})

if(MPI_COMPILE_FLAGS)
  set_target_properties(abl_mesh PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif(MPI_COMPILE_FLAGS)

if(MPI_LINK_FLAGS)
  set_target_properties(abl_mesh PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif(MPI_LINK_FLAGS)

if(MPI_COMPILE_FLAGS)
  set_target_properties(slice_mesh PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif(MPI_COMPILE_FLAGS)

if(MPI_LINK_FLAGS)
  set_target_properties(slice_mesh PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif(MPI_LINK_FLAGS)

# install(TARGETS nalumeshutils
#  EXPORT "${PROJECT_NAME}Libraries"
#  RUNTIME DESTINATION bin
#  LIBRARY DESTINATION lib
#  ARCHIVE DESTINATION lib)

install(TARGETS abl_mesh slice_mesh
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

file(GLOB MESH_HEADERS *.h)
install(FILES ${MESH_HEADERS}
  DESTINATION include)
