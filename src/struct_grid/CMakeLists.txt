
add_library(struct_grid
  StructGrid.cpp)

target_link_libraries(struct_grid
  ${MPI_LIBARIES})

if(MPI_COMPILE_FLAGS)
  set_target_properties(struct_grid PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif(MPI_COMPILE_FLAGS)

if(MPI_LINK_FLAGS)
  set_target_properties(struct_grid PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif(MPI_LINK_FLAGS)

install(TARGETS struct_grid
  EXPORT "${PROJECT_NAME}Libraries"
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)
