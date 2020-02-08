
set(STRUCT_GRID_SOURCES
  StructGrid.cpp
  )

if (HYPRE_FOUND)
  list(APPEND STRUCT_GRID_SOURCES
    HypreStructSolver.cpp)
endif()

add_library(struct_grid ${STRUCT_GRID_SOURCES})

target_compile_options(struct_grid
  PRIVATE $<$<COMPILE_LANGUAGE:CXX>:${WIND_UTILS_WARN_FLAGS}>)

target_link_libraries(struct_grid
  ${MPI_LIBARIES})

if (HYPRE_FOUND)
  target_link_libraries(struct_grid ${HYPRE_LIBRARIES})
endif()

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
