
add_library(utilscore
  CFDMesh.cpp)
target_link_libraries(utilscore ${Trilinos_LIBRARIES} ${YAML_LIBRARIES})
