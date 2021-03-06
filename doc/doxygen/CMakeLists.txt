
if (CMAKE_VERSION LESS "3.9")
  find_package(Doxygen REQUIRED)

  if (DOXYGEN_FOUND_DOT)
    set(DOXYGEN_HAVE_DOT "YES"
      CACHE STRING "Flag indicating whether Graphviz is installed")
  else()
    set(DOXYGEN_HAVE_DOT "NO"
      CACHE STRING "Flag indicating whether Graphviz is installed")
  endif()
else()
  cmake_policy(SET CMP0057 NEW)
  find_package(Doxygen
    OPTIONAL_COMPONENTS dot)
endif()

if (ENABLE_SPHINX_API_DOCS)
  set(DOXYGEN_GENERATE_XML "YES"
    CACHE STRING "Create XML tags for Sphinx Breathe")
else()
  set(DOXYGEN_GENERATE_XML "NO"
    CACHE STRING "Create XML tags for Sphinx Breathe")
endif()

configure_file(Doxyfile.in Doxyfile @ONLY)
add_custom_target(api-docs
  COMMAND ${DOXYGEN_EXECUTABLE} "${CMAKE_CURRENT_BINARY_DIR}/Doxyfile")

add_dependencies(docs api-docs)
