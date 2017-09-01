
add_custom_target(docs)

if (ENABLE_SPHINX_DOCS)
  add_subdirectory(manual)
endif()

if (ENABLE_DOXYGEN_DOCS)
  add_subdirectory(doxygen)
endif()
