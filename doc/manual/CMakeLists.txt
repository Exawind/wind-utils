
find_program(SPHINX_EXE NAMES sphinx-build
  DOC "Sphinx Documentation Generator")

if (NOT SPHINX_EXE)
  message(FATAL "Cannot find Sphinx program to generate documentation. Disable doc generation by setting -DENABLE_SPHINX_DOCS=OFF.")
endif()

find_package(LATEX COMPONENTS PDFLATEX BIBTEX)

if (ENABLE_SPHINX_API_DOCS)
  set(SPHINX_OPTS "-t" "use_breathe")
else()
  set(SPHINX_OPTS " ")
endif()

add_custom_target(sphinx-html
  COMMAND ${SPHINX_EXE} -M html
  "${CMAKE_CURRENT_SOURCE_DIR}/source" "${CMAKE_CURRENT_BINARY_DIR}"
  -c ${CMAKE_CURRENT_SOURCE_DIR}/source ${SPHINX_OPTS})
add_dependencies(sphinx-html api-docs)

if (LATEX_FOUND)
  add_custom_target(sphinx-pdf
    COMMAND ${SPHINX_EXE} -M latexpdf
    "${CMAKE_CURRENT_SOURCE_DIR}/source" "${CMAKE_CURRENT_BINARY_DIR}"
    -c ${CMAKE_CURRENT_SOURCE_DIR}/source ${SPHINX_OPTS})
  add_dependencies(sphinx-pdf api-docs)
else()
  message(WARNING "Disabling PDF generation because of missing pdflatex")
endif()

add_dependencies(docs sphinx-html)
