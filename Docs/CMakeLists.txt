cmake_minimum_required (VERSION 3.10)
project(PeleLM-Docs NONE)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

find_package(Sphinx REQUIRED)
#find_package(Doxygen)

add_custom_target(sphinx ALL
  COMMAND ${SPHINX_EXECUTABLE} -M html
  "${CMAKE_CURRENT_SOURCE_DIR}/source" "${CMAKE_CURRENT_BINARY_DIR}/source"
  -c ${CMAKE_CURRENT_SOURCE_DIR}/source)

#if(DOXYGEN_FOUND)
#  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doxygen/Doxyfile.in Doxyfile @ONLY)
#  add_custom_target(doxygen ALL
#      COMMAND ${DOXYGEN_EXECUTABLE} "${CMAKE_CURRENT_BINARY_DIR}/Doxyfile" )
#  add_dependencies(doxygen sphinx)
#endif()
