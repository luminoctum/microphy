# Installs condense library
#
#   library: condense_debug
#   library: condense_release
#
# Define the following variables
#
#   CONDENSE_INCLUDE_DIR
#   CONDENSE_LIBRARY_DEBUG
#   CONDENSE_LIBRARY_RELEASE
#   CONDENSE_FOUND

file(GLOB _src_files 
  gasphase/*.cpp
  )

foreach(build ${BUILD_TYPES})
  string(TOLOWER ${build} buildl)
  add_library(condense_${buildl} ${_src_files})

  set_target_properties(condense_${buildl}
    PROPERTIES
    COMPILE_FLAGS ${CMAKE_CXX_FLAGS_${build}}
    )

  target_include_directories(condense_${buildl}
    SYSTEM PRIVATE
    ${ATHENAPP_INCLUDE_DIR}
    ${DEBUGGER_INCLUDE_DIR}
    ${MOLECULE_INCLUDE_DIR}
    ${CMAKE_BINARY_DIR}
    )

  target_link_libraries(condense_${buildl}
    ${ATHENAPP_LIBRARY_${build}}
    ${DEBUGGER_LIBRARY_${build}}
    ${MOLECULE_LIBRARY_${build}}
    )

  set(CONDENSE_LIBRARY_${build} condense_${buildl}
    CACHE STRING 
    "gasphase library ${build}")
  mark_as_advanced(CONDENSE_LIBRARY_${build})
endforeach()

set(CONDENSE_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  CACHE PATH 
  "gasphase include directory")
set(CONDENSE_FOUND TRUE CACHE BOOL "FIND CONDENSE LIBRARY")
mark_as_advanced(CONDENSE_INCLUDE_DIR)
mark_as_advanced(CONDENSE_FOUND)
