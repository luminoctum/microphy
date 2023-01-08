# Installs microphy library
#
#   library: microphy_debug
#   library: microphy_release
#
# Define the following variables
#
#   MICROPHY_INCLUDE_DIR
#   MICROPHY_LIBRARY_DEBUG
#   MICROPHY_LIBRARY_RELEASE
#   MICROPHY_FOUND

file(GLOB _src_files 
  *.cpp
  #reaction/*.cpp
  )

foreach(build ${BUILD_TYPES})
  string(TOLOWER ${build} buildl)
  add_library(microphy_${buildl} OBJECT ${_src_files})

  set_target_properties(microphy_${buildl}
    PROPERTIES
    COMPILE_FLAGS ${CMAKE_CXX_FLAGS_${build}}
    )

  target_include_directories(microphy_${buildl}
    SYSTEM PRIVATE
    ${ATHENAPP_INCLUDE_DIR}
    ${CLIMATH_INCLUDE_DIR}
    ${DEBUGGER_INCLUDE_DIR}
    ${CLIUTILS_INCLUDE_DIR}
    ${MOLECULE_INCLUDE_DIR}
    ${CANTERA_INCLUDE_DIR}
    ${MPI_CXX_INCLUDE_PATH}
    ${CMAKE_BINARY_DIR}
    )

  target_link_directories(microphy_${buildl}
    PRIVATE
    ${CANTERA_LIBRARY_DIR}
    )

  target_link_libraries(microphy_${buildl}
    ${ATHENAPP_LIBRARY_${build}}
    ${CLIMATH_LIBRARY_${build}}
    ${DEBUGGER_LIBRARY_${build}}
    ${CLIUTILS_LIBRARY_${build}}
    ${MOLECULE_LIBRARY_${build}}
    ${CANTERA_LIBRARY}
    ${MPI_CXX_LIBRARIES}
    )

  set(MICROPHY_LIBRARY_${build} microphy_${buildl}
    CACHE STRING 
    "microphy library ${build}")
  mark_as_advanced(MICROPHY_LIBRARY_${build})
endforeach()

set(MICROPHY_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  CACHE PATH 
  "microphy include directory")
set(MICROPHY_FOUND TRUE CACHE BOOL "FIND MICROPHY LIBRARY")
mark_as_advanced(MICROPHY_INCLUDE_DIR)
mark_as_advanced(MICROPHY_FOUND)
