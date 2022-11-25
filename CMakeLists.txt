## 0. set up project ##
cmake_minimum_required(VERSION 3.18)
message(STATUS "Set up canoe project ...")

project(
  GasPhase
  VERSION 0.0.1
  LANGUAGES CXX
  )
set(CMAKE_CXX_STANDARD 14)

LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules/)

## 1. set up compiler flags ##
message(STATUS "")
message(STATUS "Set up project compiler flags ...")
message(STATUS "Include ${CMAKE_SOURCE_DIR}/cmake/setup_compiler_flags.cmake")
include(${CMAKE_SOURCE_DIR}/cmake/setup_compiler_flags.cmake)

## 2. set up project configuration ##
message(STATUS "")
message(STATUS "Set up project configuration ...")
message(STATUS "Include ${CMAKE_SOURCE_DIR}/cmake/setup_configure.cmake")
include(${CMAKE_SOURCE_DIR}/cmake/setup_configure.cmake)

## 3. set up project paths ##
message(STATUS "")
message(STATUS "Set up project paths ...")
set(EXECUTABLE_RELDIR "bin")
set(INCLUDE_RELDIR "include")
set(SOURCE_RELDIR "source")

## 4. set up project system libraries ##
message(STATUS "")
message(STATUS "Set up project system libraries ...")

if (${UseMPI})
  find_package(MPI REQUIRED)
  set(EnableMPI MPI_PARALLEL)
else()
  set(EnableMPI NOT_MPI_PARALLEL)
endif()

## 5. set up project external libraries ##
message(STATUS "")
message(STATUS "Set up project external libraries ...")

# compile climath
if (${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_SOURCE_DIR})
  add_subdirectory(climath)
endif()

# compile debugger
if (${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_SOURCE_DIR})
  add_subdirectory(debugger)
endif()

# compile athenapp
if (${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_SOURCE_DIR})
  add_subdirectory(athenapp)
endif()

# compile molecule
if (${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_SOURCE_DIR})
  add_subdirectory(molecule)
endif()

## 6. set up compilation and installation ##
message(STATUS "")
message(STATUS "Set up project complication and installation ...")

# configure file
configure_file(${CMAKE_SOURCE_DIR}/configure.hpp.in configure.hpp @ONLY)

add_subdirectory(src)

if (${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_SOURCE_DIR})
  add_subdirectory(tests)
endif()
