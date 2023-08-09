cmake_minimum_required(VERSION 3.13)

project(fort-dist LANGUAGES Fortran)

find_package(MPI REQUIRED COMPONENTS Fortran)
message(STATUS "${MPI_Fortran_LIBRARIES}")
message(STATUS ${MPI_Fortran_INCLUDE_PATH})
message(STATUS ${MPI_INCLUDE_PATH})
message(STATUS "test")

find_package(gmxfort CONFIG)

set(CMAKE_Fortran_FLAGS "-O3")

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

include(GNUInstallDirs)
include_directories(/usr/local/include)
link_directories(/usr/local/lib)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY
  ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY
  ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
  ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})

# defines targets and sources
add_subdirectory(src)

# contains an "external" library we will link to
#add_subdirectory(external)

# enable testing and define tests
#enable_testing()
#add_subdirectory(tests)