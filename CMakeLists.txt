cmake_minimum_required (VERSION 2.6)
project (bmsap Fortran)

# default build type: Release
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

# default install prefix: ~/.local
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "$ENV{HOME}/.local"
      CACHE PATH "Install path prefix, prepended onto install directories."
      FORCE)
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Install prefix: ${CMAKE_INSTALL_PREFIX}")

# set default architecture flags to optimize fot the furrent system,
# but don't overide user settings
if(NOT CMAKE_Fortran_FLAGS MATCHES "-march|-mtune")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -march=native")
endif()

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -fno-align-commons")

add_subdirectory (src)

# install dependencies
#install (FILES sample.in sample_multi.in DESTINATION ${PROJECT_NAME})
#install (DIRECTORY spectra DESTINATION ${PROJECT_NAME})
install (FILES sample.in sample_multi.in DESTINATION results)
install (DIRECTORY spectra DESTINATION share)

