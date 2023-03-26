CMAKE_MINIMUM_REQUIRED(VERSION 3.12)

project (LPM Fortran)
set ( LPM_VERSION_MAJOR 2)
set ( LPM_VERSION_MINOR 1)

FIND_PACKAGE(MPI REQUIRED)
FIND_PACKAGE(LAPACK REQUIRED)

enable_language(Fortran)
set(CMAKE_Fortran_FORMAT FREE)

# set default build type to release
if ( NOT CMAKE_BUILD_TYPE)
	set (CMAKE_BUILD_TYPE RELEASE CACHE STRING "Choose type of build. Options are : NONE, DEBUG, RELEASE." FORCE)
endif ( NOT CMAKE_BUILD_TYPE)

# default installation
#get_filename_component (default_prefix "." ABSOLUTE)
#set (CMAKE_INSTALL_PREFIX ${default_prefix}/install CACHE STRING
#      "Choose the installation directory; by default it installs in the install directory."
#      FORCE)

# FFLAGS depend on the compiler
if ( ${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
  # gfortran
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -fopenmp -ffree-form -cpp -fallow-argument-mismatch")
  set (CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-O2 -cpp -fallow-argument-mismatch")
  set (CMAKE_Fortran_FLAGS_MINSIZEREL "-O2 -cpp -fallow-argument-mismatch")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g3 -cpp -fopenmp -fbacktrace -fbounds-check -ffree-form -fallow-argument-mismatch")
elseif ( ${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
  # ifort
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -qopenmp -free -fpp")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -qopenmp -free -fpp -debug all -debug-parameters -warn all -traceback")
else ( ${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
  message ("No optimized Fortran compiler flags are known, we just try -O2...")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif ( ${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")



configure_file(cubedSphereSeed.dat ${CMAKE_BINARY_DIR}/cubedSphereSeed.dat COPYONLY)
configure_file(icosTriSeed.dat ${CMAKE_BINARY_DIR}/icosTriSeed.dat COPYONLY)
configure_file(quadRectSeed.dat ${CMAKE_BINARY_DIR}/quadRectSeed.dat COPYONLY)
configure_file(triHexSeed.dat ${CMAKE_BINARY_DIR}/triHexSeed.dat COPYONLY)
configure_file(betaPlaneSeed.dat ${CMAKE_BINARY_DIR}/betaPlaneSeed.dat COPYONLY)

configure_file(cubedSphereSeed.dat ${CMAKE_BINARY_DIR}/tests/cubedSphereSeed.dat COPYONLY)
configure_file(icosTriSeed.dat ${CMAKE_BINARY_DIR}/tests/icosTriSeed.dat COPYONLY)
configure_file(quadRectSeed.dat ${CMAKE_BINARY_DIR}/tests/quadRectSeed.dat COPYONLY)
configure_file(triHexSeed.dat ${CMAKE_BINARY_DIR}/tests/triHexSeed.dat COPYONLY)
configure_file(betaPlaneSeed.dat ${CMAKE_BINARY_DIR}/tests/betaPlaneSeed.dat COPYONLY)

configure_file(cubedSphereSeed.dat ${CMAKE_BINARY_DIR}/examples/cubedSphereSeed.dat COPYONLY)
configure_file(icosTriSeed.dat ${CMAKE_BINARY_DIR}/examples/icosTriSeed.dat COPYONLY)
configure_file(quadRectSeed.dat ${CMAKE_BINARY_DIR}/examples/quadRectSeed.dat COPYONLY)
configure_file(triHexSeed.dat ${CMAKE_BINARY_DIR}/examples/triHexSeed.dat COPYONLY)
configure_file(betaPlaneSeed.dat ${CMAKE_BINARY_DIR}/examples/betaPlaneSeed.dat COPYONLY)

ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(examples)
ENABLE_TESTING()
ADD_SUBDIRECTORY(tests)

INCLUDE_DIRECTORIES(${PROJECT_BINARY_DIR} ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR}/src)

configure_file("${CMAKE_SOURCE_DIR}/LpmConfig.h.in" "${PROJECT_BINARY_DIR}/LpmConfig.h")



