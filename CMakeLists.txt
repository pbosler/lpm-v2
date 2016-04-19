CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

project (LPM)
set ( LPM_VERSION_MAJOR 1)
set ( LPM_VERSION_MINOR 0)

FIND_PACKAGE(MPI REQUIRED)
FIND_PACKAGE(LAPACK REQUIRED)

enable_language(Fortran)

# set default build type to release
if ( NOT CMAKE_BUILD_TYPE)
	set (CMAKE_BUILD_TYPE RELEASE CACHE STRING "Choose type of build. Options are : NONE, DEBUG, RELEASE." FORCE)
endif ( NOT CMAKE_BUILD_TYPE)

# default installation
get_filename_component (default_prefix "." ABSOLUTE)
set (CMAKE_INSTALL_PREFIX ${default_prefix}/install CACHE STRING
      "Choose the installation directory; by default it installs in the install directory."
      FORCE)
      
# FFLAGS depend on the compiler
if ( ${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
  # gfortran
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2 -fopenmp")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g3 -fopenmp -fbacktrace -fbounds-check")
elseif ( ${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
  # ifort
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2 -openmp")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -openmp -debug all -debug-parameters -warn all -traceback")
else ( ${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
  message ("No optimized Fortran compiler flags are known, we just try -O2...")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif ( ${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")

add_library( lpmFortran NumberKinds.f90 OutputWriter.f90 Logger.f90 STDIntVector.f90
						SphereGeometry.f90 PlaneGeometry.f90 
						Particles.f90 Field.f90
						Edges.f90 Faces.f90 PolyMesh2d.f90 
						PlanarIncompressible.f90 PlaneIncompressibleSolver.f90
						bivar.f90 BIVARParticlesInterface.f90 BIVARRemesh.f90
						PSEDirectSum.f90 
						stripack.f ssrfpack.f ssrfpackInterface.f90 SsrfpackRemesh.f90
						MPISetup.f90 PlanarSWE.f90 SWEPlaneSolver.f90
						BetaPlane.f90 BetaPlaneSolver.f90
						Refinement.f90 SphereSWE.f90 SphereSWESolver.f90
						SphereBVE.f90 SphereBVESolver.f90 
						SphereTransport.f90)
target_link_libraries( lpmFortran ${LAPACK_LIBRARIES})						

configure_file(cubedSphereSeed.dat ${CMAKE_BINARY_DIR}/cubedSphereSeed.dat COPYONLY)
configure_file(icosTriSeed.dat ${CMAKE_BINARY_DIR}/icosTriSeed.dat COPYONLY)
configure_file(quadRectSeed.dat ${CMAKE_BINARY_DIR}/quadRectSeed.dat COPYONLY)
configure_file(triHexSeed.dat ${CMAKE_BINARY_DIR}/triHexSeed.dat COPYONLY)
configure_file(betaPlaneSeed.dat ${CMAKE_BINARY_DIR}/betaPlaneSeed.dat COPYONLY)

# UNIT TEST EXECUTABLES
ADD_EXECUTABLE(particlesAndFieldTest.exe ParticlesAndFieldTest.f90 )
TARGET_LINK_LIBRARIES(particlesAndFieldTest.exe lpmFortran)

ADD_EXECUTABLE(planarMeshTest.exe PlanarMeshTest.f90)
TARGET_LINK_LIBRARIES(planarMeshTest.exe lpmFortran)

ADD_EXECUTABLE(stdIntVectorTest.exe STDIntVectorTest.f90)
TARGET_LINK_LIBRARIES(stdIntVectorTest.exe lpmFortran)

ADD_EXECUTABLE(icosTriSphereTest.exe IcosTriMeshTester.f90)
TARGET_LINK_LIBRARIES(icosTriSphereTest.exe lpmFortran)

ADD_EXECUTABLE(cubedSphereTest.exe CubedSphereTester.f90)
TARGET_LINK_LIBRARIES(cubedSphereTest.exe lpmFortran)

ADD_EXECUTABLE( BivarConvergenceTest.exe BivarConvergenceTests.f90) # Akima version of BIVAR
TARGET_LINK_LIBRARIES(BivarConvergenceTest.exe lpmFortran)

ADD_EXECUTABLE( ssrfpackTest.exe SSRFPACKTest.f90 )
TARGET_LINK_LIBRARIES( ssrfpackTest.exe lpmFortran )

ADD_EXECUTABLE( spherePSETest.exe SpherePSEConvTest.f90 )
TARGET_LINK_LIBRARIES( spherePSETest.exe lpmFortran )

# Unit Tests
ENABLE_TESTING()
ADD_TEST( particlesAndFieldTest ${PROJECT_BINARY_DIR}/particlesAndFieldTest.exe)
ADD_TEST( stdIntVectorTest ${PROJECT_BINARY_DIR}/stdIntVectorTest.exe)
ADD_TEST( planarMeshTest ${PROJECT_BINARY_DIR}/planarMeshTest.exe 3)
ADD_TEST( icosTriSphereTest ${PROJECT_BINARY_DIR}/icosTriSphereTest.exe 2 )
ADD_TEST( cubedSphereTest ${PROJECT_BINARY_DIR}/cubedSphereTest.exe 4 )


# APPLICATION DRIVERS
ADD_EXECUTABLE( nitscheStrickland.exe NitscheStrickland.f90 )
TARGET_LINK_LIBRARIES( nitscheStrickland.exe lpmFortran)

ADD_EXECUTABLE( parabolicFlood.exe ParabolicFlood.f90)
TARGET_LINK_LIBRARIES(parabolicFlood.exe lpmFortran)

ADD_EXECUTABLE( sphereSWETC2.exe SWETestCase2.f90 )
TARGET_LINK_LIBRARIES( sphereSWETC2.exe lpmFortran )

ADD_EXECUTABLE( bveSolidBodyRotation.exe BVESolidBody.f90 )
TARGET_LINK_LIBRARIES( bveSolidBodyRotation.exe lpmFortran )

ADD_EXECUTABLE( sphereGaussianVortex.exe BVESingleGaussianVortex.f90)
TARGET_LINK_LIBRARIES( sphereGaussianVortex.exe lpmFortran)

ADD_EXECUTABLE( betaPlaneGaussianVortex.exe BetaPlaneGaussianVortex.f90)
TARGET_LINK_LIBRARIES( betaPlaneGaussianVortex.exe lpmFortran)

ADD_EXECUTABLE( collidingDipoles.exe CollidingDipoles.f90 )
TARGET_LINK_LIBRARIES( collidingDipoles.exe lpmFortran )

ADD_EXECUTABLE( sweSteadyStatePlane.exe SWESteadyState.f90)
TARGET_LINK_LIBRARIES( sweSteadyStatePlane.exe lpmFortran )

ADD_EXECUTABLE( planarGravityWave.exe PlaneSWEGravityWaves.f90)
TARGET_LINK_LIBRARIES( planarGravityWave.exe lpmFortran)

ADD_EXECUTABLE( bveRossbyHaurwitzWave54.exe RossbyHaurwitz54.f90 )
TARGET_LINK_LIBRARIES( bveRossbyHaurwitzWave54.exe lpmFortran)

# INSTALL


