INCLUDE(CTest)

macro (create_test targetFile)
    add_test(NAME ${targetFile} COMMAND $<TARGET_FILE:${targetFile}>)
endmacro (create_test)

INCLUDE_DIRECTORIES(${PROJECT_BINARY_DIR} ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR}/src)

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
create_test(particlesAndFieldTest.exe)
create_test(stdIntVectorTest.exe)
create_test(planarMeshTest.exe 3)
create_test(icosTriSphereTest.exe 2)
create_test(cubedSphereTest.exe 4)

# Regression tests
#create_test(BivarConvergenceTest.exe)
add_test(NAME ssrfpacTest COMMAND mpirun "-np" 4 $<TARGET_FILE:ssrfpackTest.exe> 3 4)
add_test(NAME spherePSETest COMMAND mpirun "-np" 4 $<TARGET_FILE:spherePSETest.exe> ${PROJECT_BINARY_DIR}/../tests/spherePSEDefault.namelist)

