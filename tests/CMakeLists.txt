INCLUDE(CTest)

macro (create_test targetFile)
    add_test(NAME ${targetFile} COMMAND $<TARGET_FILE:${targetFile}>)
    set_tests_properties(${targetFile} PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
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

ADD_EXECUTABLE(ooMeshTest.exe OOMeshTest.f90)
TARGET_LINK_LIBRARIES(ooMeshTest.exe lpmFortran)

ADD_EXECUTABLE(ooParticlesAndFieldTest.exe OOParticlesAndFieldTest.f90)
TARGET_LINK_LIBRARIES(ooParticlesAndFieldTest.exe lpmFortran)

ADD_EXECUTABLE(ooPlaneMeshTest.exe OOPlaneMeshTest.f90)
TARGET_LINK_LIBRARIES(ooPlaneMeshTest.exe lpmFortran)

ADD_EXECUTABLE(ooIcosTriMeshTest.exe OOIcosTriTests.f90)
TARGET_LINK_LIBRARIES(ooIcosTriMeshTest.exe lpmFortran)

ADD_EXECUTABLE(ooCubedSphereTest.exe OOCubedSphereTests.f90)
TARGET_LINK_LIBRARIES(ooCubedSphereTest.exe lpmFortran)

ADD_EXECUTABLE(cubicGllTest.exe CubicGLLTest.f90)
TARGET_LINK_LIBRARIES(cubicGllTest.exe lpmFortran)

# Unit Tests
create_test(particlesAndFieldTest.exe)
create_test(stdIntVectorTest.exe)
create_test(planarMeshTest.exe 3)
create_test(icosTriSphereTest.exe 2)
#create_test(cubedSphereTest.exe 4)
create_test(ooMeshTest.exe 0)
create_test(ooParticlesAndFieldTest.exe)
create_test(ooPlaneMeshTest.exe)
create_test(ooIcosTriMeshTest.exe 3)
create_test(ooCubedSphereTest.exe 3)
create_test(cubicGllTest.exe)

# Regression tests
create_test(BivarConvergenceTest.exe)
add_test(NAME ssrfpackTest COMMAND mpirun "-np" 4 $<TARGET_FILE:ssrfpackTest.exe> 3 4)
set_tests_properties(ssrfpackTest PROPERTIES PASS_REGULAR_EXPRESSION "ssrfTest_00 : interp err =   0.16686400E-03")
add_test(NAME spherePSETest COMMAND mpirun "-np" 4 $<TARGET_FILE:spherePSETest.exe> ${PROJECT_BINARY_DIR}/../tests/spherePSEDefault.namelist)
set_tests_properties(spherePSETest PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
