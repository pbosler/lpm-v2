INCLUDE_DIRECTORIES(${PROJECT_BINARY_DIR} ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR}/src)

# APPLICATION DRIVERS
ADD_EXECUTABLE( advectGaussHillsSphere.exe AdvectGaussianHills.f90)
TARGET_LINK_LIBRARIES( advectGaussHillsSphere.exe lpmFortran)

ADD_EXECUTABLE( advectGHillsDivWind.exe AdvectDivWindGaussHills.f90 )
TARGET_LINK_LIBRARIES(advectGHillsDivWind.exe lpmFortran )

ADD_EXECUTABLE( advectRossbyWave.exe AdvectRHWave.f90 )
TARGET_LINK_LIBRARIES( advectRossbyWave.exe lpmFortran)

ADD_EXECUTABLE( advectDivergentWind.exe AdvectDivergentWind.f90)
TARGET_LINK_LIBRARIES( advectDivergentWind.exe lpmFortran)

ADD_EXECUTABLE( advectSlotC.exe AdvectSlottedC.f90)
TARGET_LINK_LIBRARIES( advectSlotC.exe lpmFortran)

ADD_EXECUTABLE( nitscheStrickland.exe NitscheStrickland.f90 )
TARGET_LINK_LIBRARIES( nitscheStrickland.exe lpmFortran)

ADD_EXECUTABLE( parabolicFlood.exe ParabolicFlood.f90)
TARGET_LINK_LIBRARIES(parabolicFlood.exe lpmFortran)

ADD_EXECUTABLE( sphereSWETC2.exe SWETestCase2.f90 )
TARGET_LINK_LIBRARIES( sphereSWETC2.exe lpmFortran )

ADD_EXECUTABLE( bveSolidBodyRotation.exe BVESolidBody.f90 )
TARGET_LINK_LIBRARIES( bveSolidBodyRotation.exe lpmFortran )

ADD_EXECUTABLE( bveGaussianVortex.exe BVESingleGaussianVortex.f90)
TARGET_LINK_LIBRARIES( bveGaussianVortex.exe lpmFortran)

ADD_EXECUTABLE( betaPlaneGaussianVortex.exe BetaPlaneGaussianVortex.f90)
TARGET_LINK_LIBRARIES( betaPlaneGaussianVortex.exe lpmFortran)

ADD_EXECUTABLE( collidingDipoles.exe CollidingDipoles.f90 )
TARGET_LINK_LIBRARIES( collidingDipoles.exe lpmFortran )

ADD_EXECUTABLE( ellipticVortex.exe EllipticVortex.f90 )
TARGET_LINK_LIBRARIES( ellipticVortex.exe lpmFortran )

ADD_EXECUTABLE( sweSteadyStatePlane.exe SWESteadyState.f90)
TARGET_LINK_LIBRARIES( sweSteadyStatePlane.exe lpmFortran )

ADD_EXECUTABLE( planarGravityWave.exe PlaneSWEGravityWaves.f90)
TARGET_LINK_LIBRARIES( planarGravityWave.exe lpmFortran)

ADD_EXECUTABLE( bveRossbyHaurwitzWave54.exe RossbyHaurwitz54.f90 )
TARGET_LINK_LIBRARIES( bveRossbyHaurwitzWave54.exe lpmFortran)

ADD_EXECUTABLE( transportDemo.exe TransportDemo.f90)
TARGET_LINK_LIBRARIES( transportDemo.exe lpmFortran)
