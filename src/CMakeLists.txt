add_library( lpmFortran TypeDefs.f90 Utilities.f90 OutputWriter.f90 Logger.f90 STDIntVector.f90
						SphereGeometry.f90 PlaneGeometry.f90
						Particles.f90 Field.f90
						Edges.f90 Faces.f90 PolyMesh2d.f90 #CubicEdges.f90 CubicQuadFaces.f90
						PlanarIncompressible.f90 PlaneIncompressibleSolver.f90
						bivar.f90 BIVARParticlesInterface.f90 BIVARRemesh.f90
						PSEDirectSum.f90
						stripack.f ssrfpack.f ssrfpackInterface.f90 SsrfpackRemesh.f90
						MPISetup.f90 PlanarSWE.f90 SWEPlaneSolver.f90
						BetaPlane.f90 BetaPlaneSolver.f90
						Refinement.f90 SphereSWE.f90 SphereSWESolver.f90
						SphereBVE.f90 SphereBVESolver.f90
						SphereTransport.f90 SphereTransportSolver.f90 SphereTransportVelocities.f90
						SphereTracers.f90
						ParticlesOO.f90 EdgesOO.f90 FacesOO.f90 PolyMeshOO.f90 FieldOO.f90 CubicGLL.f90
						MPIReplicatedData.f90 MeshedParticles.f90 RK4.f90
						)
						
if ( ${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
    SET_SOURCE_FILES_PROPERTIES(stripack.f PROPERTIES COMPILE_FLAGS -ffixed-form)
    SET_SOURCE_FILES_PROPERTIES(ssrfpack.f PROPERTIES COMPILE_FLAGS -ffixed-form)
elseif ( ${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
    SET_SOURCE_FILES_PROPERTIES(stripack.f PROPERTIES COMPILE_FLAGS -nofree)
    SET_SOURCE_FILES_PROPERTIES(ssrfpack.f PROPERTIES COMPILE_FLAGS -nofree)	
endif ( ${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
						
target_link_libraries( lpmFortran ${LAPACK_LIBRARIES})

INCLUDE_DIRECTORIES(${PROJECT_BINARY_DIR} ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR}/src)
