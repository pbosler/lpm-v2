program IcosTriMeshTester

use NumberKindsModule
use OutputWriterModule
use LoggerModule
use PolyMesh2dModule
use ParticlesModule
use EdgesModule
use FacesModule
use STDIntVectorModule
use SphereGeomModule

implicit none

integer(kint) :: i, j

type(PolyMesh2d) :: sphere
integer(kint) :: initNest
integer(kint) :: maxNest
integer(kint) :: amrLimit
real(kreal) :: radius
real(kreal) ::  testLength, minEdgeLength0, avgLength
integer(kint) :: nEdges

type(Logger) :: exeLog


real(kreal) :: lons(180)
real(kreal) :: lats(91)
integer(kint) :: triFaces(91,180), nearParticles(91,180)

integer(kint) :: randomVertex
integer(kint) :: randomFace
real(kreal) :: randReal
logical(klog) :: keepGoing
type(STDIntVector) :: edgesAroundFace, adjacentFaces, verticesAroundFace, facesAroundVertex
real(kreal) :: phys(3), xyz(3)

integer(kint) :: narg
character(len=100) :: arg, filename
character(MAX_STRING_LENGTH) :: logstring

call New(exeLog, DEBUG_LOGGING_LEVEL)

call LogMessage(exelog, TRACE_LOGGING_LEVEL, "IcosTriMeshTest : ", "1. build an icosahedral triangle sphere mesh.")
call LogMessage(exelog, TRACE_LOGGING_LEVEL, "IcosTriMeshTest : ", "2. test point query algorithm.")
call LogMessage(exeLog, TRACE_LOGGING_LEVEL, "IcosTriMeshTest : ", "3. test adjacency algorithms.")
call LogMessage(exelog, TRACE_LOGGING_LEVEL, "IcosTriMeshTest : ", "4. output to paraview/vtk.")

narg = IARGC()
if ( narg == 0 ) then
	initNest = 0 
else
	call GETARG(1, arg)
	read(arg, *) initNest
endif

maxNest = initNest
amrLimit = 0
radius = 1.0_kreal

do i = 1, 180
	lons(i) = 2.0_kreal * PI / 180.0_kreal  * (i - 1)
enddo
do i = 1, 91
	lats(i) = -PI/2.0_kreal + 2.0_kreal * PI / 180.0_kreal * (i-1)
enddo

call LogMessage(exeLog, TRACE_LOGGING_LEVEL, "Test 1 : ", " constructing icosahedral triangle sphere mesh.")

	call New(sphere, ICOS_TRI_SPHERE_SEED, initNest, maxNest, amrLimit, radius )
	call LogStats(sphere, exeLog )
	do i = 1, sphere%particles%N
		if ( sphere%particles%isActive(i) .AND. sphere%particles%isPassive(i) ) then
			call LogMessage(exeLog, WARNING_LOGGING_LEVEL,"sphere Particles WARNING : both active and passive = .TRUE. at particle ", i)
			testPass = .FALSE.
		endif
	enddo
call LogMessage(exeLog, TRACE_LOGGING_LEVEL, "Test 1 : ", " complete.")

avgLength = 0.0_kreal
nEdges = 30
minEdgeLength0 = 1.0e20

if ( sphere%edges%N - count(sphere%edges%hasChildren) /= nEdges) then
	call LogMessage(exeLog, ERROR_LOGGING_LEVEL, "ERROR : nEdges = ", nEdges)
	call LogMessage(exeLog, ERROR_LOGGING_LEVEL, "edges%N_Active =  ", sphere%edges%N - count(sphere%edges%hasChildren))
	testPass = .FALSE.
endif

do i = 1, sphere%edges%N
	if ( .NOT. sphere%edges%hasChildren(i) ) then
		testLength = EdgeLength( sphere%edges, i, sphere%particles)
		avgLength = avgLength + EdgeLength( sphere%edges, i, sphere%particles)
		if ( testLength < minEdgeLength0 )  minEdgeLength0 = testLength
	endif
enddo

avgLength = avgLength / real(nEdges,kreal)

write(logstring,'(A,I3,A,E24.7)') " init nest = ", initNest, ", maxEdgeLength = ", maxEdgeLength(sphere%edges, sphere%particles)
call LogMessage(exeLog,TRACE_LOGGING_LEVEL,"MAX_EDGE_LENGTH: ", trim(logString))
write(logstring,'(A,I3,A,E24.7)') " init nest = ", initNest, ", minEdgeLength0 = ", minEdgeLength0
call LogMessage(exeLog,TRACE_LOGGING_LEVEL,"MIN_EDGE_LENGTH: ", trim(logString))
write(logstring,'(A,I3,A,E24.7)') " init nest = ", initNest, ", avgEdgeLength = ", avgLength
call LogMessage(exeLog,TRACE_LOGGING_LEVEL,"AVG_EDGE_LENGTH: ", trim(logString))




call LogMessage(exeLog, TRACE_LOGGING_LEVEL, "Test 2 : ", " correlating a uniform lat-lon grid to the icos tri mesh.")
	do j = 1, 180
		do i = 1, 91
			xyz = [ radius * cos(lats(i)) * cos(lons(j)), radius * cos(lats(i)) * sin(lons(j)), radius * sin(lats(i)) ]
			triFaces(i,j) = LocateFaceContainingPoint( sphere, xyz )
			nearParticles(i,j) = nearestParticle(sphere, xyz)
		enddo
	enddo

	write(filename,'(A,I1,A)') "icosTriSphere_", initNest, ".m"
	open(unit=WRITE_UNIT_1, file = filename, status='REPLACE',action='WRITE')

		write(WRITE_UNIT_1,'(A)',advance='NO') "lons = ["
		do i = 1, 179
			write(WRITE_UNIT_1,'(F12.6,A)',advance='NO') lons(i), ", "
		enddo
		write(WRITE_UNIT_1,'(F12.6,A)') lons(180), "];"

		write(WRITE_UNIT_1,'(A)',advance='NO') "lats = ["
		do i = 1, 90
			write(WRITE_UNIT_1,'(F12.6,A)',advance='NO') lats(i), ", "
		enddo
		write(WRITE_UNIT_1,'(F12.6,A)') lats(91), "];"

		write(WRITE_UNIT_1,'(A)',advance='NO') "triFaces = ["
		do i = 1, 90
			do j = 1, 179
				write(WRITE_UNIT_1,'(I8,A)',advance='NO') triFaces(i,j), ", "
			enddo
			write(WRITE_UNIT_1,'(I8,A)') triFaces(i,180), "; "
		enddo
		do j = 1, 179
			write(WRITE_UNIT_1,'(I8,A)',advance='NO') triFaces(91,j), ", "
		enddo
		write(WRITE_UNIT_1,'(I8,A)') triFaces(91,180), "]; "
		
		call WriteToMatlab(nearParticles, WRITE_UNIT_1, "nearestParticle")

		write(WRITE_UNIT_1,'(A,I8,A)') "figure(1); clf; hold on; contour(lons, lats, triFaces, ",&
				 sphere%faces%N_Active, "); colorbar;"
call LogMessage(exeLog, TRACE_LOGGING_LEVEL, "Test 2 : ", " complete.")

call LogMessage(exeLog, TRACE_LOGGING_LEVEL, "Test 3 : ", " testing adjacency subroutines.")

		call RANDOM_SEED()
		keepGoing = .TRUE.
		do while (keepGoing)
			call RANDOM_NUMBER(randReal)
			randomFace = 1 + floor( sphere%faces%N * randReal)
			if ( .NOT. sphere%faces%hasChildren(randomFace) ) keepGoing = .FALSE.
		enddo

		keepGoing = .TRUE.
		do while (keepGoing)
			call RANDOM_NUMBER(randReal)
			randomVertex = 1 + floor( sphere%particles%N * randReal)
			if ( sphere%particles%isPassive(randomVertex) .AND. .NOT. sphere%particles%isActive(randomVertex) ) keepGoing = .FALSE.
		enddo

		call CCWEdgesAroundFace( sphere, edgesAroundFace, randomFace)
		call CCWVerticesAroundFace( sphere, verticesAroundFace, randomFace)
		call CCWAdjacentFaces( sphere, adjacentFaces, randomFace )
		
		call CCWFacesAroundVertex( sphere, facesAroundVertex, randomVertex)
		
		phys = PhysCoord(sphere%particles, randomVertex)
		write(WRITE_UNIT_1,'(A,2(F12.6,A))') "randVert = [ ", Longitude(phys)," , ", Latitude(phys), "];"
		phys = FaceCenterPhysCoord( sphere%faces, randomFace, sphere%particles)
		write(WRITE_UNIT_1,'(A,2(F12.6,A))') "randFace = [ ", Longitude(phys)," , ", Latitude(phys), "];"
		
		write(WRITE_UNIT_1,'(A)', advance='NO') "facesNearVert = ["
		do i = 1, facesAroundVertex%N - 1
			phys = FaceCenterPhysCoord(sphere%faces, facesAroundVertex%int(i), sphere%particles )
			write(WRITE_UNIT_1,'(2(F12.6,A))') Longitude(phys), ", ", Latitude(phys), "; ..."
		enddo
		phys = FaceCenterPhysCoord(sphere%faces, facesAroundVertex%int(facesAroundVertex%N), sphere%particles )
		write(WRITE_UNIT_1,'(2(F12.6,A))') Longitude(phys), ", ", Latitude(phys), "]; "
		
		write(WRITE_UNIT_1,'(A)') "edgesAroundFaceLon = ["
		do i = 1, edgesAroundFace%N
			phys = PhysCoord(sphere%particles, sphere%edges%orig(edgesAroundFace%int(i)) )
			xyz = PhysCoord(sphere%particles, sphere%edges%dest(edgesAroundFace%int(i)) )
			write(WRITE_UNIT_1,'(2(F12.6,A))') Longitude(phys), ", ", Longitude(xyz), "; ..."
		enddo
		phys = PhysCoord(sphere%particles, sphere%edges%orig(edgesAroundFace%int(1)) )
		xyz = PhysCoord(sphere%particles, sphere%edges%dest(edgesAroundFace%int(1)) )
		write(WRITE_UNIT_1,'(2(F12.6,A))') Longitude(phys), ", ", Longitude(xyz), "];"
		
		write(WRITE_UNIT_1,'(A)') "edgesAroundFaceLat = ["
		do i = 1, edgesAroundFace%N 
			phys = PhysCoord(sphere%particles, sphere%edges%orig(edgesAroundFace%int(i)) )
			xyz = PhysCoord(sphere%particles, sphere%edges%dest(edgesAroundFace%int(i)) )
			write(WRITE_UNIT_1,'(2(F12.6,A))') Latitude(phys), ", ", Latitude(xyz), "; ..."
		enddo
		phys = PhysCoord(sphere%particles, sphere%edges%orig(edgesAroundFace%int(1)) )
		xyz = PhysCoord(sphere%particles, sphere%edges%dest(edgesAroundFace%int(1)) )
		write(WRITE_UNIT_1,'(2(F12.6,A))') Latitude(phys), ", ", Latitude(xyz), "];"
		
		write(WRITE_UNIT_1,'(A)',advance='NO') "vertsAroundFace = ["
		do i = 1, verticesAroundFace%N - 1
			phys = PhysCoord(sphere%particles, verticesAroundFace%int(i))
			write(WRITE_UNIT_1,'(2(F12.6,A))') Longitude(phys), ", ", Latitude(phys), "; ..."
		enddo
		phys = PhysCoord(sphere%particles, verticesAroundFace%int(verticesAroundFace%N))
		write(WRITE_UNIT_1,'(2(F12.6,A))') Longitude(phys), ", ", Latitude(phys), "];"
		
		write(WRITE_UNIT_1,'(A)',advance='NO') "adjFaces = ["
		do i = 1, adjacentFaces%N - 1
			phys = FaceCenterPhysCoord(sphere%faces, adjacentFaces%int(i), sphere%particles)
			write(WRITE_UNIT_1,'(2(F12.6,A))') Longitude(phys), ", ", Latitude(phys), "; ..."
		enddo
		phys = FaceCenterPhysCoord(sphere%faces, adjacentFaces%int(adjacentFaces%N), sphere%particles)
		write(WRITE_UNIT_1,'(2(F12.6,A))') Longitude(phys), ", ", Latitude(phys), "]; "
		
		write(WRITE_UNIT_1,'(A)') "plot(randFace(1), randFace(2),'ksq','MarkerSize',14)"
		write(WRITE_UNIT_1,'(A)') "plot(edgesAroundFaceLon,edgesAroundFaceLat,'k-','LineWidth',2)"
		write(WRITE_UNIT_1,'(A)') "plot(vertsAroundFace(:,1), vertsAroundFace(:,2),'ko','MarkerSize',14)"
		write(WRITE_UNIT_1,'(A)') "plot(randVert(1), randVert(2),'k*','MarkerSize',14)"
		write(WRITE_UNIT_1,'(A)') "plot(facesNearVert(:,1), facesNearVert(:,2),'ksq','MarkerSize',14)"
		
	close(WRITE_UNIT_1)
call LogMessage(exeLog, TRACE_LOGGING_LEVEL, "Test 3 : ", " complete.")
	

call LogMessage(exeLog, TRACE_LOGGING_LEVEL, "Test 4 : ", " writing .vtk output.")

write(filename,'(A,I1,A)') "icosTriSphere_", initNest, ".vtk"
open(unit=WRITE_UNIT_1, file=filename, status='REPLACE',action='WRITE')
call WriteMeshToVTKPolyData( sphere, WRITE_UNIT_1, 'icosTriTest')

call LogMessage(exeLog, TRACE_LOGGING_LEVEL, "Test 4 : ", " complete.")
close(WRITE_UNIT_1)

call Delete( sphere )

if (testPass) call LogMessage(exeLog, TRACE_LOGGING_LEVEL, "Test result: " , "PASSED.")
call Delete(exeLog)

end program
