program cubedSphereMeshTest

use NumberKindsModule
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

type(Logger) :: exeLog

real(kreal) :: lons(180)
real(kreal) :: lats(91)
integer(kint) :: quadfaces(91,180)

integer(kint) :: randomVertex
integer(kint) :: randomFace
real(kreal) :: randReal
logical(klog) :: keepGoing
type(STDIntVector) :: edgesAroundFace, adjacentFaces, verticesAroundFace, facesAroundVertex
real(kreal) :: phys(3), xyz(3)

integer(kint) :: narg
character(len=100) :: arg, filename

!
!----------------
! PROGRAM START
!----------------
!

call New(exeLog, DEBUG_LOGGING_LEVEL)

call LogMessage(exelog, TRACE_LOGGING_LEVEL, "CubedSphereMeshTest : ", "1. build an cubed sphere mesh.")
call LogMessage(exelog, TRACE_LOGGING_LEVEL, "CubedSphereMeshTest : ", "2. test point query algorithm.")
call LogMessage(exeLog, TRACE_LOGGING_LEVEL, "CubedSphereMeshTest : ", "3. test adjacency algorithms.")
call LogMessage(exelog, TRACE_LOGGING_LEVEL, "CubedSphereMeshTest : ", "4. output to paraview/vtk.")


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

!
!----------------
! Test 1
!----------------
!
call LogMessage(exeLog, TRACE_LOGGING_LEVEL, "Test 1 : ", " constructing cubed sphere mesh.")

	call New(sphere, CUBED_SPHERE_SEED, initNest, maxNest, amrLimit, radius )
	call LogStats(sphere, exeLog )
	do i = 1, sphere%particles%N
		if ( sphere%particles%isActive(i) .AND. sphere%particles%isPassive(i) ) then
			call LogMessage(exeLog, WARNING_LOGGING_LEVEL,"sphere Particles WARNING : both active and passive = .TRUE. at particle ", i)
		endif
	enddo
call LogMessage(exeLog, TRACE_LOGGING_LEVEL, "Test 1 : ", " complete.")

!
!----------------
! Test 2
!----------------
!
call LogMessage(exeLog, TRACE_LOGGING_LEVEL, "Test 2 : ", " correlating a uniform lat-lon grid to the icos tri mesh.")
	do j = 1, 180
		do i = 1, 91
			xyz = [ radius * cos(lats(i)) * cos(lons(j)), radius * cos(lats(i)) * sin(lons(j)), radius * sin(lats(i)) ]
			quadFaces(i,j) = LocateFaceContainingPoint( sphere, xyz )
		enddo
	enddo

	write(filename,'(A,I1,A)') "cubedSphere_", initNest, ".m"
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

		write(WRITE_UNIT_1,'(A)',advance='NO') "quadFaces = ["
		do i = 1, 90
			do j = 1, 179
				write(WRITE_UNIT_1,'(I8,A)',advance='NO') quadFaces(i,j), ", "
			enddo
			write(WRITE_UNIT_1,'(I8,A)') quadFaces(i,180), "; "
		enddo
		do j = 1, 179
			write(WRITE_UNIT_1,'(I8,A)',advance='NO') quadFaces(91,j), ", "
		enddo
		write(WRITE_UNIT_1,'(I8,A)') quadFaces(91,180), "]; "

		write(WRITE_UNIT_1,'(A,I8,A)') "figure(1); clf; hold on; contour(lons, lats, quadFaces, ",&
				 sphere%faces%N_Active, "); colorbar;"
call LogMessage(exeLog, TRACE_LOGGING_LEVEL, "Test 2 : ", " complete.")


!
!----------------
! Test 3
!----------------
!
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

!
!----------------
! Test 4
!----------------
!
call LogMessage(exeLog, TRACE_LOGGING_LEVEL, "Test 4 : ", " writing .vtk output.")

write(filename,'(A,I1,A)') "cubedSphere_", initNest, ".vtk"
open(unit=WRITE_UNIT_1, file=filename, status='REPLACE',action='WRITE')
call WriteMeshToVTKPolyData( sphere, WRITE_UNIT_1, 'icosTriTest')

call LogMessage(exeLog, TRACE_LOGGING_LEVEL, "Test 4 : ", " complete.")
close(WRITE_UNIT_1)

!
!----------------
! PROGRAM END
!----------------
!
call Delete(sphere)
call Delete(exeLog)

end program 