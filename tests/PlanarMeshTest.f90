program planarMeshTester

use NumberKindsModule
use LoggerModule
use PolyMesh2dModule
use ParticlesModule
use EdgesModule
use FacesModule
use STDIntVectorModule

implicit none

integer(kint) :: i, j

type(PolyMesh2d) :: triMesh
type(PolyMesh2d) :: quadMesh
integer(kint) :: initNest
integer(kint) :: maxNest
integer(kint) :: amrLimit
real(kreal) :: ampFactor

type(Logger) :: exeLog

integer(kint), parameter :: nn = 101
real(kreal), parameter :: dx = 0.1
real(kreal), parameter :: xmin = -5.0_kreal
real(kreal), parameter :: xmax = 5.0_kreal
real(kreal), parameter :: ymin = xmin
real(kreal), parameter :: ymax = xmax
real(kreal) :: x(nn), xVec(3)
real(kreal) :: y(nn)
integer(kint) :: triFaces(nn,nn) = 0
integer(kint) :: quadFaces(nn,nn) = 0

integer(kint) :: randomVertex, randomFace
real(kreal) :: randReal
logical(klog) :: keepGoing
type(STDIntVector) :: edgesAroundFace, adjacentFaces, verticesAroundFace, facesAroundVertex
real(kreal) :: phys(3)

integer(kint) :: narg
character(len=100) :: arg, filename

call New(exeLog, DEBUG_LOGGING_LEVEL )
call LogMessage(exeLog, TRACE_LOGGING_LEVEL, "PlanarMeshTest: ", "build a planar mesh of triangles and quadrilaterals.")
call LogMessage(exeLog, TRACE_LOGGING_LEVEL, "PlanarMeshTest: ", "test point query algorithm in both meshes.")
call LogMessage(exeLog, TRACE_LOGGING_LEVEL, "PlanarMeshTest: ", "output to Matlab.")

narg = IARGC()
if (  narg == 0 ) then
	initNest = 0
else 
	call GETARG(1, arg)
	read(arg,*) initNest
endif

maxNest = initNest
amrLimit = 0
ampFactor = 3.0_kreal

do i = 1, nn
	x(i) = xmin + dx * (i-1)
	y(i) = ymin + dx * (i-1)
enddo

call LogMessage(exeLog, DEBUG_LOGGING_LEVEL, "PlanarMeshTest: ", "BUILDING TRIANGULAR MESH.")
call New(triMesh, TRI_HEX_SEED, initnest, maxnest, amrlimit, ampFactor)
call LogStats(triMesh, exeLog)

do i = 1, triMesh%particles%N
	if ( triMesh%particles%isActive(i) .AND. triMesh%particles%isPassive(i) ) then
		call LogMessage(exeLog, WARNING_LOGGING_LEVEL,"TriMesh Particles WARNING : both active and passive = .TRUE. at particle ", i)
	endif
enddo

!if ( procRank == 0 .AND. initNest <= 2 ) then
!	print *, "DEBUG : PRINTING ALL MESH INFO  "
!	call PrintDebugInfo( triMesh )
!endif
write(filename,'(A,I1,A)')"planeTriMeshTestMatlab",initNest,".m"
open(unit=WRITE_UNIT_1,file=filename,status='REPLACE')

	write(WRITE_UNIT_1,'(A)') "unifx = -5:0.1:5;"
	write(WRITE_UNIT_1,'(A)') "unify = unifx;"
	write(WRITE_UNIT_1,'(A)') "[xx, yy] = meshgrid(unifx, unify);"
	
	call WriteMeshToMatlab(triMesh, WRITE_UNIT_1)

	write(WRITE_UNIT_1,*) "figure(1);clf;hold on;"
	write(WRITE_UNIT_1,*) "plot(x,y,'k+');"
	write(WRITE_UNIT_1,*) "for i = 1:", triMesh%edges%N
	write(WRITE_UNIT_1,*) "		eX = [ x(edgeVerts(i,1)), x(edgeVerts(i,2))];"
	write(WRITE_UNIT_1,*) "     eY = [ y(edgeVerts(i,1)), y(edgeVerts(i,2))];"
	write(WRITE_UNIT_1,*) "     if edgeHasChildren(i) > 0 "
	write(WRITE_UNIT_1,*) "			plot(eX,eY,'b:x','LineWidth',2);"
	write(WRITE_UNIT_1,*) "     else"
	write(WRITE_UNIT_1,*) "     	plot(eX,eY,'k-o');"
	write(WRITE_UNIT_1,*) "		end "
	write(WRITE_UNIT_1,*) "end"
	write(WRITE_UNIT_1,*) "title('edges and particles');"
	write(WRITE_UNIT_1,*) "xlim([-5.1,5.1]);"
	write(WRITE_UNIT_1,*) "ylim([-5.1,5.1]);"

	write(WRITE_UNIT_1,*) "figure(2);clf;hold on;"
	write(WRITE_UNIT_1,*) "plot(x,y,'k+');"
	write(WRITE_UNIT_1,*) "for i = 1:", triMesh%faces%N
	write(WRITE_UNIT_1,*) "		fX = [ x(faceVerts(i,1)), x(faceVerts(i,2)), ..."
	write(WRITE_UNIT_1,*) "            x(faceVerts(i,3)), x(faceVerts(i,1))];"
	write(WRITE_UNIT_1,*) "     fY = [ y(faceVerts(i,1)), y(faceVerts(i,2)), ..."
	write(WRITE_UNIT_1,*) "            y(faceVerts(i,3)), y(faceVerts(i,1))];"
	write(WRITE_UNIT_1,*) "     if faceHasChildren(i) > 0 "
	write(WRITE_UNIT_1,*) "			plot(fX,fY,'b:x','LineWidth',2);"
	write(WRITE_UNIT_1,*) "     else"
	write(WRITE_UNIT_1,*) "     	plot(fX,fY,'k-o');"
	write(WRITE_UNIT_1,*) "		end "
	write(WRITE_UNIT_1,*) "		if faceCenterParticle(i) > 0 "
	write(WRITE_UNIT_1,*) "			plot(x(faceCenterParticle(i)),y(faceCenterParticle(i)),'ro');"
	write(WRITE_UNIT_1,*) "     end"
	write(WRITE_UNIT_1,*) "end"
	write(WRITE_UNIT_1,*) "title('faces and particles');"
	write(WRITE_UNIT_1,*) "xlim([-5.1,5.1]);"
	write(WRITE_UNIT_1,*) "ylim([-5.1,5.1]);"
	
	xvec = 0.0_kreal
	
	do i = 1, nn
		xvec(1) = x(i)
		do j = 1, nn
			xvec(2) = y(j)
			if ( pointIsOutsideMesh( triMesh, xvec ) ) &
				triFaces(j,i) = 1
		enddo
	enddo
	call LogMessage(exeLog, TRACE_LOGGING_LEVEL, "n points outside mesh = ", sum(triFaces))
	
	write(WRITE_UNIT_1,'(A)',advance='NO') "outsideTriMesh = ["
	do i = 1, nn - 1
		do j = 1, nn - 1
			write(WRITE_UNIT_1,'(I6,A)',advance='NO') triFaces(i,j), ", "
		enddo
		write(WRITE_UNIT_1,'(I6,A)') triFaces(i,nn), "; ..."
	enddo
	do j = 1, nn-1
		write(WRITE_UNIT_1,'(I6,A)',advance='NO') triFaces(nn,j), ", "
	enddo
	write(WRITE_UNIT_1,'(I6,A)') triFaces(nn,nn), "]; "
	
	write(WRITE_UNIT_1,*) "figure(3);clf;hold on;"
	write(WRITE_UNIT_1,*) "plot(x,y,'k+');"
	write(WRITE_UNIT_1,*) "for i = 1:", triMesh%edges%N
	write(WRITE_UNIT_1,*) "		eX = [ x(edgeVerts(i,1)), x(edgeVerts(i,2))];"
	write(WRITE_UNIT_1,*) "     eY = [ y(edgeVerts(i,1)), y(edgeVerts(i,2))];"
	write(WRITE_UNIT_1,*) "     if edgeHasChildren(i) > 0 "
	write(WRITE_UNIT_1,*) "			plot(eX,eY,'b:x','LineWidth',2);"
	write(WRITE_UNIT_1,*) "     else"
	write(WRITE_UNIT_1,*) "     	plot(eX,eY,'k-o');"
	write(WRITE_UNIT_1,*) "		end "
	write(WRITE_UNIT_1,*) "end"
	write(WRITE_UNIT_1,*) "scatter(reshape(xx, 101*101, 1), reshape(yy, 101*101, 1),..."
	write(WRITE_UNIT_1,*) "         12, reshape(outsideTriMesh,101*101,1));"
	write(WRITE_UNIT_1,*) "title('points outside mesh');"
	write(WRITE_UNIT_1,*) "xlim([-5.1,5.1]);"
	write(WRITE_UNIT_1,*) "ylim([-5.1,5.1]);"
	
	triFaces = 0
	do i = 1, nn
		xVec(1) = x(i)
		do j = 1, nn
			xVec(2) = y(j)
			triFaces(j,i) = locateFaceContainingPoint(triMesh, xVec)
		enddo
	enddo

	write(WRITE_UNIT_1,'(A)',advance='NO') "nearestFaceTriMesh = ["
	do i = 1, nn - 1
		do j = 1, nn - 1
			write(WRITE_UNIT_1,'(I6,A)',advance='NO') triFaces(i,j), ", "
		enddo
		write(WRITE_UNIT_1,'(I6,A)') triFaces(i,nn), "; ..."
	enddo
	do j = 1, nn-1
		write(WRITE_UNIT_1,'(I6,A)',advance='NO') triFaces(nn,j), ", "
	enddo
	write(WRITE_UNIT_1,'(I6,A)') triFaces(nn,nn), "]; "

	write(WRITE_UNIT_1,*) "figure(4);clf;hold on;"
	write(WRITE_UNIT_1,*) "plot(x,y,'k+');"
	write(WRITE_UNIT_1,*) "for i = 1:", triMesh%edges%N
	write(WRITE_UNIT_1,*) "		eX = [ x(edgeVerts(i,1)), x(edgeVerts(i,2))];"
	write(WRITE_UNIT_1,*) "     eY = [ y(edgeVerts(i,1)), y(edgeVerts(i,2))];"
	write(WRITE_UNIT_1,*) "     if edgeHasChildren(i) > 0 "
	write(WRITE_UNIT_1,*) "			plot(eX,eY,'b:x','LineWidth',2);"
	write(WRITE_UNIT_1,*) "     else"
	write(WRITE_UNIT_1,*) "     	plot(eX,eY,'k-o');"
	write(WRITE_UNIT_1,*) "		end "
	write(WRITE_UNIT_1,*) "end"
	write(WRITE_UNIT_1,*) "contourf(unifx, unify, nearestFaceTriMesh,", triMesh%faces%N_Active, ");"
	write(WRITE_UNIT_1,*) "title('nearest face');"
	write(WRITE_UNIT_1,*) "xlim([-5.1,5.1]);"
	write(WRITE_UNIT_1,*) "ylim([-5.1,5.1]);"
	
	keepGoing = .TRUE.
	call RANDOM_SEED()
	do while (keepGoing)
		call RANDOM_NUMBER(randReal)
		randomVertex = 1 + floor(triMesh%particles%N * randReal)
		if ( triMesh%particles%isPassive(randomVertex) .AND. .NOT. triMesh%particles%isActive(randomVertex) ) keepGoing = .FALSE.
	enddo
	
	keepGoing = .TRUE.
	do while (keepGoing)
		call RANDOM_NUMBER(randReal)
		randomFace = 1 + floor(triMesh%faces%N * randReal)
		if ( .NOT. triMesh%faces%hasChildren(randomFace) ) keepGoing = .FALSE.
	enddo
	
	call initialize(edgesAroundFace)
	call initialize(verticesAroundFace)
	call initialize(adjacentFaces)
	call initialize(facesAroundVertex)
	
	call CCWEdgesAroundFace( triMesh, edgesAroundFace, randomFace )
	
	call CCWVerticesAroundFace(triMesh,verticesAroundFace, randomFace)
	
	call CCWAdjacentFaces(triMesh, adjacentFaces, randomFace )
	
	call CCWFacesAroundVertex( triMesh, facesAroundVertex, randomVertex)	
	
	phys = PhysCoord(triMesh%particles, randomVertex)
	
	write(WRITE_UNIT_1,*) "randVert = [", phys(1), ", ", phys(2), "];"
	write(WRITE_UNIT_1,'(A)',advance='NO') "facesNearVert = ["
	do i = 1, facesAroundVertex%N -1
		if ( facesAroundVertex%int(i) > 0 ) then
		write(WRITE_UNIT_1,*) triMesh%particles%x( triMesh%faces%centerParticle(facesAroundVertex%int(i))), ", ",&
							  triMesh%particles%y( triMesh%faces%centerParticle(facesAroundVertex%int(i))), "; ..."
		endif
	enddo
	if ( facesAroundVertex%int(facesAroundVertex%N) > 0 ) then
	write(WRITE_UNIT_1,*) triMesh%particles%x( triMesh%faces%centerParticle(facesAroundVertex%int(facesAroundVertex%N))), ", ",&
						  triMesh%particles%y( triMesh%faces%centerParticle(facesAroundVertex%int(facesAroundVertex%N))), "];"
	endif
	
	phys = FaceCenterPhysCoord(triMesh%faces, randomFace, triMesh%particles)
	write(WRITE_UNIT_1,*) "randFace = [", phys(1), ", ", phys(2), "];"
	write(WRITE_UNIT_1,'(A)',advance='NO') "edgesAroundFaceX = ["
	do i = 1, edgesAroundFace%N - 1
		write(WRITE_UNIT_1,*) triMesh%particles%x( triMesh%edges%orig(edgesAroundFace%int(i))), ", ", &
							  triMesh%particles%x( triMesh%edges%dest(edgesAroundFace%int(i))), "; ..."
	enddo
	write(WRITE_UNIT_1,*) triMesh%particles%x( triMesh%edges%orig(edgesAroundFace%int(edgesAroundFace%N))), ", ", &
						  triMesh%particles%x( triMesh%edges%dest(edgesAroundFace%int(edgesAroundFace%N))), "];"
	write(WRITE_UNIT_1,'(A)',advance='NO') "edgesAroundFaceY = ["
	do i = 1, edgesAroundFace%N - 1
		write(WRITE_UNIT_1,*) triMesh%particles%y( triMesh%edges%orig(edgesAroundFace%int(i))), ", ", &
							  triMesh%particles%y( triMesh%edges%dest(edgesAroundFace%int(i))), "; ..."
	enddo
	write(WRITE_UNIT_1,*) triMesh%particles%y( triMesh%edges%orig(edgesAroundFace%int(edgesAroundFace%N))), ", ", &
						  triMesh%particles%y( triMesh%edges%dest(edgesAroundFace%int(edgesAroundFace%N))), "];"
	
	write(WRITE_UNIT_1,'(A)',advance='NO') "vertsAroundFace = ["
	do i = 1, verticesAroundFace%N - 1
		write(WRITE_UNIT_1,*) triMesh%particles%x( verticesAroundFace%int(i)), ", ",&
							  triMesh%particles%y(verticesAroundFace%int(i)), "; ..."
	enddo
	write(WRITE_UNIT_1,*) triMesh%particles%x( verticesAroundFace%int(verticesAroundFace%N)), ", ",&
						  triMesh%particles%y(verticesAroundFace%int(verticesAroundFace%N)), "];"
	write(WRITE_UNIT_1,'(A)',advance='NO') "adjacentFaces = ["
	do i = 1, adjacentFaces%N - 1
		if ( adjacentFaces%int(i) > 0 ) then
			phys = FaceCenterPhysCoord(triMesh%faces, adjacentFaces%int(i), triMesh%particles)
			write(WRITE_UNIT_1,*) phys(1), ", ", phys(2), "; ..."
		endif
	enddo
	if ( adjacentFaces%int(adjacentFaces%N) > 0 ) then
		phys = FaceCenterPhysCoord(triMesh%faces, adjacentFaces%int(adjacentFaces%N), triMesh%particles)
		write(WRITE_UNIT_1,*) phys(1), ", ", phys(2), "];"
	endif

close(WRITE_UNIT_1)

call Delete(triMesh)

call LogMessage(exeLog, DEBUG_LOGGING_LEVEL, "PlanarMeshTest: ", "BUILDING QUADRILATERAL MESH.")
call New(quadMesh, QUAD_RECT_SEED, initNest, maxNest, amrLimit, ampFactor)
call LogStats(quadMesh, exeLog)

do i = 1, quadMesh%particles%N
	if ( quadMesh%particles%isActive(i) .AND. quadMesh%particles%isPassive(i) ) then
		call LogMessage(exeLog, WARNING_LOGGING_LEVEL,"QuadMesh Particles WARNING : both active and passive = .TRUE. at particle ", i)
	endif
enddo

!if ( procRank == 0 .AND. initNest <= 2 ) then
!	print *, "DEBUG : PRINTING ALL MESH INFO  "
!	call PrintDebugInfo( quadMesh )
!endif

write(filename,'(A,I1,A)')"planeQuadMeshTestMatlab",initNest,".m"
open(unit=WRITE_UNIT_1,file=filename,status='REPLACE')

	call WriteMeshToMatlab(quadMesh, WRITE_UNIT_1)

	write(WRITE_UNIT_1,*) "figure(5);clf;hold on;"
	write(WRITE_UNIT_1,*) "plot(x,y,'k+');"
	write(WRITE_UNIT_1,*) "for i = 1:", quadMesh%edges%N
	write(WRITE_UNIT_1,*) "		eX = [ x(edgeVerts(i,1)), x(edgeVerts(i,2))];"
	write(WRITE_UNIT_1,*) "     eY = [ y(edgeVerts(i,1)), y(edgeVerts(i,2))];"
	write(WRITE_UNIT_1,*) "     if edgeHasChildren(i) > 0 "
	write(WRITE_UNIT_1,*) "			%plot(eX,eY,'b:x','LineWidth',2);"
	write(WRITE_UNIT_1,*) "     else"
	write(WRITE_UNIT_1,*) "     	plot(eX,eY,'k-o');"
	write(WRITE_UNIT_1,*) "		end "
	write(WRITE_UNIT_1,*) "end"
	write(WRITE_UNIT_1,*) "title('edges and particles');"
	write(WRITE_UNIT_1,*) "xlim([-5.1,5.1]);"
	write(WRITE_UNIT_1,*) "ylim([-5.1,5.1]);"

	write(WRITE_UNIT_1,*) "figure(6);clf;hold on;"
	write(WRITE_UNIT_1,*) "plot(x,y,'k+');"
	write(WRITE_UNIT_1,*) "for i = 1:", quadMesh%faces%N
	write(WRITE_UNIT_1,*) "		fX = [ x(faceVerts(i,1)), x(faceVerts(i,2)), x(faceVerts(i,3)),..."
	write(WRITE_UNIT_1,*) "            x(faceVerts(i,4)), x(faceVerts(i,1))];"
	write(WRITE_UNIT_1,*) "     fY = [ y(faceVerts(i,1)), y(faceVerts(i,2)), y(faceVerts(i,3)),..."
	write(WRITE_UNIT_1,*) "            y(faceVerts(i,4)), y(faceVerts(i,1))];"
	write(WRITE_UNIT_1,*) "     if faceHasChildren(i) > 0 "
	write(WRITE_UNIT_1,*) "			%plot(fX,fY,'b:x','LineWidth',2);"
	write(WRITE_UNIT_1,*) "     else"
	write(WRITE_UNIT_1,*) "     	plot(fX,fY,'k-o');"
	write(WRITE_UNIT_1,*) "		end "
	write(WRITE_UNIT_1,*) "		if faceCenterParticle(i) > 0 "
	write(WRITE_UNIT_1,*) "			plot(x(faceCenterParticle(i)),y(faceCenterParticle(i)),'ro');"
	write(WRITE_UNIT_1,*) "     end"
	write(WRITE_UNIT_1,*) "end"
	write(WRITE_UNIT_1,*) "title('faces and particles');"
	write(WRITE_UNIT_1,*) "xlim([-5.1,5.1]);"
	write(WRITE_UNIT_1,*) "ylim([-5.1,5.1]);"
	
	do i = 1, nn
		xvec(1) = x(i)
		do j = 1, nn
			xvec(2) = y(j)
			if ( pointIsOutsideMesh( quadMesh, xvec ) ) &
				quadFaces(j,i) = 1
		enddo
	enddo
	call LogMessage(exeLog, TRACE_LOGGING_LEVEL, "n points outside mesh = ", sum(quadFaces))
	
	write(WRITE_UNIT_1,'(A)',advance='NO') "outsideQuadMesh = ["
	do i = 1, nn - 1
		do j = 1, nn - 1
			write(WRITE_UNIT_1,'(I6,A)',advance='NO') quadFaces(i,j), ", "
		enddo
		write(WRITE_UNIT_1,'(I6,A)') quadFaces(i,nn), "; ..."
	enddo
	do j = 1, nn-1
		write(WRITE_UNIT_1,'(I6,A)',advance='NO') quadFaces(nn,j), ", "
	enddo
	write(WRITE_UNIT_1,'(I6,A)') quadFaces(nn,nn), "]; "

	write(WRITE_UNIT_1,*) "figure(7);clf;hold on;"
	write(WRITE_UNIT_1,*) "plot(x,y,'k+');"
	write(WRITE_UNIT_1,*) "for i = 1:", quadMesh%edges%N
	write(WRITE_UNIT_1,*) "		eX = [ x(edgeVerts(i,1)), x(edgeVerts(i,2))];"
	write(WRITE_UNIT_1,*) "     eY = [ y(edgeVerts(i,1)), y(edgeVerts(i,2))];"
	write(WRITE_UNIT_1,*) "     if edgeHasChildren(i) > 0 "
	write(WRITE_UNIT_1,*) "			plot(eX,eY,'b:x','LineWidth',2);"
	write(WRITE_UNIT_1,*) "     else"
	write(WRITE_UNIT_1,*) "     	plot(eX,eY,'k-o');"
	write(WRITE_UNIT_1,*) "		end "
	write(WRITE_UNIT_1,*) "end"
	write(WRITE_UNIT_1,*) "scatter(reshape(xx, 101*101, 1), reshape(yy, 101*101, 1),..."
	write(WRITE_UNIT_1,*) "         12, reshape(outsideQuadMesh,101*101,1));"
	write(WRITE_UNIT_1,*) "title('points outside mesh');"
	write(WRITE_UNIT_1,*) "xlim([-5.1,5.1]);"
	write(WRITE_UNIT_1,*) "ylim([-5.1,5.1]);"
	
	quadFaces = 0
	do i = 1, nn
		xVec(1) = x(i)
		do j = 1, nn
			xVec(2) = y(j)
			quadFaces(j,i) = locateFaceContainingPoint(quadMesh, xVec)
		enddo
	enddo

	write(WRITE_UNIT_1,'(A)',advance='NO') "nearestFaceQuadMesh = ["
	do i = 1, nn - 1
		do j = 1, nn - 1
			write(WRITE_UNIT_1,'(I6,A)',advance='NO') quadFaces(i,j), ", "
		enddo
		write(WRITE_UNIT_1,'(I6,A)') quadFaces(i,nn), "; ..."
	enddo
	do j = 1, nn-1
		write(WRITE_UNIT_1,'(I6,A)',advance='NO') quadFaces(nn,j), ", "
	enddo
	write(WRITE_UNIT_1,'(I6,A)') quadFaces(nn,nn), "]; "

	write(WRITE_UNIT_1,*) "figure(8);clf;hold on;"
	write(WRITE_UNIT_1,*) "plot(x,y,'k+');"
	write(WRITE_UNIT_1,*) "for i = 1:", quadMesh%edges%N
	write(WRITE_UNIT_1,*) "		eX = [ x(edgeVerts(i,1)), x(edgeVerts(i,2))];"
	write(WRITE_UNIT_1,*) "     eY = [ y(edgeVerts(i,1)), y(edgeVerts(i,2))];"
	write(WRITE_UNIT_1,*) "     if edgeHasChildren(i) > 0 "
	write(WRITE_UNIT_1,*) "			plot(eX,eY,'b:x','LineWidth',2);"
	write(WRITE_UNIT_1,*) "     else"
	write(WRITE_UNIT_1,*) "     	plot(eX,eY,'k-o');"
	write(WRITE_UNIT_1,*) "		end "
	write(WRITE_UNIT_1,*) "end"
	write(WRITE_UNIT_1,*) "contourf(unifx, unify, nearestFaceQuadMesh,", quadMesh%faces%N_Active, ");"
	write(WRITE_UNIT_1,*) "title('nearest face');"
	write(WRITE_UNIT_1,*) "xlim([-5.1,5.1]);"
	write(WRITE_UNIT_1,*) "ylim([-5.1,5.1]);"
	
close(WRITE_UNIT_1)

call Delete(quadMesh)


call Delete(exeLog)
end program 