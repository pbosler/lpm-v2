program cubedSphereTest

use NumberKindsModule
use LoggerModule
use PolyMeshOOModule
use ParticlesOOModule
use EdgesOOModule
use FacesOOModule
use SphereGeomModule
use STDIntVectorModule
use OutputWriterModule

#include "LpmConfig.h"

type(PolyMesh2d) :: mesh
integer(kint) :: initNest
integer(kint), parameter :: amrLimit = 0
real(kreal), parameter :: sphRadius = one
integer(kint) :: nedges

type(Logger) :: exeLog

real(kreal), dimension(180) :: lons
real(kreal), dimension(91) :: lats
real(kreal), parameter :: dlambda = 2.0_kreal * PI / 180.0_kreal
real(kreal), dimension(3) :: qpt

integer(kint), dimension(91,180) :: inFace

integer(kint) :: randomVertex
integer(kint) :: randomFace
logical(klog) :: keepGoing
real(kreal) :: randReal
type(STDIntVector) :: edgesAroundFace, adjacentFaces
real(kreal) :: xyz(3)

integer(kint) :: i, j

character(len=MAX_STRING_LENGTH) :: logString
character(len=MAX_STRING_LENGTH) :: fname

call New(exeLog, DEBUG_LOGGING_LEVEL)
call getInput(initNest)

call LogMessage(exelog, TRACE_LOGGING_LEVEL, "CubedSphereMeshTest part 1:", " build an icosahedral triangle sphere mesh.")
call LogMessage(exelog, TRACE_LOGGING_LEVEL, "CubedSphereMeshTest part 2:", " test point query algorithm.")
call LogMessage(exeLog, TRACE_LOGGING_LEVEL, "CubedSphereMeshTest part 3:", " test adjacency algorithms.")
call LogMessage(exelog, TRACE_LOGGING_LEVEL, "CubedSphereMeshTest part 4:", " output to paraview/vtk.")

call mesh%init("cubed_sphere", initNest, initNest, amrLimit, sphRadius)
call mesh%logStats(exeLog)

!call mesh%edges%writeMatlab(6)

do j=1,180
    lons(j) = dlambda *(j-1)
enddo
do i=1,91
    lats(i) = -0.5_kreal * PI + dlambda*(i-1)
enddo

do j=1,180
    do i=1,91
        qpt = [cos(lats(i))*cos(lons(j)), cos(lats(i))*sin(lons(j)), sin(lats(i))]
        inFace(i,j) = mesh%locatePointInFace(qpt)
    enddo
enddo
print *, "Uniform mesh points located."

qpt = [cos(-0.25_kreal)*cos(2.95_kreal),cos(-0.25_kreal)*sin(2.95_kreal), sin(-0.25_kreal)]
randomFace = mesh%locatePointInFace(qpt)
print *, "query point located in face ", randomFace

edgesAroundFace = mesh%ccwEdgesAroundFace(randomFace)
print *, "edgesAroundFace ", randomFace, ":"
call edgesAroundFace%print()
adjacentFaces = mesh%ccwAdjacentFaces(randomFace)
print *, "adjacentFaces to ", randomFace, ":"
call adjacentFaces%print()

write(fname,'(A,I1,A)') "cubedSphereTest", initNest, ".m"
open(unit=WRITE_UNIT_1, file=trim(fname), action='write', status='replace')

call writeToMatlab(lons, WRITE_UNIT_1, 'lons')
call writeToMatlab(lats, WRITE_UNIT_1, 'lats')
call writeToMatlab(inFace, WRITE_UNIT_1, 'inFace')
write(WRITE_UNIT_1, '(A,2G15.8,A)') "queryPt = [", 2.95_kreal, -0.25_kreal, "];"
write(WRITE_UNIT_1, '(A,2G15.8,A)') "qInFace = [", longitude(mesh%faces%physCentroid(randomFace, mesH%particles)), &
    latitude(mesH%faces%physCentroid(randomFace, mesh%particles)), "];"
write(WRITE_UNIT_1, '(A)') "edgesAroundFaceLon = ["
do i=1, edgesAroundFace%N
    xyz = mesh%particles%physCoord(mesh%edges%orig(edgesAroundFace%int(i)))
    qpt = mesh%particles%physCoord(mesh%edges%dest(edgesAroundFace%int(i)))
    write(WRITE_UNIT_1,'((G15.8,A))') Longitude(xyz), ", ..." 
    write(WRITE_UNIT_1,'(G15.8,A)') longitude(qpt), ", ..."
enddo
xyz = mesh%particles%physCoord(mesh%edges%orig(edgesAroundFace%int(1)))
qpt = mesh%particles%physCoord(mesh%edges%dest(edgesAroundFace%int(1)))
write(WRITE_UNIT_1,'((G15.8,A))') Longitude(xyz), ", ... "
write(WRITE_UNIT_1,'(G15.8,A)') longitude(qpt), "]; "
write(WRITE_UNIT_1, '(A)') "edgesAroundFaceLat = ["
do i=1, edgesAroundFace%N
    xyz = mesh%particles%physCoord(mesh%edges%orig(edgesAroundFace%int(i)))
    qpt = mesh%particles%physCoord(mesh%edges%dest(edgesAroundFace%int(i)))
    write(WRITE_UNIT_1,'((G15.8,A))') latitude(xyz), ", ..." 
    write(WRITE_UNIT_1,'(G15.8,A)') latitude(qpt), ", ..."
enddo
xyz = mesh%particles%physCoord(mesh%edges%orig(edgesAroundFace%int(1)))
qpt = mesh%particles%physCoord(mesh%edges%dest(edgesAroundFace%int(1)))
write(WRITE_UNIT_1,'((G15.8,A))') latitude(xyz), ", ..." 
write(WRITE_UNIT_1,'(G15.8,A)') latitude(qpt), "];"
write(WRITE_UNIT_1, '(A)', advance='no') "adjFaces = ["
do i=1, adjacentFaces%N-1
    xyz = mesh%faces%physCentroid(adjacentFaces%int(i), mesh%particles)
    write(WRITE_UNIT_1,'(2(G15.8,A))') Longitude(xyz), ", ", latitude(xyz), "; ..."
enddo
xyz = mesh%faces%physCentroid(adjacentFaces%int(adjacentFaces%N), mesh%particles)
write(WRITE_UNIT_1,'(2(G15.8,A))') Longitude(xyz), ", ", latitude(xyz), "];"

write(WRITE_UNIT_1,'(A,I8,A)') "figure(1);clf; hold on; contour(lons, lats, inFace,", mesh%faces%N_Active, "); colorbar;"
write(WRITE_UNIT_1,'(A)') "plot(queryPt(1), queryPt(2),'k.','MarkerSize',20)"
write(WRITE_UNIT_1,'(A)') "plot(qInFace(1), qInFace(2), 'bo', 'MarkerSize',10)"
write(WRITE_UNIT_1,'(A)') "plot(edgesAroundFaceLon, edgesAroundFaceLat,'r-','Linewidth',2)"

close(WRITE_UNIT_1)

write(fname,'(A,I1,A)') "cubedSphereTest", initNest, ".vtp"
open(unit=WRITE_UNIT_1, file=trim(fname), action='write', status='replace')
call mesh%writeVTKSerialStartXML(WRITE_UNIT_1)
call mesh%writeVTKSerialEndXML(WRITE_UNIT_1)
close(WRITE_UNIT_1)

call LogMessage(exeLog, TRACE_LOGGING_LEVEL, "result: ", "PASSED")

contains

subroutine getInput(initNest)
    integer(kint), intent(out) :: initNest
    integer(kint) :: argc
    character(len=MAX_STRING_LENGTH) :: argv
    argc = IARGC()
    if (argc==0) then
        initNest = 0
    else
        call GETARG(1, argv)
        read(argv,*) initNest
    endif
end subroutine

end program