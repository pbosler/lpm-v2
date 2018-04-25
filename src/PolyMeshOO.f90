module PolyMeshOOModule
!> @file PolyMesh2d.f90
!> Data structure and methods to manage 2d polygonal meshes in the plane or on 2d manifolds embedded in R3.
!>
!> @author Peter Bosler, Sandia National Laboratories Center for Computing Research
!>
!>
!> @defgroup PolyMesh2d PolyMesh2d
!> @brief Each polygonal mesh is made up of @ref Particles, @ref Edges, and @ref Faces.
!>
!> All meshes start with an initial "seed" that is recursively refined to generate the @ref Faces quadtree, its @ref Edges,
!> and the associated @ref Particles.
!> Mesh seeds are supplied in  `*Seed.dat` files.
!>
!>
!>  Current seeds for __planar meshes__ include a mesh of quadrilateral faces from a square or a mesh of triangular faces from
!>	a hexagon, both with free boundaries, and a quadrilateral seed that is periodic in the x-direction and free in the y-direction.
!>	Users choose the initial mesh from one of the predefined seeds, defined to approximately cover the domain
!>	@f$ D = \{ (x,y) : x,y\in [-1,1] \} @f$, depending on the chosen seed.
!>	The seed is then multiplied by a user-supplied scalar `maxR` to cover larger subsets of the plane.
!>
!>  Seeds for __spherical meshes__ include a mesh of quadrilaterals from the cubed sphere or a mesh of triangles from
!>	an icosahedral triangulation of the sphere.
!>	Users choose one of the spherical seeds (which discretize a unit sphere) and may choose a different radius.
!>
!>  Seeds are specified in the code by using the correct identifying integer parameter defined in @ref NumberKinds.
!>
!>	The mesh seed is recursively refined until a user-supplied initial depth is achieved in the @ref Faces quadtree.
!>	This refinement increases spatial resolution.
!>
!>	Adaptive refinement depends on the context of each problem and is handled by subclasses that include relevant @ref Field objects
!>	for storing and handling physical data.
!>  These methods are provided in the @ref Refinement module.
!>
!>	Some currently implemented mesh seeds are illustrated below. @n
!>	Particles are drawn as circles and annotated by @f$p_i@f$ where the subscript i indicates
!>	the index of that particle in the Particles structure.@n
!>	Edges are illustrated as arrows @f$ e_i @f$ and Faces by @f$ F_i @f$.
!>
!>	Seeds are recursively refined by either facesmodule::dividetriface or facesmodule::dividequadface
!>	until a desired spatial resolution is achieved.
!>
!>  Point-location queries are provided that use a tree search to accelerate a walk search using the polymesh2dmodule::locatefacecontainingpoint function.
!>
!>	@image html QuadRectSeed.png "Seed for an initially square mesh of quadrilateral faces."
!>	@image latex QuadRectSeed.eps "Seed for an initially square mesh of quadrilateral faces." width=5in
!>
!>	@image html TriHexSeed.png "Seed for an initially hexagonal mesh of triangular faces."
!>	@image latex TriHexSeed.eps "Seed for an initially hexagonal mesh of triangular faces." width=5in
!>
!>  @image html cubedSphereSeed.png "Seed for a cubed sphere mesh of quadrilateral faces."
!>	@image latex cubedSphereSeed.eps "Seed for a cubed sphere mesh of quadrilateral faces." width=5in
!>
!>	@image html icosTriSeed.png "Seed for an icosahedral mesh of triangular faces on the sphere."
!>	@image latex icosTriSeed.eps "Seed for an icosahedral mesh of triangular faces on the sphere." width=5in
!>
!> @{

use NumberKindsModule
use UtilitiesModule
use STDIntVectorModule
use OutputWriterModule
use LoggerModule
use ParticlesOOModule
use EdgesOOModule
use FacesOOModule
use PlaneGeomModule
use SphereGeomModule
use FieldOOModule
use CubicGLLModule

#include "LpmConfig.h"

implicit none
private

public PolyMesh2d, meshKindFromString

!> @brief Combines mesh primitives into one complete mesh object. Extend this type by adding fields for applications.
type PolyMesh2d
    integer(kint) :: meshSeed !< geometry kind (e.g., planar or spherical) as defined in @ref NumberKinds.  Seeds define the initial, unrefined mesh
    class(Particles), pointer :: particles => null() !< @ref Particles defining vertices, edges, and interiors of mesh faces and edges
    class(Edges), pointer :: edges => null() !< @ref Edges connecting origin and destination particles for each edge, and each edge's left and right face
    class(Faces), pointer :: faces => null() !< @ref Faces
    integer(kint) :: faceKind = 0 !< face kind (e.g., triangular or quadrilateral) as defined in @ref NumberKinds
    integer(kint) :: geomKind = 0 !< geometry kind (e.g., planar or spherical) as defined in @ref NumberKinds
    integer(kint) :: initNest = 0 !< initial refinement level (in terms of the @ref Faces quadtree) for a mesh
    integer(kint) :: maxNest = 0 !< maximum tree Faces tree depth allowed. Used for memory allocation.
    integer(kint) :: amrLimit = 0 !< Number of times a face may be divided beyond the initial refinement level
    integer(kint) :: nRootFaces = 0 !< Number of faces in the mesh seed
    real(kreal) :: t = dzero !< time value represented by the mesh
    real(kreal) :: maxRadius = 1.0_kreal

    contains
        procedure :: init
        procedure :: logStats
        procedure :: getSeed
        procedure, private :: readSeedFile
        procedure, private :: nVerticesInMesh
        procedure, private :: nEdgesInMesh
        procedure, private :: nFacesInMesh
        procedure :: writeMatlab
        procedure :: writeVTKSerialStartXML
        procedure :: writePlanarVTKStartWithHeightFromField
        procedure :: writeLagCoordsToVTKPointDataXML
        procedure :: writeVTKSerialEndXML
        procedure :: integrateScalar
        procedure :: interpolateScalar
        procedure :: setLagCoordFromOther
        procedure :: interpolateLagCoord
        procedure :: locatePointInFace
        procedure, private :: locatePointTree
        procedure, private :: locatePointWalk
        procedure, private :: nearestRootFace
        procedure :: ccwEdgesAroundFace
        procedure :: ccwAdjacentFaces
        procedure :: pointIsOutsideMesh
        procedure :: copy
        procedure :: meshSize
!        procedure :: refine

        final :: deleteMesh
end type

logical(klog), save :: logInit = .FALSE.
type(Logger) :: log
character(len=28), save :: logKey = 'PolyMesh2d'
integer(kint), parameter :: logLevel = DEBUG_LOGGING_LEVEL
character(len=MAX_STRING_LENGTH) :: logstring

contains

!@>brief Returns the integer defined in @ref NumberKindsModule that corresponds to a mesh_type string.
function meshKindFromString(str)
    integer(kint) :: meshKindFromString
    character(len=*), intent(in) :: str

    meshKindFromString = 0
    if (trim(str) == "planar_tri") then
        meshKindFromString = TRI_HEX_SEED
    elseif (trim(str) == "planar_quad") then
        meshKindFromString = QUAD_RECT_SEED
    elseif (trim(str) == "planar_cubic_quad") then
        meshKindFromString = CUBIC_PLANE_SEED
    elseif (trim(str) == "cubed_sphere") then
        meshKindFromString = CUBED_SPHERE_SEED
    elseif (trim(str) == "icos_tri_sphere") then
        meshKindFromString = ICOS_TRI_SPHERE_SEED
    else
        call LogMessage(log, ERROR_LOGGING_LEVEL, "invalid mesh_type string: ", str)
    endif
end function

!> @brief Allocates memory for a new mesh (@ref Particles, @ref Edges, and @ref Faces).
!> Initializes a new mesh using a "mesh seed," which is recursively refined until a desired spatial resolution is achieved.
!> The recursive refinement defines the @ref Faces quadtree and the @ref Edges binary tree.
subroutine init(self, mesh_type, initNest, maxNest, amrLimit, maxRadius)
    class(PolyMesh2d), intent(inout) :: self
    character(len=*), intent(in) ::  mesh_type
    integer(kint), intent(in) :: initNest, maxNest, amrLimit
    real(kreal), intent(in) :: maxRadius
    !
    integer(kint) :: i,j, nFacesOld, startIndex

    if ( .NOT. logInit) call InitLogger(log, procRank)

    if (associated(self%particles)) then
        call LogMessage(log, ERROR_LOGGING_LEVEL, trim(logkey)//" ERROR:", &
            " init called on mesh with particles already associated.")
        return
    endif
    if (associated(self%edges)) then
        call LogMessage(log, ERROR_LOGGING_LEVEL, trim(logkey)//" ERROR:", &
            " init called on mesh with particles already associated.")
        return
    endif
    if (associated(self%faces)) then
        call LogMessage(log, ERROR_LOGGING_LEVEL, trim(logkey)//" ERROR:", &
            " init called on mesh with faces already associated.")
        return
    endif

    self%meshSeed = meshKindFromString(mesh_type)
    self%initNest = initNest
    self%maxNest = maxNest
    self%amrLimit = amrLimit
    self%maxRadius = maxRadius

    self%t = dzero

    call self%getSeed()

    startIndex = 1
    do i=1, initNest
        nFacesOld = self%faces%N
        do j=startIndex, nFacesOld
            call self%faces%divide(j, self%particles, self%edges)
        enddo
        startIndex = nFacesOld + 1
    enddo
end subroutine

subroutine copy(self, other)
    class(PolyMesh2d), intent(inout) :: self
    class(PolyMesh2d), intent(in) :: other
    
    self%meshSeed = other%meshSeed
    self%initNest = other%initNest
    self%maxNest = other%maxNest
    self%amrLimit = other%amrLimit
    self%maxRadius = other%maxRadius
    self%t = other%t
    self%nRootFaces = other%nRootFaces
    self%geomKind = other%geomKind
    self%faceKind = other%faceKind
    
    call self%particles%copy(other%particles)
    call self%edges%copy(other%edges)
    call self%faces%copy(other%faces)
end subroutine

!> @brief Integrates a scalar field over a mesh
!> @param self Target PolyMesh2d object
!> @param scalar Field to be integrated (from @ref FieldsOOModule)
pure function integrateScalar(self, scalar)
    real(kreal) :: integrateScalar
    class(PolyMesh2d), intent(in) ::  self
    class(Field), intent(in) :: scalar
    !
    integer(kint) :: i, j, pindex

    integrateScalar = dzero
    if (self%faceKind == TRI_PANEL .or. self%faceKind == QUAD_PANEL) then
        do i=1, self%faces%N
            if (.not. self%faces%hasChildren(i)) then
                pindex = self%faces%centerParticles(1,i)
                integrateScalar = integrateScalar + self%particles%weight(pindex) * scalar%comp1(pindex)
            endif
        enddo
    elseif (self%faceKind == QUAD_CUBIC_PANEL) then
        do i=1, self%faces%N
            if (.not. self%faces%hasChildren(i)) then
                do j=1,12
                    integrateScalar = integrateScalar + self%particles%weight(self%faces%vertices(j,i)) * &
                        scalar%comp1(self%faces%vertices(j,i))
                enddo
                do j=1,4
                    integrateScalar = integrateScalar + self%particles%weight(self%faces%centerParticles(j,i)) * &
                        scalar%comp1(self%faces%centerParticles(j,i))
                enddo
            endif
        enddo
    endif
end function

function pointIsOutsideMesh(self, xyz, faceIndex)
    logical(klog) :: pointIsOutsideMesh
    class(PolyMesh2d), intent(in) :: self
    real(kreal), dimension(3), intent(in) :: xyz
    integer(kint), intent(in), optional :: faceIndex
    !
    type(STDIntVector) :: faceEdges
    type(STDIntVector) :: boundaryEdges
    integer(kint) :: i, faceId
    real(kreal), dimension(3) :: centroid, p, q, v1, reflection
    real(kreal) :: insideDist, outsideDist
    
    pointIsOutsideMesh = .FALSE.
    if (self%geomKind == PLANAR_GEOM) then
        if (present(faceIndex)) then
            faceId = faceIndex
        else
            faceId = self%locatePointInFace(xyz)
        endif
        faceEdges = self%ccwEdgesAroundFace(faceId)
        call initialize(boundaryEdges)
        do i=1, faceEdges%N
            if (self%edges%onBoundary(faceEdges%int(i))) &
                call boundaryEdges%pushBack(faceEdges%int(i))
        enddo
        if (boundaryEdges%n > 0) then
            v1 = dzero
            centroid = self%faces%physCentroid(faceId, self%particles)
            insideDist = sqrt(sum((centroid-xyz)*(centroid-xyz)))
            do i=1, boundaryEdges%n
                q = self%edges%vector(boundaryEdges%int(i), self%particles)
                if (self%edges%positiveOrientation(boundaryEdges%int(i), faceId)) then
                    v1 = self%particles%physCoord(self%edges%orig(boundaryEdges%int(i)))
                else
                    v1 = self%particles%physCoord(self%edges%dest(boundaryEdges%int(i)))
                    q = -q
                endif
                q = q/sqrt(sum(q*q))
                p = centroid - v1
                reflection = centroid - 2.0_kreal * (p - sum(p*q)*q)
                outsideDist = sqrt(sum((reflection-xyz)*(reflection-xyz)))
                if (outsideDist < insideDist) &
                    pointIsOutsideMesh = .TRUE.
            enddo
        endif
    endif
end function

function interpolateScalar(self, xyz, scalar)
    real(kreal) :: interpolateScalar
    class(PolyMesh2d), intent(in)  :: self
    real(kreal), dimension(3), intent(in) :: xyz
    class(Field), intent(in) :: scalar
    !
    integer(kint) :: faceId, i, j
    real(kreal), dimension(3,4) :: faceCorners
    real(kreal), dimension(2) :: refCrd
    real(kreal), dimension(12) :: vertvals
    real(kreal), dimension(4) :: cntrVals

    interpolateScalar = FILL_VALUE
    if (self%faceKind == QUAD_CUBIC_PANEL) then
        if (self%geomKind == PLANAR_GEOM) then
            faceId = self%locatePointInFace(xyz)
            if (self%pointIsOutsideMesh(xyz, faceId)) then
                interpolateScalar=dzero
            else
                do i=1,4
                    faceCorners(:,i) = self%particles%physCoord(self%faces%vertices(mod(3*i+9,12)+1,faceId))
                enddo

                refCrd = invertBilinear(faceCorners, xyz)

                do i=1,12
                    vertvals(i) = scalar%comp1(self%faces%vertices(i,faceId))
                enddo
                do i=1,4
                    cntrVals(i) = scalar%comp1(self%faces%centerParticles(i,faceId))
                enddo
                interpolateScalar = bivariateScalarInterpolation(refCrd, vertVals, cntrVals)
            endif
        endif
    endif
end function

subroutine setLagCoordFromOther(self, other)
    class(PolyMesh2d), intent(inout) :: self
    class(PolyMesh2d), intent(in) :: other
    !
    real(kreal), dimension(3) :: lagX, physX
    integer(kint) :: i
    
    do i=1, self%particles%n
        physX = self%particles%physCoord(i)
        lagX = other%interpolateLagCoord(physX)
        self%particles%x0(i) = lagX(1)
        self%particles%y0(i) = lagX(2)
        if (self%geomKind == SPHERE_GEOM) &
            self%particles%z0(i) = lagX(3)
    enddo
end subroutine


function interpolateLagCoord(self, xyz)
    real(kreal), dimension(3) :: interpolateLagCoord
    class(PolyMesh2d), intent(in) :: self
    real(kreal), dimension(3), intent(in) :: xyz
    !
    integer(kint) :: faceId, i, j
    real(kreal), dimension(3,4) :: faceCorners
    real(kreal), dimension(2) :: refCrd
    real(kreal), dimension(12) :: vertX0, vertY0
    real(kreal), dimension(4) :: cntrX0, cntrY0
    
    interpolateLagCoord = FILL_VALUE
    if (self%faceKind == QUAD_CUBIC_PANEL) then
        if (self%geomKind == PLANAR_GEOM) then
            faceId = self%locatePointInFace(xyz)
            interpolateLagCoord(3) = dzero
            if (.not. self%pointIsOutsideMesh(xyz, faceId)) then
                do i=1,4
                    faceCorners(:,i) = self%particles%physCoord(self%faces%vertices(mod(3*i+9,12)+1,faceId))
                enddo
                refCrd = invertBilinear(faceCorners, xyz)

                do i=1,12
                    vertX0(i) = self%particles%x0(self%faces%vertices(i,faceId))
                    vertY0(i) = self%particles%y0(self%faces%vertices(i,faceId))
                enddo
                do i=1,4
                    cntrX0(i) = self%particles%x0(self%faces%centerParticles(i,faceId))
                    cntrY0(i) = self%particles%y0(self%faces%centerParticles(i,faceId))
                enddo
                interpolateLagCoord(1) = bivariateScalarInterpolation(refCrd, vertX0, cntrX0)
                interpolateLagCoord(2) = bivariateScalarInterpolation(refCrd, vertY0, cntrY0)
            endif
        endif
    endif
end function

pure function meshSize(self)
    real(kreal) :: meshSize
    class(PolyMesh2d), intent(in) :: self
    meshSize = self%edges%avgLength(self%particles)
end function

!> @brief Outputs basic mesh info to a @ref Logger object.
!> @param[in] self Target mesh
!> @param[inout] aLog logger
subroutine logStats(self, aLog)
    class(PolyMesh2d), intent(in) :: self
    type(Logger), intent(inout) :: aLog
    !
    call LogMessage(aLog, TRACE_LOGGING_LEVEL, logKey, " PolyMesh2d stats:")
    call StartSection(aLog)
    call LogMessage(aLog, TRACE_LOGGING_LEVEL, "meshSeed = ", self%meshSeed)
    call LogMessage(aLog, TRACE_LOGGING_LEVEL, "initNest = ", self%initNest)
    call LogMessage(aLog, TRACE_LOGGING_LEVEL, "maxNest = ", self%maxNest)
    call LogMessage(alog, TRACE_LOGGING_LEVEL, "amrLimit = ", self%amrLimit)
    call LogMessage(alog, TRACE_LOGGING_LEVEL, "maxRadius = ", self%maxRadius)
    call LogMessage(aLog, TRACE_LOGGING_LEVEL, "geomKind = ", geomString(self%geomKind))
    call LogMessage(aLog, TRACE_LOGGING_LEVEL, "faceKind = ", faceString(self%faceKind))
    call self%particles%logStats(aLog)
    call self%edges%logStats(alog)
    call self%faces%logStats(aLog)
    call EndSection(aLog)
end subroutine

!> @brief Writes a mesh to a Matlab-readable .m file
!>
!> @param[in] self Mesh to output
!> @param[in] fileunit Integer unit of .m file
subroutine writeMatlab(self, fileunit)
    class(PolyMesh2d), intent(in) :: self
    integer(kint), intent(in) :: fileunit

    call self%particles%writeMatlab(fileunit)
    call self%edges%writeMatlab(fileunit)
    call self%faces%writeMatlab(fileunit)
end subroutine

!@>brief Writes mesh data to a serial VTK .xml polydata file
!> The file is incomplete after this subroutine finishes, to allow
!> point data from @ref Field to be added.
!> Once all data is added, call writeVTKSerialEndXML to finish the file.
subroutine writeVTKSerialStartXML(self, fileunit)
    class(PolyMesh2d), intent(in) :: self
    integer(kint), intent(in) :: fileunit
    !
    integer(kint) :: i, j

    write(fileunit,'(A)') '<?xml version="1.0"?>'
    write(fileunit,'(A)') '<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian">'
    write(fileunit,'(A)') ' <PolyData>'
    write(fileunit,'(A,I8,A)',advance='no') '   <Piece NumberOfPoints="', self%particles%N, &
        '" NumberOfVerts="0" NumberOfLines="0" '
    if (self%faceKind == TRI_PANEL) then
        write(fileunit,'(A,I8,A)') 'NumberOfStrips="0" NumberOfPolys="', 3*self%faces%N_Active, '">'
    elseif(self%faceKind == QUAD_PANEL) then
        write(fileunit,'(A,I8,A)') 'NumberOfStrips="0" NumberOfPolys="', 4*self%faces%N_Active, '">'
    elseif(self%faceKind == QUAD_CUBIC_PANEL) then
        write(fileunit,'(A,I8,A)') 'NumberOfStrips="0" NumberOfPolys="', 9*self%faces%N_Active, '">'
    endif

    write(fileunit,'(A)') '    <Points>'
    if (self%geomKind == PLANAR_GEOM) then
        write(fileunit,'(A)') '      <DataArray type="Float32" NumberOfComponents="3" format="ascii">'
        do i=1, self%particles%N
            !write(fileunit,'(3(F24.9,A))',advance='no') self%particles%x(i), ' ', self%particles%y(i), ' ', dzero, ' '
            write(fileunit,'(3(F24.9,A))') self%particles%x(i), ' ', self%particles%y(i), ' ', dzero, ' '
        enddo
    elseif (self%geomKind == SPHERE_GEOM) then
        write(fileunit,'(A)') '      <DataArray type="Float32" NumberOfComponents="3" format="ascii">'
        do i=1, self%particles%N
            !write(fileunit,'(3(F24.9,A))',advance='no') self%particles%x(i), ' ', self%particles%y(i), ' ', &
            !    self%particles%z(i), ' '
            write(fileunit,'(3(F24.9,A))') self%particles%x(i), ' ', self%particles%y(i), ' ', &
                self%particles%z(i), ' '
        enddo
    endif
    !write(fileunit,'(A)') ''
    write(fileunit,'(A)') '      </DataArray>'
    write(fileunit,'(A)') '    </Points>'
    write(fileunit,'(A)') '    <Polys>'
    write(fileunit,'(A)') '      <DataArray type="Int32" Name="connectivity" format="ascii">'
    if (self%faceKind == TRI_PANEL) then
        do i=1, self%faces%N
            if (.not. self%faces%hasChildren(i)) then
                do j=1,3
                    write(fileunit,'(3(I20,A))') self%faces%vertices(j,i)-1, ' ', &
                        self%faces%vertices(mod(j,3)+1,i)-1, ' ', self%faces%centerParticles(1,i)-1, ' '
                enddo
            endif
        enddo
        !write(fileunit,'(A)') ''
        write(fileunit,'(A)') '      </DataArray>'
        write(fileunit,'(A)') '      <DataArray type="Int32" Name="offsets" format="ascii">'
        do i = 1, 3*self%faces%N_Active
            write(fileunit,'(I20)') 3*i
        enddo
    elseif (self%faceKind == QUAD_PANEL) then
        do i=1, self%faces%N
            if (.not. self%faces%hasChildren(i)) then
                do j=1,4
                    write(fileunit,'(4(I20,A))') self%faces%vertices(j,i)-1, ' ', &
                        self%faces%vertices(mod(j,4)+1,i)-1, ' ', self%faces%centerParticles(1,i)-1, ' '
                enddo
            endif
        enddo
        !write(fileunit,'(A)') ''
        write(fileunit,'(A)') '      </DataArray>'
        write(fileunit,'(A)') '      <DataArray type="Int32" Name="offsets" format="ascii">'
        do i = 1, 4*self%faces%N_Active
            write(fileunit,'(I20)') 3*i
        enddo
    elseif (self%faceKind == QUAD_CUBIC_PANEL) then
        do i=1,self%faces%N
            if (.not. self%faces%hasChildren(i)) then
                write(fileunit,'(4(I20,A))') self%faces%vertices(1,i)-1, ' ', self%faces%vertices(2,i)-1, ' ', &
                    self%faces%centerParticles(1,i)-1, ' ', self%faces%vertices(12,i)-1, ' '
                write(fileunit,'(4(I20,A))') self%faces%vertices(2,i)-1, ' ', self%faces%vertices(3,i)-1, ' ', &
                    self%faces%centerParticles(2,i)-1, ' ', self%faces%centerParticles(1,i)-1, ' '
                write(fileunit,'(4(I20,A))') self%faces%vertices(3,i)-1, ' ', self%faces%vertices(4,i)-1, ' ', &
                    self%faces%vertices(5,i)-1, ' ', self%faces%centerParticles(2,i)-1, ' '
                write(fileunit,'(4(I20,A))') self%faces%centerParticles(2,i)-1, ' ', &
                    self%faces%vertices(5,i)-1, ' ', self%faces%vertices(6,i)-1, ' ', self%faces%centerParticles(3,i)-1, ' '
                write(fileunit,'(4(I20,A))') self%faces%centerParticles(3,i)-1, ' ', &
                    self%faces%vertices(6,i)-1, ' ', self%faces%vertices(7,i)-1, ' ', self%faces%vertices(8,i)-1, ' '
                write(fileunit,'(4(I20,A))') self%faces%centerParticles(4,i)-1, ' ', &
                    self%faces%centerParticles(3,i)-1, ' ', self%faces%vertices(8,i)-1, ' ', self%faces%vertices(9,i)-1, ' '
                write(fileunit,'(4(I20,A))') self%faces%vertices(11,i)-1, ' ', &
                    self%faces%centerParticles(4,i)-1, ' ', self%faces%vertices(9,i)-1, ' ', self%faces%vertices(10,i)-1, ' '
                write(fileunit,'(4(I20,A))') self%faces%vertices(12,i)-1, ' ', &
                    self%faces%centerParticles(1,i)-1, ' ', self%faces%centerParticles(4,i)-1, ' ', self%faces%vertices(11,i)-1, ' '
                do j=1,4
                    write(fileunit,'(I20,A)') self%faces%centerParticles(j,i)-1, ' '
                enddo
            endif
        enddo
        !write(fileunit,'(A)') ''
        write(fileunit,'(A)') '      </DataArray>'
        write(fileunit,'(A)') '      <DataArray type="Int32" Name="offsets" format="ascii">'
        do i = 1, 9*self%faces%N_Active
            write(fileunit,'(I20,A)') 4*i
        enddo
    endif
    write(fileunit,'(A)') ''
    write(fileunit,'(A)') '      </DataArray>'
    write(fileunit,'(A)') '    </Polys>'
end subroutine

subroutine writeLagCoordsToVTKPointDataXML(self, fileunit)
    class(PolyMesh2d), intent(in) :: self
    integer(kint), intent(in) :: fileunit
    !
    integer(kint) :: i,j
    
    write(fileunit,'(A)') '      <DataArray type="Float32" Name="LagCoords" NumberOfComponents="3" format="ascii">'
    if (self%geomKind == PLANAR_GEOM) then
        do i=1,self%particles%N
            write(fileunit,'(3(G15.8,A))') self%particles%x0(i), " ", self%particles%y0(i), " ", dzero, ""
        enddo    
    else
        do i=1,self%particles%N
            write(fileunit,'(3(G15.8,A))') self%particles%x0(i), " ", self%particles%y0(i), " ", self%particles%z0(i), ""
        enddo
    endif
    write(fileunit,'(A)') '     </DataArray>'
end subroutine

subroutine writePlanarVTKStartWithHeightFromField(self, hfield, fileunit)
    class(PolyMesh2d), intent(in) :: self
    class(Field), intent(in) :: hField
    integer(kint), intent(in) :: fileunit
    !
    integer(kint) :: i, j

    if (self%geomKind /= PLANAR_GEOM) then
        call LogMessage(log, ERROR_LOGGING_LEVEL, trim(logKey)//" writePlanarVTKStartWithHeightFromField error: ", &
            "cannot use height field on 3d data.")
        return
    endif
    if (hfield%nDim /= 1) then
        call LogMessage(log, ERROR_LOGGING_LEVEL, trim(logKey)//" writePlanarVTKStartWithHeightFromField error: ", &
            "only scalar fields can be used for height")
        return
    endif
    if (hfield%N /= self%particles%N) then
        write(logstring,'(2(A,I8))') "particles%N = ", self%particles%N, ", hField%n = ", hField%n
        call LogMessage(log, ERROR_LOGGING_LEVEL, trim(logKey)//&
            " writePlanarVTKStartWithHeightFromField error: field size mismatch", logstring)
        return
    endif

    write(fileunit,'(A)') '<?xml version="1.0"?>'
    write(fileunit,'(A)') '<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian">'
    write(fileunit,'(A)') ' <PolyData>'
    write(fileunit,'(A,I8,A)',advance='no') '   <Piece NumberOfPoints="', self%particles%N, &
        '" NumberOfVerts="0" NumberOfLines="0" '
    if (self%faceKind == TRI_PANEL) then
        write(fileunit,'(A,I8,A)') 'NumberOfStrips="0" NumberOfPolys="', 3*self%faces%N_Active, '">'
    elseif(self%faceKind == QUAD_PANEL) then
        write(fileunit,'(A,I8,A)') 'NumberOfStrips="0" NumberOfPolys="', 4*self%faces%N_Active, '">'
    elseif(self%faceKind == QUAD_CUBIC_PANEL) then
        write(fileunit,'(A,I8,A)') 'NumberOfStrips="0" NumberOfPolys="', 9*self%faces%N_Active, '">'
    endif
    write(fileunit,'(A)') '    <Points>'

    write(fileunit,'(A)') '      <DataArray type="Float32" NumberOfComponents="3" format="ascii">'
    do i=1, self%particles%N
        write(fileunit,'(3(F24.9,A))') self%particles%x(i), ' ', self%particles%y(i), ' ', hfield%comp1(i), ' '
    enddo
    write(fileunit,'(A)') '      </DataArray>'
    write(fileunit,'(A)') '    </Points>'
    write(fileunit,'(A)') '    <Polys>'
    write(fileunit,'(A)') '      <DataArray type="Int32" Name="connectivity" format="ascii">'
    if (self%faceKind == TRI_PANEL) then
        do i=1, self%faces%N
            if (.not. self%faces%hasChildren(i)) then
                do j=1,3
                    write(fileunit,'(3(I20,A))') self%faces%vertices(j,i)-1, ' ', &
                        self%faces%vertices(mod(j,3)+1,i)-1, ' ', self%faces%centerParticles(1,i)-1, ' '
                enddo
            endif
        enddo
        !write(fileunit,'(A)') ''
        write(fileunit,'(A)') '      </DataArray>'
        write(fileunit,'(A)') '      <DataArray type="Int32" Name="offsets" format="ascii">'
        do i = 1, 3*self%faces%N_Active
            write(fileunit,'(I20)') 3*i
        enddo
    elseif (self%faceKind == QUAD_PANEL) then
        do i=1, self%faces%N
            if (.not. self%faces%hasChildren(i)) then
                do j=1,4
                    write(fileunit,'(4(I20,A))') self%faces%vertices(j,i)-1, ' ', &
                        self%faces%vertices(mod(j,4)+1,i)-1, ' ', self%faces%centerParticles(1,i)-1, ' '
                enddo
            endif
        enddo
        !write(fileunit,'(A)') ''
        write(fileunit,'(A)') '      </DataArray>'
        write(fileunit,'(A)') '      <DataArray type="Int32" Name="offsets" format="ascii">'
        do i = 1, 4*self%faces%N_Active
            write(fileunit,'(I20)') 3*i
        enddo
    elseif (self%faceKind == QUAD_CUBIC_PANEL) then
        do i=1,self%faces%N
            if (.not. self%faces%hasChildren(i)) then
                write(fileunit,'(4(I20,A))') self%faces%vertices(1,i)-1, ' ', self%faces%vertices(2,i)-1, ' ', &
                    self%faces%centerParticles(1,i)-1, ' ', self%faces%vertices(12,i)-1, ' '
                write(fileunit,'(4(I20,A))') self%faces%vertices(2,i)-1, ' ', self%faces%vertices(3,i)-1, ' ', &
                    self%faces%centerParticles(2,i)-1, ' ', self%faces%centerParticles(1,i)-1, ' '
                write(fileunit,'(4(I20,A))') self%faces%vertices(3,i)-1, ' ', self%faces%vertices(4,i)-1, ' ', &
                    self%faces%vertices(5,i)-1, ' ', self%faces%centerParticles(2,i)-1, ' '
                write(fileunit,'(4(I20,A))') self%faces%centerParticles(2,i)-1, ' ', &
                    self%faces%vertices(5,i)-1, ' ', self%faces%vertices(6,i)-1, ' ', self%faces%centerParticles(3,i)-1, ' '
                write(fileunit,'(4(I20,A))') self%faces%centerParticles(3,i)-1, ' ', &
                    self%faces%vertices(6,i)-1, ' ', self%faces%vertices(7,i)-1, ' ', self%faces%vertices(8,i)-1, ' '
                write(fileunit,'(4(I20,A))') self%faces%centerParticles(4,i)-1, ' ', &
                    self%faces%centerParticles(3,i)-1, ' ', self%faces%vertices(8,i)-1, ' ', self%faces%vertices(9,i)-1, ' '
                write(fileunit,'(4(I20,A))') self%faces%vertices(11,i)-1, ' ', &
                    self%faces%centerParticles(4,i)-1, ' ', self%faces%vertices(9,i)-1, ' ', self%faces%vertices(10,i)-1, ' '
                write(fileunit,'(4(I20,A))') self%faces%vertices(12,i)-1, ' ', &
                    self%faces%centerParticles(1,i)-1, ' ', self%faces%centerParticles(4,i)-1, ' ', self%faces%vertices(11,i)-1, ' '
                do j=1,4
                    write(fileunit,'(I20,A)') self%faces%centerParticles(j,i)-1, ' '
                enddo
            endif
        enddo
        !write(fileunit,'(A)') ''
        write(fileunit,'(A)') '      </DataArray>'
        write(fileunit,'(A)') '      <DataArray type="Int32" Name="offsets" format="ascii">'
        do i = 1, 9*self%faces%N_Active
            write(fileunit,'(I20,A)') 4*i
        enddo
    endif
    write(fileunit,'(A)') ''
    write(fileunit,'(A)') '      </DataArray>'
    write(fileunit,'(A)') '    </Polys>'
end subroutine

!> @brief Completes the .xml file
subroutine writeVTKSerialEndXML(self, fileunit)
    class(PolyMesh2d), intent(in) :: self
    integer(kint), intent(in) :: fileunit
    write(fileunit,'(A)') '    </Piece>'
    write(fileunit,'(A)') '  </PolyData>'
    write(fileunit,'(A)') '</VTKFile>'
end subroutine

!> @brief Computes the number of vertices in a uniformly refined mesh.
!> Used to compute memory allocation requirements.
!> @param[in] self Target mesh
!> @param[in] initNest initial refinement level
!> @return Number of vertices in a mesh
pure function nVerticesInMesh(self, nestLevel)
    integer(kint) :: nVerticesInMesh
    class(PolyMesh2d), intent(in) :: self
    integer(kint), intent(in) :: nestLevel
    !
    integer(kint) :: i
    nVerticesInMesh = 0
    select case (self%meshSeed)
        case (TRI_HEX_SEED)
            do i = 2**nestLevel +1, 2**(nestLevel+1)
                nVerticesInMesh = nVerticesInMesh + i
            enddo
            nVerticesInMesh = 2*nVerticesInMesh + 2**(nestLevel+1) + 1
        case (QUAD_RECT_SEED)
            nVerticesInMesh = 3
            do i=1,nestLevel
                nVerticesInMesh = nVerticesInMesh + 2**i
            enddo
            nVerticesInMesh = nVerticesInMesh * nVerticesInMesh
        case (CUBIC_PLANE_SEED)
            nVerticesInMesh = 49
            do i=1,nestLevel
                nVerticesInMesh = (2*sqrt(real(nVerticesInMesh,kreal))-1)**2
            enddo
        case (ICOS_TRI_SPHERE_SEED)
            nVerticesInMesh = 2 + 10*4**nestLevel
        case (CUBED_SPHERE_SEED)
            nVerticesInMesh = 2 + 6*4**nestLevel
    end select
end function

!> @brief Computes the number of faces in a uniformly refined mesh.
!> Used to compute memory allocation requirements.
!> @param[in] self Target mesh
!> @param[in] initNest initial refinement level
!> @return Number of faces in a mesh
pure function nFacesInMesh(self, nestLevel)
    integer(kint) :: nFacesInMesh
    class(PolyMesh2d),intent(in) :: self
    integer(kint), intent(in) :: nestLevel
    !
    integer(kint) :: i

    nFacesInMesh = 0
    select case (self%meshSeed)
        case (TRI_HEX_SEED)
            nFacesInMesh = 6*4**nestLevel
        case (QUAD_RECT_SEED)
            nFacesInMesh = 4*4**nestLevel
        case (CUBIC_PLANE_SEED)
            nFacesInMesh = 4*4**nestLevel
        case (ICOS_TRI_SPHERE_SEED)
            nFacesInMesh = 20*4**nestLevel
        case (CUBED_SPHERE_SEED)
            nFacesInMesh = 6*4**nestLevel
    end select
end function

!> @brief Uses the number of vertices and the number of faces (via Euler's formulas) to determine the the number of edges in a uniformly refined mesh.
!> Used to compute memory allocation requirements.
!> @param[in] self Target mesh
!> @param[in] nVertices number of vertices in a mesh
!> @param[in] nFaces number of faces in a mesh
!> @return Number of edges in a mesh
pure function nEdgesInMesh(self, nVerts, nFaces, nestLevel)
    integer(kint) :: nEdgesInMesh
    class(PolyMesh2d), intent(in) :: self
    integer(kint), intent(in) :: nVerts, nFaces
    integer(kint), intent(in), optional :: nestLevel
    !
    integer(kint) :: nv, nf, i

    nEdgesInMesh = 0
    select case (self%meshSeed)
        case (QUAD_RECT_SEED)
            nEdgesInMesh = nFaces + nVerts -1
        case (TRI_HEX_SEED)
            nEdgesInMesh = nFaces + nVerts -1
        case (CUBED_SPHERE_SEED)
            nEdgesInMesh = nFaces + nVerts -2
        case (ICOS_TRI_SPHERE_SEED)
            nEdgesInMesh = nFaces + nVerts -2
        case (CUBIC_PLANE_SEED)
            nf = 4*4**nestLevel
            nv = 3
            do i=1,nestLevel
                nv = nv + 2**i
            enddo
            nv = nv*nv
            nEdgesInMesh = nf + nv - 1
    end select
end function

function ccwAdjacentFaces(self, faceIndex)
    type(STDIntVector) :: ccwAdjacentFaces
    class(PolyMesh2d), intent(in) :: self
    integer(kint), intent(in) :: faceIndex
    !
    type(STDIntVector) :: leafEdges
    integer(kint) :: i

    leafEdges = self%ccwEdgesAroundFace(faceIndex)
    call initialize(ccwAdjacentFaces)
    do i=1,leafEdges%n
        if (self%edges%positiveOrientation(leafEdges%int(i), faceIndex)) then
            call ccwAdjacentFaces%pushBack(self%edges%rightFace(leafEdges%int(i)))
        else
            call ccwAdjacentFaces%pushBack(self%edges%leftFace(leafEdges%int(i)))
        endif
    enddo
end function

function ccwEdgesAroundFace(self, faceIndex)
    type(STDIntVector) :: ccwEdgesAroundFace
    class(PolyMesh2d), intent(in) :: self
    integer(kint), intent(in) :: faceIndex
    !
    integer(kint) :: i, j, ne
    type(STDIntVector), dimension(4) :: edgeLeaves

    if (self%faceKind == TRI_PANEL) then
        ne = 3
    else
        ne = 4
    endif
    do i=1, ne
        edgeLeaves(i) = self%edges%getLeafEdges(self%faces%edges(i,faceIndex))
    enddo

    call initialize(ccwEdgesAroundFace)
    do i=1,ne
        do j=1, edgeLeaves(i)%n
            call ccwEdgesAroundFace%pushBack(edgeLeaves(i)%int(j))
        enddo
    enddo
end function

function locatePointInFace(self, queryPt)
    integer(kint) :: locatePointinFace
    class(PolyMesh2d), intent(in) :: self
    real(kreal), dimension(3), intent(in) :: queryPt
    !
    integer(kint) :: treeStart, walkStart

    locatePointInFace = 0
    treeStart = self%nearestRootFace(queryPt)
    walkStart = self%locatePointTree(queryPt, treeStart)
    locatePointInFace = self%locatePointWalk(queryPt, walkStart)
end function

function nearestRootFace(self, queryPt)
    integer(kint) :: nearestRootFace
    real(kreal), dimension(3), intent(in) :: queryPt
    class(PolyMesh2d), intent(in) :: self
    !
    integer(kint) :: i
    real(kreal) :: dist, testDist
    real(kreal), dimension(3) :: centroid

    nearestRootFace = 1
    centroid = self%faces%physCentroid(1, self%particles)
    if (self%geomKind == SPHERE_GEOM) then
        dist = SphereDistance(centroid, queryPt)
        do i=2, self%nRootFaces
            centroid = self%faces%physcentroid(i, self%particles)
            testDist = SphereDistance(centroid, queryPt)
            if (testDist < dist) then
                dist = testDist
                nearestRootFace = i
            endif
        enddo
    else
        dist = sqrt(sum((queryPt-centroid)*(queryPt-centroid)))
        do i=2, self%nRootFaces
            centroid = self%faces%physCentroid(i, self%particles)
            testDist = sqrt(sum((queryPt-centroid)*(queryPt-centroid)))
            if (testDist < dist) then
                dist = testDist
                nearestRootFace = i
            endif
        enddo
    endif
end function

recursive function locatePointWalk(self, queryPt, startIndex) result(faceIndex)
    integer(kint) :: faceIndex
    class(PolyMesh2d), intent(in) :: self
    real(kreal), dimension(3), intent(in) :: queryPt
    integer(kint), intent(in) :: startIndex
    !
    real(kreal) :: dist, testDist
    real(kreal), dimension(3) :: centroid
    integer(kint) :: currentFace, i
    type(STDIntVector) :: adjFaces

    if (startIndex == 0) then
        call LogMessage(log, ERROR_LOGGING_LEVEL, trim(logkey)//" locatePointWalkSearch ERROR:"," invalid startIndex.")
		return
    endif
    if (self%faces%hasChildren(startIndex)) then
        call LogMessage(log, ERROR_LOGGING_LEVEL, trim(logkey)//" locatePointWalkSearch ERROR:"," expected a leaf face.")
		return
    endif

    adjFaces = ccwAdjacentFaces(self, startIndex)

    faceIndex = 0
    currentFace = startIndex
    centroid = self%faces%physCentroid(startIndex, self%particles)
    if (self%geomKind == SPHERE_GEOM) then
        dist = SphereDistance(centroid, queryPt)
    else
        dist = sqrt(sum((centroid-queryPt)*(centroid-queryPt)))
    endif
    do i=1, adjFaces%n
        if (adjFaces%int(i) > 0 ) then ! skip faces adjacent to domain boundary
            centroid = self%faces%physCentroid(adjFaces%int(i), self%particles)
            if (self%geomKind == SPHERE_GEOM) then
                testDist = SphereDistance(centroid, queryPt)
            else
                testDist = sqrt(sum((centroid-queryPt)*(centroid-queryPt)))
            endif
            if (testDist < dist) then
                dist = testDist
                currentFace = adjFaces%int(i)
            endif
        endif
    enddo

    if (currentFace == startIndex) then
        faceIndex = currentFace
    else
        faceIndex = locatePointWalk(self, queryPt, currentFace)
    endif
end function

recursive function locatePointTree(self, queryPt, startIndex) result(faceIndex)
    integer(kint) :: faceIndex
    class(PolyMesh2d), intent(in) :: self
    real(kreal), dimension(3), intent(in) :: queryPt
    integer(kint), intent(in) :: startIndex
    !
    real(kreal) :: dist, testDist
    real(kreal), dimension(3) :: centroid
    integer(kint) :: i, nearestChild

    faceIndex = 0
    nearestChild = 0
    centroid = dzero
    dist = 999.0d20
    if (self%faces%hasChildren(startIndex)) then
        do i=1,4
            centroid = self%faces%physCentroid(self%faces%children(i,startIndex), self%particles)
            if (self%geomKind == SPHERE_GEOM) then
                testDist = SphereDistance(queryPt, centroid)
            else
                testDist = sqrt(sum( (queryPt-centroid)*(queryPt-centroid)))
            endif
            if (testDist < dist) then
                nearestChild = self%faces%children(i, startIndex)
                dist = testDist
            endif
        enddo
        faceIndex = locatePointTree(self, queryPt, nearestChild)
    else
        faceIndex = startIndex
    endif
end function

!> @brief Initializes a mesh using a mesh seed from a file to define the root set of @ref Particles, @ref Edges, and @ref Faces.
!> Mesh seeds are defined in both @ref NumberKinds and via the corresponding.  See meshSeeds.py for complete details.
!> Multiplies particle postions by a specified amount to change the spatial domain (all mesh seeds approximately cover
!> the unit sphere of their embedded space).
!>
!> @param[inout] self Target mesh to be initialized
subroutine getSeed(self)
    class(PolyMesh2d), intent(inout) :: self
    !
    character(len=128) :: seedfilename
    integer(kint) :: i, j, nSeedParticles, nSeedEdges, nSeedFaces, nSeedVerts, nv, nc
    integer(kint), allocatable :: seedEdgeOrigs(:), seedEdgeDests(:), seedEdgeLefts(:), seedEdgeRights(:)
    integer(kint), allocatable :: seedEdgeInts(:,:), seedFaceEdges(:,:), seedFaceVerts(:,:), seedFaceCenters(:,:)
    real(kreal), allocatable :: seedXYZ(:,:)
    real(kreal) :: jac(16)
    integer(kint) :: nMaxParticles, nMaxEdges, nMaxFaces

!    call LogMessage(log, DEBUG_LOGGING_LEVEL, trim(logkey)//" entering : ", "getSeed")
!    if (associated(self%particles)) deallocate(self%particles)
!    if (associated(self%edges)) deallocate(self%edges)
!    if (associated(self%faces)) deallocate(self%faces)

    allocate( Particles :: self%particles)

    if (self%meshSeed == TRI_HEX_SEED) then
         allocate( LinearEdges :: self%edges)
         allocate( TriLinearFaces :: self%faces)

        seedfilename = MESH_SEED_DIR//"/triHex.namelist"

        nSeedParticles = 13
        nSeedEdges = 12
        nSeedFaces = 6
        nSeedVerts = 7

        self%faceKind = TRI_PANEL
        self%geomKind = PLANAR_GEOM
        nv = 3
        nc=1
    elseif (self%meshSeed == QUAD_RECT_SEED) then
        allocate(LinearEdges :: self%edges)
        allocate(QuadLinearFaces :: self%faces)

        seedfilename = MESH_SEED_DIR//"/quadRect.namelist"
        nSeedParticles = 13
        nSeedEdges = 12
        nSeedFaces = 4
        nv = 4
        nc=1
        self%faceKind = QUAD_PANEL
        self%geomKind = PLANAR_GEOM
    elseif (self%meshSeed == CUBIC_PLANE_SEED) then
        allocate(CubicEdges :: self%edges)
        allocate(QuadCubicFaces :: self%faces)
        seedfilename = MESH_SEED_DIR//"/quadCubic.namelist"
        nSeedParticles = 49
        nSeedEdges = 12
        nSeedFaces = 4
        nv=12
        nc=4
        self%faceKind= QUAD_CUBIC_PANEL
        self%geomKind = PLANAR_GEOM
    elseif (self%meshSeed == ICOS_TRI_SPHERE_SEED) then
        allocate(LinearEdges :: self%edges)
        allocate(TriLinearFaces :: self%faces)

        seedfilename = MESH_SEED_DIR//"/icosTri.namelist"
        nSeedParticles = 32
        nSeedEdges = 30
        nSeedFaces = 20
        nv=3
        nc=1
        self%faceKind = TRI_PANEL
        self%geomKind = SPHERE_GEOM
    elseif (self%meshSeed == CUBED_SPHERE_SEED) then
        allocate(LinearEdges :: self%edges)
        allocate(QuadLinearFaces :: self%faces)

        seedfilename = MESH_SEED_DIR//"/cubedSphere.namelist"
        nSeedParticles = 14
        nSeedEdges = 12
        nSeedFaces = 6
        nv = 4
        nc=1
        self%faceKind = QUAD_PANEL
        self%geomKind = SPHERE_GEOM
!    elseif (self%mesh_type == "cubedSphereCubic") then
!        allocate(CubicEdges :: self%edges)
!        allocate(QuadCubicFaces :: self%faces)
!
!        self%faceKind = QUAD_CUBIC_PANEL
!        self%geomKind = SPHERE_GEOM
!        nv = 12
!        nc = 4
    endif

    if (self%meshSeed == CUBIC_PLANE_SEED) then
        nMaxParticles = self%nVerticesInMesh(self%maxNest)
    else
        nMaxParticles = self%nVerticesInMesh(self%maxNest) + self%nFacesInMesh(self%maxNest)
    endif
    nMaxFaces = 0
    nMaxEdges = 0
    do i=0, self%maxNest
        nMaxFaces = nMaxFaces + self%nFacesInMesh(i)
        nMaxEdges = nMaxEdges + self%nEdgesInMesh(self%nVerticesInMesh(i), self%nFacesInMesh(i),i)
    enddo

    self%nRootFaces = nSeedFaces
    call self%particles%init(nMaxParticles, self%geomKind)
    call self%edges%init(nMaxEdges)
    call self%faces%init(self%faceKind, nMaxFaces)


    allocate(seedXYZ(3,nSeedParticles))
    allocate(seedEdgeOrigs(nSeedEdges))
    allocate(seedEdgeDests(nSeedEdges))
    allocate(seedEdgeLefts(nSeedEdges))
    allocate(seedEdgeRights(nSeedEdges))
    if (self%faceKind == TRI_PANEL) then
        allocate(seedFaceVerts(3,nSeedFaces))
        allocate(seedFaceEdges(3,nSeedFaces))
        allocate(seedFaceCenters(1,nSeedFaces))
        allocate(seedEdgeInts(1,1))
    elseif (self%faceKind == QUAD_PANEL) then
        allocate(seedFaceVerts(4,nSeedFaces))
        allocate(seedFaceEdges(4,nSeedFaces))
        allocate(seedFaceCenters(1,nSeedFaces))
        allocate(seedEdgeInts(1,1))
    elseif (self%faceKind == QUAD_CUBIC_PANEL) then
        allocate(seedFaceVerts(12,nSeedFaces))
        allocate(seedFaceEdges(4,nSeedFaces))
        allocate(seedFaceCenters(4,nSeedFaces))
        allocate(seedEdgeInts(2,nSeedEdges))
    endif

    call self%readSeedFile(seedfilename, seedXYZ, seedEdgeOrigs, seedEdgeDests, seedEdgeLefts, seedEdgeRights, &
        seedEdgeInts, seedFaceVerts, seedFaceEdges, seedFaceCenters)

    seedXYZ = self%maxRadius * seedXYZ

    ! initialize mesh
    do i=1, nSeedParticles
        call self%particles%insert(seedXYZ(:,i), seedXYZ(:,i))
        !print *, "new particles ", i, " = (", self%particles%physCoord(i), ")"
    enddo

    if ( self%particles%N /= nSeedParticles ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL, logkey//" initMeshFromSeed ERROR : "," particles%N.")
	endif

	do i=1, nSeedEdges
	    call self%edges%insert(seedEdgeOrigs(i), seedEdgeDests(i), seedEdgeLefts(i), seedEdgeRights(i), seedEdgeInts(:,i))
	enddo

	if ( self%edges%N /= nSeedEdges ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL, logkey//" initMeshFromSeed ERROR : "," edges%N.")
	endif

	do i=1, nSeedFaces
	    call self%faces%insert(seedFaceCenters(:,i), seedFaceVerts(:,i), seedFaceEdges(:,i))
	enddo
	!
	!   set face areas and particle weights
	!
	self%faces%N_Active = nSeedFaces
	self%faces%area = dzero
    do i=1, nSeedFaces
        self%faces%area(i) = self%faces%setArea(i, self%particles)
    enddo


    if (self%faceKind == QUAD_PANEL .or. self%faceKind == TRI_PANEL) then
        do i=1, nSeedFaces
            self%particles%weight(self%faces%centerParticles(1,i)) = self%faces%area(i)
        enddo
    elseif (self%faceKind == QUAD_CUBIC_PANEL) then
        do i=1, nSeedFaces
            seedXYZ(:,1) = self%particles%physCoord(self%faces%vertices(1,i))
            seedXYZ(:,2) = self%particles%physCoord(self%faces%vertices(4,i))
            seedXYZ(:,3) = self%particles%physCoord(self%faces%vertices(7,i))
            seedXYZ(:,4) = self%particles%physCoord(self%faces%vertices(10,i))
            do j=1,12
                jac(j) = bilinearPlaneJacobian(seedXYZ(:,1:4), quad16_vertex_qp(1,j), quad16_vertex_qp(2,j))
            enddo
            do j=1,4
                jac(12+j) = bilinearPlaneJacobian(seedXYZ(:,1:4), quad16_center_qp(1,j), quad16_center_qp(2,j))
            enddo
            do j=1,12
                self%particles%weight(self%faces%vertices(j,i)) = quad16_vertex_qw(j) * jac(j)
            enddo
            do j=1,4
                self%particles%weight(self%faces%centerParticles(j,i)) = quad16_center_qw(j)*jac(12+j)
            enddo
        enddo
    endif
end subroutine

!> @brief Reads data from a mesh seed file.
!>
!> @param[in] self Target mesh
!> @param[in] seedFilename filename
!> @param[out] xyz position vectors of mesh seed particles
!> @param[out] origs indices of origin particles for each edge
!> @param[out] dests indices of destination particles for each edge
!> @param[out] lefts indices of left faces for each edge
!> @param[out] rights indices of right faces for each edge
!> @param[out] faceVerts indices of vertex particles for each face
!> @param[out] faceEdges indices of edges for each face
!> @param[out] faceCenters indices of faces' internal particles
subroutine readSeedFile(self, seedfilename, xyz, origs, dests, lefts, rights, ints, faceverts, faceedges, facecenters)
    class(PolyMesh2d), intent(in) :: self
    character(len=*), intent(in) :: seedfilename
    real(kreal), dimension(:,:), intent(inout) :: xyz
    integer(kint), dimension(:), intent(inout) :: origs, dests, lefts, rights
    integer(kint), dimension(:,:), intent(inout) :: ints
    integer(kint), dimension(:,:), intent(inout) :: faceverts, faceedges, facecenters
    !
    integer(kint) :: readStat
    integer(kint) :: nParticles, nEdges, nFaces, nVertsPerFace, nCentersPerFace, nIntsPerEdge, i
    namelist /seed/ xyz, origs, dests, lefts, rights, ints, faceverts, faceedges, facecenters

!    call LogMessage(log, DEBUG_LOGGING_LEVEL,trim(logkey)//" entering : ", "readSeedFile")
!    print *, "nParticles = ", size(xyz,1)
!    print *, "size(origs) = ", size(origs)
!    print *, "size(dests) = ", size(dests)
!    print *, "size(faceVerts,1)", size(faceverts,1)
!    print *, "size(faceverts,2)", size(faceverts,2)
!    print *, "size(ints)", size(ints,1), size(ints,2)

    xyz = dzero
    origs = 0
    dests = 0
    lefts = 0
    rights = 0
    ints = 0
    faceverts = 0
    faceedges = 0
    facecenters = 0
    open(unit=READ_UNIT, file=seedfilename, action='read', status='old', iostat=readStat)
    if (readStat /= 0) then
        call LogMessage(log, ERROR_LOGGING_LEVEL, "readSeedFile error : cannot open file ", seedfilename)
        return
    endif
    read(READ_UNIT, nml=seed, iostat=readstat)
    if (readStat /= 0) then
        call LogMessage(log, ERROR_LOGGING_LEVEL, "readSeedFile error : cannot read namelist ", seedfilename)
        return
    endif    
    close(READ_UNIT)
!    print *, "readstat ", readstat
!    print *, "Particles in file : "
!    do i=1, nParticles
!        print *, i, ": (", xyz(:,i), ")"
!    enddo
    
    origs = origs+1
    dests = dests+1
    lefts = lefts+1
    rights= rights+1
    ints = ints+1
    faceverts = faceverts+1
    faceedges = faceedges+1
    facecenters = facecenters+1
end subroutine

!> @brief Deletes and frees memory associated with a @ref PolyMesh2d object.
!> @param[inout] self Target mesh
subroutine deleteMesh(self)
    type(PolyMesh2d), intent(inout) :: self
    deallocate(self%particles)
    nullify(self%particles)
    deallocate(self%edges)
    nullify(self%edges)
    deallocate(self%faces)
    nullify(self%faces)
end subroutine

!> @brief Initializes a logger for the Faces module
!>
!> Output is controlled both by message priority and by MPI Rank
!> @param[out] aLog Target Logger object
!> @param[in] rank Rank of this processor
subroutine InitLogger(aLog,rank)
! Initialize a logger for this module and processor
	type(Logger), intent(out) :: aLog
	integer(kint), intent(in) :: rank
	write(logKey,'(A,A,I0.2,A)') trim(logKey),'_',rank,' : '
	if ( rank == 0 ) then
		call New(aLog,logLevel)
	else
		call New(aLog,ERROR_LOGGING_LEVEL)
	endif
	logInit = .TRUE.
end subroutine

!> @}
end module

