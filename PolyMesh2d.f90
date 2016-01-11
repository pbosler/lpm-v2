module PolyMesh2dModule
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
use STDIntVectorModule
use OutputWriterModule
use LoggerModule
use ParticlesModule
use EdgesModule
use FacesModule
use PlaneGeomModule
use SphereGeomModule

implicit none

private
public PolyMesh2d
public New, Delete, Copy
public locateFaceContainingPoint, nearestParticle
public pointIsOutsideMesh
public CCWEdgesAroundFace, CCWVerticesAroundFace, CCWAdjacentFaces
public CCWFacesAroundVertex
public LogStats, PrintDebugInfo
public WriteMeshToMatlab, WriteMeshToVTKPolyData
public MeshSeedString
public ProjectParticlesToSphere
public MaxLagrangianVariationPerFace

!> @brief Combines mesh primitives into one complete mesh object.
type PolyMesh2d
	type(Particles) :: particles !< @ref Particles defining the centers and vertices of each mesh face
	type(Edges) :: edges !< @ref Edges connecting origin and destination particles for each edge, and each edge's left and right face
	type(Faces) :: faces !< @ref Faces 
	integer(kint) :: faceKind = 0 !< face kind (e.g., triangular or quadrilateral) as defined in @ref NumberKinds
	integer(kint) :: geomKind = 0 !< geometry kind (e.g., planar or spherical) as defined in @ref NumberKinds
	integer(kint) :: meshSeed = 0 !< mesh seed as defined in @ref NumberKinds
	integer(kint) :: initNest = -1 !< initial refinement level (in terms of the @ref Faces quadtree) for a mesh
	integer(kint) :: amrLimit = 0 !< Number of times a face may be divided beyond the initial refinement level
	real(kreal) :: t = 0.0_kreal !< time value represented by the mesh
	
	contains
		final :: deletePrivate
end type

interface New
	module procedure newPrivate
end interface

interface Delete
	module procedure deletePrivate
end interface

interface Copy
	module procedure copyPrivate
end interface

interface LogStats
	module procedure logStatsPrivate
end interface

interface PrintDebugInfo
	module procedure printDebugPrivate
end interface

!
!----------------
! Logging
!----------------
!
logical(klog), save :: logInit = .FALSE.
type(Logger) :: log
character(len=28), save :: logKey = 'PolyMesh2d'
integer(kint), parameter :: logLevel = TRACE_LOGGING_LEVEL
character(len=MAX_STRING_LENGTH) :: logstring

contains

!> @brief Allocates memory for a new mesh (@ref Particles, @ref Edges, and @ref Faces).
!> Initializes a new mesh using a "mesh seed," which is recursively refined until a desired spatial resolution is achieved.
!> The recursive refinement defines the @ref Faces quadtree and the @ref Edges binary tree.
subroutine newPrivate(self, meshSeed, initNest, maxNest, amrLimit, ampFactor )
	type(PolyMesh2d), intent(out) :: self
	integer(kint), intent(in) :: meshSeed
	integer(kint), intent(in) :: maxNest
	integer(kint), intent(in) :: initNest
	integer(kint), intent(in) :: amrLimit
	real(kreal), intent(in) :: ampFactor
	!
	integer(kint) :: nMaxParticles, nMaxVertices, nMaxFaces, nMaxEdges
	integer(kint) :: i, j, startIndex, nFacesOld

	if ( .NOT. logInit) call InitLogger(log, procRank)
	
	self%meshSeed = meshSeed
	self%initNest = initNest
	self%amrLimit = amrLimit
	if ( meshSeed < TRI_HEX_SEED .OR. meshSeed > BETA_PLANE_SEED ) then
		call LogMessage(log,ERROR_LOGGING_LEVEL, logkey," NewMesh ERROR : invalid meshSeed.")
		return
	endif
	self%faceKind = FaceKindFromSeed(meshSeed)
	self%geomKind = GeomKindFromSeed(meshSeed)
	
	nMaxParticles = nVerticesInMesh(self, maxNest) + nFacesInMesh(self, maxNest)
	nMaxVertices = 0
	nMaxFaces = 0
	nMaxEdges = 0
	do i = 0, maxNest
		nMaxFaces = nMaxFaces + nFacesInMesh(self, i)
		nMaxEdges = nMaxEdges + nEdgesInMesh(self, nVerticesInMesh(self,i), nFacesInMesh(self,i) )
	enddo
	
	call LogMessage(log, DEBUG_LOGGING_LEVEL,trim(logkey)//" allocating memory for nParticles = ", nMaxParticles )
	call LogMessage(log, DEBUG_LOGGING_LEVEL,trim(logkey)//" allocating memory for nEdges = ", nMaxEdges )
	call LogMessage(log, DEBUG_LOGGING_LEVEL,trim(logkey)//" allocating memory for nFaces = ", nMaxFaces )
	call LogMessage(log, DEBUG_LOGGING_LEVEL,trim(logkey)//" faceKind = ", self%faceKind )
	call LogMessage(log, DEBUG_LOGGING_LEVEL,trim(logkey)//" geomKind = ", self%geomKind )
	
	call New(self%particles, nMaxParticles, self%geomKind )
	call New(self%edges, nMaxEdges)
	call New(self%faces, self%faceKind, nMaxFaces )
	
	call initializeMeshFromSeed(self, ampFactor)
	
	startIndex = 1
	do i = 1, initNest
		nFacesOld = self%faces%N
		do j = startIndex, nFacesOld
			if ( self%faceKind == TRI_PANEL ) then
				call DivideTriFace(self%faces, j, self%particles, self%edges)
			elseif ( self%faceKind == QUAD_PANEL ) then
				call DivideQuadFace(self%faces,j, self%particles, self%edges)
			endif
		enddo
		startIndex = nFacesOld + 1
	enddo
	
	do i = 1, self%particles%N
		if ( self%particles%isPassive(i) ) call SortIncidentEdgesAtParticle( self%particles, i )
	enddo
end subroutine

!> @brief Deletes and frees memory associated with a @ref PolyMesh2d object.
!> @param[inout] self Target mesh
subroutine deletePrivate(self)
	type(PolyMesh2d), intent(inout) :: self
	call Delete(self%faces)
	call Delete(self%edges)
	call Delete(self%particles)
!	self%meshSeed = 0
!	self%initNest = 0
!	self%amrLimit = 0
!	self%geomKind = 0
!	self%faceKind = 0
end subroutine

!> @brief Peforms a deep copy of one mesh to another (all mesh components are deep copied).
!> Both meshes must be allocated appropriately prior to calling this subroutine.
!> @param[inout] self Target mesh
!> @param[in] other Source mesh
subroutine copyPrivate( self, other )
	type(PolyMesh2d), intent(inout) :: self
	type(PolyMesh2d), intent(in) :: other
	self%meshSeed = other%meshSeed
	self%initNest = other%initNest
	self%amrLimit = other%amrLimit
	self%geomKind = other%geomKind
	self%faceKind = other%faceKind
	call Copy(self%particles, other%particles)
	call Copy(self%edges, other%edges)
	call Copy(self%faces, other%faces)
end subroutine

!> @brief Computes the maximum variation of the Lagrangian parameter over all faces in a mesh.
!> @param[in] self Target mesh
!> @return Maximum over all faces of the variation of the Lagrangian parameter per face
function MaxLagrangianVariationPerFace( self )
	real(kreal) :: MaxLagrangianVariationPerFace
	type(PolyMesh2d), intent(in) :: self
	!
	integer(kint) :: i, j
	integer(kint) :: particleIndex
	real(kreal), dimension(3) :: maxX0, minX0, x0
	real(kreal) :: testVar
	type(STDIntVector) :: faceVerts
	
	MaxLagrangianVariationPerFace = 0.0_kreal
	
	do i = 1, self%faces%N
		if ( .NOT. self%faces%hasChildren(i) ) then
			particleIndex = self%faces%centerParticle(i)
			maxX0 = LagCoord(self%particles, particleIndex)
			minX0 = LagCoord(self%particles, particleIndex)
			
			!call CCWVerticesAroundFace( self, faceVerts, i )
			do j = 1, self%faceKind
				x0 = LagCoord(self%particles, self%faces%vertices(j,i))
				if ( x0(1) < minX0(1) ) minX0(1) = x0(1)
				if ( x0(2) < minX0(2) ) minX0(2) = x0(2)
				if ( x0(3) < minX0(3) ) minX0(3) = x0(3)
				if ( x0(1) > maxX0(1) ) maxX0(1) = x0(1)
				if ( x0(2) > maxX0(2) ) maxX0(2) = x0(2)
				if ( x0(3) > maxX0(3) ) maxX0(3) = x0(3)
			enddo
			
			testVar = sum(maxX0 - minX0)
			if ( testVar > MaxLagrangianVariationPerFace ) MaxLagrangianVariationPerFace = testVar
		endif
	enddo
end function 

!> @brief Locates the face containing an arbitrary query point.
!> This function assumes faces are convex polygons.
!> @todo This function should handle points outside the mesh appropriately.
!> 
!> @param[in] self Target mesh
!> @param[in] queryPt point in the spatial domain
!> @return Integer index of face (in a @ref Faces object) containing the query point.
function LocateFaceContainingPoint(self, queryPt)
	integer(kint) :: locateFaceContainingPoint
	type(PolyMesh2d), intent(in) :: self
	real(kreal), intent(in) :: queryPt(:)
	!
	integer(kint) :: treeStart
	integer(kint) :: walkStart
!	call LogMessage(log, DEBUG_LOGGING_LEVEL, trim(logkey)//" LocateFaceContainingPoint:", " entering.")
	LocateFaceContainingPoint = 0
	treeStart = nearestRootFace(self, queryPt)
	walkStart = locatePointTreeSearch(self, queryPt, treeStart)
	locateFaceContainingPoint = locatePointWalkSearch(self, queryPt, walkStart)
	if ( LocateFaceContainingPoint == 0 .OR. LocateFaceContainingPoint > self%faces%N ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL, trim(logkey)//"LocateFaceContainingPoint ERROR : ", "bad output.")
	endif
end function

!> @brief Finds the particle closest to an arbitrary query point.
!> @todo This function should handle points outside the mesh appropriately.
!>
!> @param[in] self Target mesh
!> @param[in] queryPt point in the spatial domain
!> @return Integer index of particle (in a @ref Particles object) closest to the query point.
function nearestParticle(self, queryPt)
	integer(kint) :: nearestParticle
	type(PolyMesh2d), intent(in) :: self
	real(kreal), intent(in) :: queryPt(:)
	!
	integer(kint) :: faceIndex, i
	type(STDIntVector) :: faceVerts
	real(kreal) :: dist, testDist
	
	faceIndex = LocateFaceContainingPoint(self, queryPt)
	call initialize(faceVerts)
	call CCWVerticesAroundFace( self, faceVerts, faceIndex)
	nearestParticle = self%faces%centerParticle(faceIndex)
	
	if ( self%geomKind == SPHERE_GEOM ) then
		dist  = SphereDistance( PhysCoord(self%particles, nearestParticle), queryPt)
		do i = 1, faceVerts%N
			testDist = SphereDistance( PhysCoord(self%particles, faceVerts%int(i)), queryPt)
			if ( testDist < dist ) then
				dist = testDist
				nearestParticle = faceVerts%int(i)
			endif
		enddo
	else
		dist = ChordDistance( PhysCoord(self%particles, nearestParticle), queryPt)
		do i = 1, faceVerts%N
			testDist = ChordDistance( PhysCoord( self%particles, faceVerts%int(i)), queryPt)
			if ( testDist < dist ) then
				dist = testDist
				nearestParticle = faceVerts%int(i)
			endif
		enddo
	endif
end function

!> @brief Returns an @ref STDIntVector containing indices of @ref Edges around a face, listed in counter-clockwise order.
!> 
!> @param[in] self Target mesh
!> @param[out] leafEdges list of CCW edges around face
!> @param[in] faceIndex index to a face in @ref Faces
subroutine CCWEdgesAroundFace( self, leafEdges, faceIndex )
	type(PolyMesh2d), intent(in) :: self
	type(STDIntVector), intent(out) :: leafEdges
	integer(kint), intent(in) :: faceIndex
	!
	integer(kint) :: i, j, nParentEdges
	type(STDIntVector) :: edgeLeaves(4)
	
!	call LogMessage(log,DEBUG_LOGGING_LEVEL,trim(logKey)//" CCWEdgesAroundFace : ", "entering.")
	
	if ( self%faceKind == TRI_PANEL ) then
		nParentEdges = 3
	elseif ( self%faceKind == QUAD_PANEL ) THEN
		nParentEdges = 4
	endif
	
	do i = 1, nParentEdges
		call GetLeafEdgesFromParent( self%edges, self%faces%edges(i,faceIndex), edgeLeaves(i))
	enddo
	
	call initialize(leafEdges)
	do i = 1, nParentEdges
		do j = 1, edgeLeaves(i)%N
			call leafEdges%pushBack(edgeLeaves(i)%int(j))
		enddo
	enddo
end subroutine

!> @brief Returns an @ref STDIntVector containing indices of vertex @ref Particles around a face, listed in counter-clockwise order.
!> 
!> @param[in] self Target mesh
!> @param[out] verts list of vertices around face, in CCW order
!> @param[in] faceIndex index to a face in @ref Faces
subroutine CCWVerticesAroundFace( self, verts, faceIndex )
	type(PolyMesh2d), intent(in) :: self
	type(STDIntVector), intent(out) :: verts
	integer(kint), intent(in) :: faceIndex
	!
	type(STDIntVector) :: leafEdges
	integer(kint) :: i
	
!	call LogMessage(log,DEBUG_LOGGING_LEVEL,trim(logKey)//" CCWVerticesAroundFace : ", "entering.")
	
	call CCWEdgesAroundFace(self, leafEdges, faceIndex)
	call initialize(verts)
	
	do i = 1, leafEdges%N
		if ( positiveEdge(self%edges, faceIndex, leafEdges%int(i)) ) then
			call verts%pushBack(self%edges%dest(leafEdges%int(i)))
		else
			call verts%pushBack(self%edges%orig(leafEdges%int(i)))
		endif
	enddo	
end subroutine

!> @brief Projects all particles in a mesh to the surface of a sphere with a user-specified radius
!> 
!> @param[inout] self Target mesh
!> @param[in] radius radius of sphere
subroutine ProjectParticlesToSphere( self, radius )
	type(PolyMesh2d), intent(inout) :: self
	real(kreal), intent(in) :: radius
	integer(kint) :: i
	real(kreal) :: scale
	
	if ( self%geomKind /= SPHERE_GEOM ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL, trim(logKey)//" ProjectToSphere ERROR : ", "invalid geomKind.")
		return
	endif
	
	do i = 1, self%particles%n
		scale = radius / sqrt( self%particles%x(i) * self%particles%x(i) + self%particles%y(i) * self%particles%y(i) + &
					 self%particles%z(i) * self%particles%z(i))
		self%particles%x(i) = self%particles%x(i) * scale
		self%particles%y(i) = self%particles%y(i) * scale
		self%particles%z(i) = self%particles%z(i) * scale
	enddo	
end subroutine

!> @brief Returns .TRUE. if a query point is outside the boundaries of a planar mesh.
!> 
!> @param[in] self Target mesh
!> @param[in] xVec position vector of query point
!> @return .TRUE. if query point is outside mesh boundaries
function pointIsOutsideMesh(self, xVec )
	logical(klog) :: pointIsOutsideMesh
	type(PolyMesh2d), intent(in) :: self
	real(kreal), intent(in) :: xVec(:)
	!
	integer(kint) :: faceIndex
	type(STDIntVector) :: faceEdges
	type(STDIntVector) :: boundaryEdges
	integer(kint) :: i
	real(kreal) :: centroid(3), v1(3), v2(3), p(3), q(3), reflection(3)
	real(kreal) :: interiorDist, exteriorDist
	
!	call LogMessage(log,DEBUG_LOGGING_LEVEL,trim(logKey)//" pointIsOutsideMesh : ", "entering.")
	
	pointIsOutsideMesh = .FALSE.
	faceIndex = locateFaceContainingPoint(self, xVec)

	call CCWEdgesAroundFace(self, faceEdges, faceIndex)
	call initialize(boundaryEdges)
	do i = 1, faceEdges%N
		if ( onBoundary(self%edges, faceEdges%int(i) ) ) then
			call boundaryEdges%pushback(faceEdges%int(i)) 
		endif
	enddo
!	call LogMessage(log,DEBUG_LOGGING_LEVEL,"pointIsOutsideMesh: n boundary edges = ", boundaryEdges%N)
	v1 = 0.0_kreal
	v2 = 0.0_kreal
	if ( boundaryEdges%N > 0 ) then
		centroid = FaceCentroid(self%faces, faceIndex, self%particles)
		interiorDist = sqrt(sum( (centroid-xVec)*(centroid-xVec)))
		do i = 1, boundaryEdges%N
			v1 = PhysCoord(self%particles, self%edges%orig( boundaryEdges%int(i) ))
			v2 = PhysCoord(self%particles, self%edges%dest( boundaryEdges%int(i) ))
			q = v2 - v1
			q = q/sqrt(sum(q*q))
			p = centroid - v1
			reflection = centroid - 2.0_kreal * (p - sum(p*q)*q)
			exteriorDist = sqrt(sum( (reflection-xVec)*(reflection-xVec)))
			if ( exteriorDist < interiorDist ) & 
				pointIsOutsideMesh = .TRUE.
				
!			print *, "centroid = ", centroid, ", reflection = ", reflection
		enddo		
	endif
end function

!> @brief Returns an @ref STDIntVector containing indices of @ref Faces adjacent to a face, listed in counter-clockwise order.
!> 
!> @param[in] self Target mesh
!> @param[out] adjFaces list of adjacent faces, in CCW order
!> @param[in] faceIndex index to a face in @ref Faces
subroutine CCWAdjacentFaces( self, adjFaces, faceIndex )
	type(PolyMesh2d), intent(in) :: self
	type(STDIntVector), intent(out) :: adjFaces
	integer(kint), intent(in) :: faceIndex
	!
	type(STDIntVector) :: leafEdges
	integer(kint) :: i
	
!	call LogMessage(log,DEBUG_LOGGING_LEVEL,trim(logKey)//" CCWAdjacentFaces : ", "entering.")
	
	call CCWEdgesAroundFace(self, leafEdges, faceIndex )
	
	call initialize(adjFaces)
	
	do i = 1, leafEdges%N
		if ( positiveEdge(self%edges, faceIndex, leafEdges%int(i)) ) then
			call adjFaces%pushBack(self%edges%rightFace(leafEdges%int(i)))
		else
			call adjFaces%pushBack(self%edges%leftFace(leafEdges%int(i)))
		endif
	enddo
end subroutine

!> @brief Returns an @ref STDIntVector containing indices of @ref Faces incident to a vertex, listed in counter-clockwise order.
!> @todo Verify this subroutine works correctly 
!>
!> @param[in] self Target mesh
!> @param[out] adjFaces list of adjacent faces, in CCW order
!> @param[in] vertexIndex index to a particle in @ref Particles object
subroutine CCWFacesAroundVertex( self, adjFaces, vertexIndex )
	type(PolyMesh2d), intent(in) :: self
	type(STDIntVector), intent(out) :: adjFaces
	integer(kint), intent(in) :: vertexIndex
	!
	integer(kint) :: i
	type(STDIntVector) :: incEdges
	
!	call LogMessage(log,DEBUG_LOGGING_LEVEL,trim(logKey)//" CCWFacesAroundVertex : ", "entering.")
	
	if ( self%particles%area(vertexIndex) > 0.0_kreal .OR. self%particles%isActive(vertexIndex) ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL, logkey//" CCWFacesAroundVertex ERROR : ", " vertices should have 0 area.")
		return
	endif
	
	call initialize(adjFaces)
	call initialize(incEdges,self%particles%incidentEdges(1:self%particles%nEdges(vertexIndex), vertexIndex))
	do i = 1, incEdges%N
		if ( self%edges%orig(incEdges%int(i)) == vertexIndex ) then
			call adjFaces%pushBack( self%edges%leftFace( incEdges%int(i)))	
		elseif ( self%edges%dest(incEdges%int(i)) == vertexIndex) then
			call adjFaces%pushBack( self%edges%rightFace(incEdges%int(i)))
		else
			write(logstring,'(A,I8,A)') " connectivity ERROR at vertex ", vertexIndex, ": has incident edge that does not connect."
			call LogMessage(log, ERROR_LOGGING_LEVEL,trim(logKey)//" CCWFacesAroundVertex :",logstring)
		endif
	enddo
end subroutine

!> @brief Finds the nearest root face (faces included in the mesh seed, with no parents) to an arbitrary query point.
!> Used to initialize a tree search of a @ref Faces data structure.
!> 
!> @param[in] self Target mesh
!> @param[in] queryPt position vector of query point
!> @return Index to @ref Faces of nearest root face
function nearestRootFace(self, queryPt)
	integer(kint) :: nearestRootFace
	type(PolyMesh2d), intent(in) :: self
	real(kreal), intent(in) :: queryPt(:)
	!
	integer(kint) :: i, nRootFaces
	real(kreal) :: dist, testDist, cntd(3)
	
	nearestRootFace = 1
	
!	call LogMessage(log,DEBUG_LOGGING_LEVEL,trim(logKey)//" nearestRootFace : ", "entering.")
	
	select case (self%meshSeed)
		case ( TRI_HEX_SEED, CUBED_SPHERE_SEED )
			nRootFaces = 6
		case ( QUAD_RECT_SEED, QUAD_RECT_PERIODIC_SEED, POLAR_DISC_SEED, BETA_PLANE_SEED)
			nRootFaces = 4
		case ( ICOS_TRI_SPHERE_SEED)
			nRootFaces = 20
	end select
	
	cntd = faceCentroid(self%faces, 1, self%particles)
	if ( self%geomKind == SPHERE_GEOM ) then
		dist = SphereDistance( cntd, queryPt)
		do i = 2, nRootFaces
			cntd = faceCentroid(self%faces, i, self%particles)
			testDist = SphereDistance( cntd, queryPt )
			if ( testDist < dist ) then
				dist = testDist
				nearestRootFace = i
			endif
		enddo
	else
		dist = sqrt(sum( (cntd-queryPt)*(cntd-queryPt)))
		do i = 2, nRootFaces
			cntd = faceCentroid(self%faces, i, self%particles)
			testDist = sqrt(sum( (cntd-queryPt)*(cntd-queryPt)))
			if ( testDist < dist ) then
				dist = testDist
				nearestRootFace = i
			endif
		enddo		
	endif
end function

!> @brief Determines the appropriate value of faceKind (as defined by @ref NumberKinds) based on a specified mesh seed.
!> @param[in] seed as defined in @ref NumberKinds
!> @return faceKind as defined in @ref NumberKinds
function FaceKindFromSeed( seed )
	integer(kint) :: FaceKindFromSeed
	integer(kint), intent(in) :: seed
	FaceKindFromSeed = 0
	select case (seed)
		case ( TRI_HEX_SEED, ICOS_TRI_SPHERE_SEED)
			FaceKindFromSeed = TRI_PANEL
		case ( QUAD_RECT_SEED, QUAD_RECT_PERIODIC_SEED, POLAR_DISC_SEED, CUBED_SPHERE_SEED, BETA_PLANE_SEED)
			FaceKindFromSeed = QUAD_PANEL
	end select
end function

!> @brief Determines the geometry (as defined by @ref NumberKinds) based on a specified mesh seed.
!> @param[in] seed as defined in @ref NumberKinds
!> @return geomKind as defined in @ref NumberKinds
function GeomKindFromSeed( seed )
	integer(kint) :: GeomKindFromSeed
	integer(kint), intent(in) :: seed
	GeomKindFromSeed = 0
	select case (seed)
		case (TRI_HEX_SEED, QUAD_RECT_SEED, QUAD_RECT_PERIODIC_SEED, POLAR_DISC_SEED, BETA_PLANE_SEED)
			GeomKindFromSeed = PLANAR_GEOM
		case (ICOS_TRI_SPHERE_SEED, CUBED_SPHERE_SEED)
			GeomKindFromSeed = SPHERE_GEOM
	end select
end function

!> @brief Performs a tree search of a @ref Faces object to locate the face containing a query point.
!> Used to initialize a walk search.
!> @param[in] self Target mesh
!> @param[in] queryPt position vector of query point
!> @param[in] index index (to @ref Faces) to start search
!> @return faceIndex index of face containing query point, as determined by a tree search.  Note that child particles
!> may move outside the boundaries of their parent faces, so this result may sometimes be incorrect.  
!> A walk search provides the correct result always.
recursive function locatePointTreeSearch( self, queryPt, index ) result(faceIndex)
	integer(kint) :: faceIndex
	type(PolyMesh2d), intent(in) :: self
	real(kreal), intent(in) :: queryPt(:)
	integer(kint), intent(in) :: index
	!
	real(kreal) :: dist, testDist, cntd(3)
	integer(kint) :: i, nearestChild
	
	faceIndex = 0
	
!	call LogMessage(log,DEBUG_LOGGING_LEVEL,trim(logKey)//" locatePointTreeSearch : ", "entering.")
	
	dist = 1.0d20
	nearestChild = 0
	cntd = 0.0_kreal
	if ( self%faces%hasChildren(index) ) then
		do i = 1, 4
			cntd = FaceCentroid(self%faces, self%faces%children(i,index), self%particles)
			if ( self%geomKind == SPHERE_GEOM ) then
				testDist = SphereDistance( queryPt, cntd )
			else
				testDist = sqrt(sum( (queryPt-cntd)*(queryPt-cntd)))
			endif
			if ( testDist < dist ) then
				nearestChild = self%faces%children(i,index)
				dist = testDist
			endif
		enddo
		faceIndex = locatePointTreeSearch( self, queryPt, nearestChild)
	else
		faceIndex = index
		return
	endif	
end function

!> @brief Performs a walk search of a @ref Faces object to locate the face containing a query point.
!> Accelerated by using the output of a tree search for its first guess.
!> @param[in] self Target mesh
!> @param[in] queryPt position vector of query point
!> @param[in] index index (to @ref Faces) to start search
!> @return faceIndex index of face containing query point
recursive function locatePointWalkSearch(self, queryPt, index) result(faceIndex)
	integer(kint) :: faceIndex
	type(PolyMesh2d), intent(in) :: self
	real(kreal), intent(in) :: queryPt(:)
	integer(kint), intent(in) :: index
	!
	real(kreal) :: dist, testDist, cntd(3)
	type(STDIntVector) :: adjFaces
	integer(kint) :: currentFace, i
	
!	call LogMessage(log,DEBUG_LOGGING_LEVEL,trim(logKey)//" locatePointWalkSearch : ", "entering.")
	
	if ( self%faces%hasChildren(index) ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL, logkey//" locatePointWalkSearch ERROR:"," expected a leaf face.")
		return
	endif

	faceIndex = 0
	currentFace = index
	cntd = FaceCentroid(self%faces, index, self%particles)	
	if ( self%geomKind == SPHERE_GEOM ) then
		dist = SphereDistance( cntd, queryPt)
	else
		dist = sqrt(sum( (cntd - queryPt)*(cntd-queryPt)))
	endif
	
	call CCWAdjacentFaces(self, adjFaces, currentFace)
	
	do i = 1, adjFaces%N
		if ( adjFaces%int(i) > 0 ) then
			cntd = FaceCentroid( self%faces, adjFaces%int(i), self%particles )
			if ( self%GeomKind == SPHERE_GEOM ) then
				testDist = SphereDistance( cntd, queryPt)
			else
				testDist = sqrt(sum( (cntd - queryPt)*(cntd-queryPt) ))
			endif	
		
			if ( testDist < dist ) then
				dist = testDist
				currentFace = adjFaces%int(i)
			endif
		endif
	enddo
	
	if ( currentFace == index ) then
		faceIndex = currentFace
		return
	else
		faceIndex = locatePointWalkSearch( self, queryPt, currentFace )
	endif	
end function

!> @brief Computes the number of vertices in a uniformly refined mesh.  
!> Used to compute memory allocation requirements.
!> @param[in] self Target mesh
!> @param[in] initNest initial refinement level
!> @return Number of vertices in a mesh
function nVerticesInMesh( self, initNest )
	integer(kint) :: nVerticesInMesh
	type(PolyMesh2d), intent(in) :: self
	integer(kint), intent(in) :: initNest
	!
	integer(kint) :: i, result
	
	result = 0
	select case ( self%meshSeed)
 		case ( TRI_HEX_SEED ) 
			do i = 2**initNest + 1, 2**(initNest+1)
				result = result + i
			enddo
			result = 2*result + 2**(initNest+1)+1
		case ( QUAD_RECT_SEED, QUAD_RECT_PERIODIC_SEED, BETA_PLANE_SEED) 
			result = 3
			do i = 1, initNest
				result = result + 2**i
			enddo
			result = result*result
		case ( ICOS_TRI_SPHERE_SEED) 
			result = 2 + 10 * 4 ** initNest
		case ( CUBED_SPHERE_SEED ) 
			result = 2 + 6 * 4 ** initNest
	end select
	nVerticesInMesh = result
end function

!> @brief Computes the number of faces in a uniformly refined mesh.  
!> Used to compute memory allocation requirements.
!> @param[in] self Target mesh
!> @param[in] initNest initial refinement level
!> @return Number of faces in a mesh
function nFacesInMesh( self, initNest )
	integer(kint) :: nFacesInMesh
	type(PolyMesh2d), intent(in) :: self
	integer(kint), intent(in) :: initNest
	
	nFacesInMesh = 0
	select case (self%meshSeed)
		case (TRI_HEX_SEED)
			nFacesInMesh = 6 * 4 ** initNest
		case (QUAD_RECT_SEED, QUAD_RECT_PERIODIC_SEED, BETA_PLANE_SEED)
			nFacesInMesh = 4 * 4 ** initNest
		case (ICOS_TRI_SPHERE_SEED)
			nFacesInMesh = 20 * 4 ** initNest
		case (CUBED_SPHERE_SEED)
			nFacesInMesh = 6 * 4 ** initNest
	end select
end function

!> @brief Uses the number of vertices and the number of faces (via Euler's formulas) to determine the the number of edges in a uniformly refined mesh.  
!> Used to compute memory allocation requirements.
!> @param[in] self Target mesh
!> @param[in] nVertices number of vertices in a mesh
!> @param[in] nFaces number of faces in a mesh
!> @return Number of edges in a mesh
function nEdgesInMesh( self, nVertices, nFaces )
	integer(kint) :: nEdgesInMesh
	type(PolyMesh2d), intent(in) :: self
	integer(kint), intent(in) :: nVertices
	integer(kint), intent(in) :: nFaces
	nEdgesInMesh = 0
	select case (self%geomKind)
		case (PLANAR_GEOM)
			nEdgesInMesh = nFaces + nVertices - 1
		case (SPHERE_GEOM)
			nEdgesInMesh = nFaces + nVertices - 2
	end select
end function

!> @brief Initializes a mesh using a mesh seed .dat file to define the root set of @ref Particles, @ref Edges, and @ref Faces.
!> Mesh seeds are defined in both @ref NumberKinds and via the corresponding .dat file.
!> Multiplies particle postions by a specified amount to change the spatial domain (all mesh seeds approximately cover 
!> the unit sphere of their embedded space).  
!> 
!> @param[inout] self Target mesh to be initialized
!> @param[in] ampFactor 
subroutine initializeMeshFromSeed( self, ampFactor )
	type(PolyMesh2d), intent(inout) :: self
	real(kreal), intent(in) :: ampFactor
	!
	character(len=56) :: seedFilename
	integer(kint) :: i
	integer(kint) :: nSeedParticles, nSeedEdges, nSeedFaces, nSeedVerts
	real(kreal), allocatable :: seedXYZ(:,:)
	integer(kint), allocatable :: seedEdgeOrigs(:), seedEdgeDests(:), seedEdgeLefts(:), seedEdgeRights(:)
	integer(kint), allocatable :: seedFaceVerts(:,:), seedFaceEdges(:,:), seedVertexDegree(:)
	
	select case ( self%meshSeed )
		case (TRI_HEX_SEED)
			nSeedParticles = 13
			nSeedEdges = 12
			nSeedFaces = 6
			nSeedVerts = 7
!			self%faceKind = TRI_PANEL
!			self%geomKind = PLANAR_GEOM
			seedFilename = "triHexSeed.dat"
		case (QUAD_RECT_SEED)
			nSeedParticles = 13
			nSeedEdges = 12
			nSeedFaces = 4
			nSeedVerts = 9
!			self%faceKind = QUAD_PANEL
!			self%geomKind = PLANAR_GEOM
			seedFilename = "quadRectSeed.dat"
		case (POLAR_DISC_SEED)
			call LogMessage(log, ERROR_LOGGING_LEVEL,logkey," initMeshFromSeed ERROR : seed not implemented.")
			return
		case (QUAD_RECT_PERIODIC_SEED)
			nSeedParticles = 13
			nSeedEdges = 12
			nSeedFaces = 4
			nSeedVerts = 9
			call LogMessage(log, ERROR_LOGGING_LEVEL,logkey," initMeshFromSeed ERROR : seed not implemented.")
			return
		case (ICOS_TRI_SPHERE_SEED)
			nSeedParticles = 32
			nSeedEdges = 30
			nSeedFaces = 20
			nSeedVerts = 12
!			self%faceKind = TRI_PANEL
!			self%geomKind = SPHERE_GEOM
			seedFilename = "icosTriSeed.dat"
		case (CUBED_SPHERE_SEED)
			nSeedParticles = 14
			nSeedEdges = 12
			nSeedFaces = 6
			nSeedVerts = 8
!			self%faceKind = QUAD_PANEL
!			self%geomKind = SPHERE_GEOM
			seedFilename = "cubedSphereSeed.dat"
		case (BETA_PLANE_SEED)
			nSeedParticles = 13
			nSeedEdges = 12
			nSeedFaces = 4
			nSeedVerts = 9
			seedFilename = "betaPlaneSeed.dat"
	end select
	
	
	allocate(seedXYZ(3,nSeedParticles))
	allocate(seedEdgeOrigs(nSeedEdges))
	allocate(seedEdgeDests(nSeedEdges))
	allocate(seedEdgeLefts(nSeedEdges))
	allocate(seedEdgeRights(nSeedEdges))
	allocate(seedVertexDegree(nSeedVerts))
	if ( self%faceKind == TRI_PANEL ) then
		allocate(seedFaceVerts(3,nSeedFaces))
		allocate(seedFaceEdges(3,nSeedFaces))
	elseif ( self%faceKind == QUAD_PANEL ) then
		allocate(seedFaceVerts(4,nSeedFaces))
		allocate(seedFaceEdges(4,nSeedFaces))
	endif
	
	call readSeedFile(self, seedFilename, nSeedParticles, nSeedEdges, nSeedFaces, seedXYZ, seedEdgeOrigs, seedEdgeDests, &
					  seedEdgeLefts, seedEdgeRights, seedFaceVerts, seedFaceEdges, seedVertexDegree)
					  
	if ( self%meshSeed == ICOS_TRI_SPHERE_SEED ) then
		do i = 1, 20
			seedXYZ(:,12+i) = SphereTriCenter( seedXYZ(:,seedFaceVerts(1,i)), seedXYZ(:,seedFaceVerts(2,i)), &
											   seedXYZ(:,seedFaceVerts(3,i)) )
		enddo
	elseif ( self%meshSeed == BETA_PLANE_SEED ) then
		do i = 1, 13
			seedXYZ(1,i) = 0.5_kreal * seedXYZ(1,i) + 0.5_kreal
			seedXYZ(2,i) = 0.5_kreal * seedXYZ(2,i) 
		enddo		
	endif
	
	seedXYZ = ampFactor * seedXYZ
	
	!
	!	initialize mesh with seed data
	!
	do i = 1, nSeedParticles
		call InsertParticle( self%particles, seedXYZ(:,i), seedXYZ(:,i))
	enddo
	if ( self%particles%N /= nSeedParticles ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL, logkey//" initMeshFromSeed ERROR : "," particles%N.")
	endif
	do i = 1, nSeedEdges
		call InsertEdge( self%Edges, self%particles, seedEdgeOrigs(i), seedEdgeDests(i), seedEdgeLefts(i), seedEdgeRights(i))
	enddo
	if ( self%edges%N /= nSeedEdges ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL, logkey//" initMeshFromSeed ERROR : "," edges%N.")
	endif
	do i = 1, nSeedFaces
		call InsertFace( self%Faces, nSeedVerts + i, seedFaceVerts(:,i), seedFaceEdges(:,i))
		call MakeParticleActive( self%particles, nSeedVerts +i)
	enddo
	if ( self%faces%N /= nSeedFaces ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL, logkey//" initMeshFromSeed ERROR : "," faces%N.")
	endif
	self%faces%N_Active = nSeedFaces
	if ( self%faceKind == TRI_PANEL ) then
		do i = 1, nSeedFaces
			self%particles%area( self%faces%centerParticle(i) ) = TriFaceArea( self%faces, i, self%particles)
		enddo
	elseif ( self%faceKind == QUAD_PANEL ) then
		do i = 1, nSeedFaces
			self%particles%area( self%faces%centerParticle(i) ) = QuadFaceArea( self%faces, i, self%particles)
		enddo
	endif
	!
	!	initialize dual information
	!
	do i = 1, nSeedVerts
		if ( self%particles%nEdges(i) /= seedVertexDegree(i) ) then
			write(logstring,*) " dual connectivity error at vertex ", i
			call LogMessage(log, WARNING_LOGGING_LEVEL,logkey//" initFromSeed WARNING: ",logstring)
		endif
		call SortIncidentEdgesAtParticle( self%particles, i )
	enddo

	deallocate(seedXYZ)
	deallocate(seedEdgeOrigs)
	deallocate(seedEdgeDests)
	deallocate(seedEdgeLefts)
	deallocate(seedEdgeRights)
	deallocate(seedFaceVerts)
	deallocate(seedFaceEdges)
end subroutine

!> @brief Reads data from a mesh seed .dat file.
!> 
!> @param[in] self Target mesh
!> @param[in] seedFilename filename of .dat file
!> @param[in] nParticles number of particles in mesh seed
!> @param[in] nEdges number of edges in mesh seed
!> @param[in] nFaces number of faces in mesh seed
!> @param[out] seedXYZ position vectors of mesh seed particles
!> @param[out] seedEdgeOrigs indices of origin particles for each edge
!> @param[out] seedEdgeDests indices of destination particles for each edge
!> @param[out] seedEdgeLefts indices of left faces for each edge
!> @param[out] seedEdgeRights indices of right faces for each edge
!> @param[out] seedFaceVerts indices of vertex particles for each face
!> @param[out] seedFaceEdges indices of edges for each face
!> @param[out] seedVertexDegree degree (number of incident edges) of each vertex in mesh seed
subroutine readSeedFile( self, seedFilename, nParticles, nEdges, nFaces, seedXYZ, &
						 seedEdgeOrigs, seedEdgeDests, seedEdgeLefts, seedEdgeRights, &
						 seedFaceVerts, seedFaceEdges, seedVertexDegree )
	type(PolyMesh2d), intent(in) :: self
	character(len=*), intent(in) :: seedFilename
	integer(kint), intent(in) :: nParticles, nEdges, nFaces
	real(kreal), intent(out) :: seedXYZ(3,nParticles)
	integer(kint), intent(out) :: seedEdgeOrigs(nEdges)
	integer(kint), intent(out) :: seedEdgeDests(nEdges)
	integer(kint), intent(out) :: seedEdgeLefts(nEdges)
	integer(kint), intent(out) :: seedEdgeRights(nEdges)
	integer(kint), intent(inout) :: seedFaceVerts(:,:)
	integer(kint), intent(inout) :: seedFaceEdges(:,:)
	integer(kint), intent(out) :: seedVertexDegree(nParticles - nFaces)
	!
	character(len=56) :: tempString							 
	integer(kint) :: i, readStat
	
	open(unit=READ_UNIT,file=seedFilename,STATUS='OLD',ACTION='READ',iostat=readStat)
	if ( readStat /= 0 ) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,logKey," cannot open seed file.")
		return
	endif
	!
	!	read root particle positions
	!
	read(READ_UNIT,*) tempString
	if ( self%geomKind == PLANAR_GEOM ) then
		do i = 1, nParticles
			read(READ_UNIT,*) seedXYZ(1,i), seedXYZ(2,i)
		enddo
	else
		if ( self%meshSeed == ICOS_TRI_SPHERE_SEED ) then
			do i = 1, 12
				read(READ_UNIT,*) seedXYZ(1,i), seedXYZ(2,i), seedXYZ(3,i)
			enddo
		else
			do i = 1, nParticles
				read(READ_UNIT,*) seedXYZ(1,i), seedXYZ(2,i), seedXYZ(3,i)
			enddo
		endif
	endif
	!
	!	read root edges
	!
	read(READ_UNIT,*) tempString
	do i = 1, nEdges
		read(READ_UNIT,*) seedEdgeOrigs(i), seedEdgeDests(i), seedEdgeLefts(i), seedEdgeRights(i)
	enddo
	! 	adjust to fortran (base index = 1, not 0)
	seedEdgeOrigs = seedEdgeOrigs + 1
	seedEdgeDests = seedEdgeDests + 1
	seedEdgeLefts = seedEdgeLefts + 1
	seedEdgeRights = seedEdgeRights + 1
	!
	!	read root faces
	!
	read(READ_UNIT,*) tempString
	if ( self%faceKind == TRI_PANEL ) then
		do i = 1, nFaces
			read(READ_UNIT,*) seedFaceVerts(1,i), seedFaceVerts(2,i), seedFaceVerts(3,i)
		enddo
		read(READ_UNIT,*) tempString
		do i =1, nFaces
			read(READ_UNIT,*) seedFaceEdges(1,i), seedFaceEdges(2,i), seedFaceEdges(3,i)
		enddo
	elseif ( self%faceKind == QUAD_PANEL ) then
		do i = 1, nFaces
			read(READ_UNIT,*) seedFaceVerts(1,i), seedFaceVerts(2,i), seedFaceVerts(3,i), seedFaceVerts(4,i)
		enddo
		read(READ_UNIT,*) tempString
		do i =1, nFaces
			read(READ_UNIT,*) seedFaceEdges(1,i), seedFaceEdges(2,i), seedFaceEdges(3,i), seedFaceEdges(4,i)
		enddo
	endif
	!	adjust to fortran (base index = 1)
	seedFaceEdges = seedFaceEdges + 1
	seedFaceVerts = seedFaceVerts + 1
	!
	!	read vertex degree
	!
	read(READ_UNIT,*) tempString
	do i = 1, nParticles - nFaces
		read(READ_UNIT,*) seedVertexDegree(i)
	enddo
	
	close(READ_UNIT)
end subroutine

!> @brief Writes a mesh to a Matlab-readable .m file
!> 
!> @param[in] self Mesh to output
!> @param[in] fileunit Integer unit of .m file
subroutine WriteMeshToMatlab( self, fileunit )
	type(PolyMesh2d), intent(in) :: self
	integer(kint), intent(in) :: fileunit
	
	call WriteParticlesToMatlab(self%particles, fileunit)
	call WriteEdgesToMatlab(self%edges, fileunit)
	call WriteFacesToMatlab(self%faces, fileunit)
end subroutine

!> @brief Outputs basic mesh info to a @ref Logger object.
!> @param[in] self Target mesh
!> @param[inout] aLog logger
subroutine logStatsPrivate(self, aLog)
	type(PolyMesh2d), intent(in) :: self
	type(Logger), intent(inout) :: aLog
	call LogMessage(aLog, TRACE_LOGGING_LEVEL, logKey, "PolyMesh Stats:" )
	call StartSection(aLog)
	call MeshSeedString(logString, self%meshSeed)
	call LogMessage(aLog, TRACE_LOGGING_LEVEL, "meshSeed = ", trim(logString))
	if ( self%geomKind == PLANAR_GEOM ) then
		call LogMessage(aLog, TRACE_LOGGING_LEVEL, "geomKind = ", "PLANAR_GEOM")
	elseif ( self%geomKind == SPHERE_GEOM ) then
		call LogMessage(aLog, TRACE_LOGGING_LEVEL, "geomKind = ", "SPHERE_GEOM")
	elseif ( self%geomKind == EUCLIDEAN_3D ) then
		call LogMessage(aLog, TRACE_LOGGING_LEVEL, "geomKind = ", "EUCLIDEAN_3D")
	else
		call LogMessage(aLog, TRACE_LOGGING_LEVEL, "geomKind = ", "invalid geomKind.")
	endif
	call LogMessage(aLog, TRACE_LOGGING_LEVEL, "initNest = ", self%initNest)
	call LogMessage(aLog, TRACE_LOGGING_LEVEL, "amrLimit = ", self%amrLimit)
	call LogMessage(aLog, TRACE_LOGGING_LEVEL, "max Edge len = ", MaxEdgeLength(self%edges, self%particles) )
	call LogMessage(aLog, TRACE_LOGGING_LEVEL, "min Edge len = ", MinEdgeLength(self%edges, self%particles) )
	call LogMessage(aLog, TRACE_LOGGING_LEVEL, "avg Edge len = ", AvgEdgeLength(self%edges, self%particles) )
	call LogStats(self%particles, aLog)
	call LogStats(self%edges, aLog)
	call LogStats(self%faces, aLog)
	call EndSection(aLog)
end subroutine

!> @brief Outputs a mesh to a legacy formatted .vtk ASCII file.
!> Does not output any @ref Field data.
!> @param[in] self Target mesh
!> @param[in] fileunit Integer unit of .vtk file
!> @param[in] title Optional title for vtk output
subroutine WriteMeshToVTKPolyData( self, fileunit, title )
	type(PolyMesh2d), intent(in) :: self
	integer(kint), intent(in) :: fileunit
	character(len=*), intent(in), optional :: title
	
	if ( present(title) ) then
		call WriteVTKPoints(self%particles, fileunit, title)
	else
		call WriteVTKPoints(self%particles, fileunit)
	endif
	
	call WriteFacesToVTKPolygons( self%faces, fileunit )
	call WriteVTKPointDataSectionHeader(fileunit, self%particles%N)
	call WriteVTKLagCoords(self%particles, fileunit )
	
	call WriteFaceAreaToVTKCellData(self%faces, self%particles, fileunit)	
end subroutine

!> @brief Prints detailed mesh information to console.
!> Used for debugging.
!> @param[in] self Target mesh 
subroutine PrintDebugPrivate( self ) 
	type(PolyMesh2d), intent(in) :: self
	print *, "PolyMesh2d DEBUG info : "
	print *, "meshSeed = ", self%meshSeed
	print *, "geomKind = ", self%geomKind
	print *, "initNest = ", self%initNest
	print *, "amrLimit = ", self%amrLimit
	call PrintDebugInfo( self%particles )
	call PrintDebugInfo( self%edges )
	call PrintDebugInfo( self%faces )
end subroutine

!> @brief Converts the integer definitions of mesh seeds (from @ref NumberKinds) to their corresponding strings.
!> @param[out] seedString
!> @param[in] msInt mesh seed integer 
subroutine MeshSeedString( seedString, msInt )
	character(len=*), intent(inout) :: seedString
	integer(kint), intent(in) :: msInt
	select case (msInt)
		case (TRI_HEX_SEED )
			write(seedString,'(A)') "TRI_HEX_SEED"
		case (QUAD_RECT_SEED)
			write(seedString,'(A)') "QUAD_RECT_SEED"
		case (QUAD_RECT_PERIODIC_SEED)
			write(seedString,'(A)') "QUAD_RECT_PERIODIC_SEED"
		case (POLAR_DISC_SEED)
			write(seedString,'(A)') "POLAR_DISC_SEED"
		case (ICOS_TRI_SPHERE_SEED)
			write(seedString,'(A)') "ICOS_TRI_SPHERE_SEED"
		case (CUBED_SPHERE_SEED)
			write(seedString,'(A)') "CUBED_SPHERE_SEED"
		case (BETA_PLANE_SEED)
			write(seedString,'(A)') "BETA_PLANE_SEED"
		case default
			write(seedString,'(A)') "invalid seed"
	end select
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