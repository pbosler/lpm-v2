module FacesModule
!> @file Faces.f90
!> Provides a primitive data structure and methods for creating faces of polyhedral meshes.
!> @author Peter Bosler, Sandia National Laboratories Center for Computing Research
!>
!>
!> @defgroup Faces Faces
!> @brief Faces of polyhedral meshes have centers and vertices via the @ref Particles and edges via the @ref Edges.
!> 
!> The faces data structures is a "structure of arrays," so that all information about face i is located at index i in the 
!> relevant array.  
!> For example, the index to the particle at the center of face i `aFaces%%centerParticle(i)` 
!> and the vertices of face i can be accessed by `aFaces%%vertices(:,i)`. 
!> 
!> The faces data structure maintains arrays of pointers (integer indices) into @ref Particles and @ref Edges objects
!> to avoid reproducing the data contained by those structures. 
!> 
!> In addition the faces are organized into a quadtree, defined by the `hasChildren`, `children`, and `parent` arrays,
!> to facilitate faster searching through the data structure.
!> 
!> @{
use NumberKindsModule
use STDIntVectorModule
use LoggerModule
use ParticlesModule
use EdgesModule
use PlaneGeomModule
use SphereGeomModule

implicit none

private
public Faces
public New, Delete, Copy
public FaceArea, FaceCenterPhysCoord, FaceCenterLagCoord
public DivideQuadFace, DivideTriFace, InsertFace
public WriteFacesToVTKPolygons, WriteFacesToMatlab, WriteFaceAreaToVTKCellData
public positiveEdge
public FaceCentroid, TriFaceArea, QuadFaceArea
public LogStats, PrintDebugInfo
public SharedEdge, ParticleOppositeTriEdge
public CountParents

!> @brief Each faces knows the indices (to @ref Particles) of its vertices and center and the indices (to @ref Edges)
!> of its edges.  
!> 
!> Designed to accomodate quadrilateral and triangular faces.  
!> Faces are recursively divided to provide desired spatial resolution, and stored in the induced quadtree.
!>
type Faces
	integer(kint), allocatable :: centerParticle(:) !< indices to faces' center particles in a particlesmodule::particles object
	integer(kint), allocatable :: vertices(:,:)  !< indices to faces' vertex particles in a particlesmodule::particles object
	integer(kint), allocatable :: edges(:,:) !< indices to faces' edges in an edgesmodule::edges object
	logical(klog), allocatable :: hasChildren(:) !< hasChildren(i) is .TRUE. if face i has been divided
	integer(kint), allocatable :: children(:,:) !< indices to child faces in a faces object
	integer(kint), allocatable :: parent(:) !< indices to parent faces in a faces object
	integer(kint) :: faceKind !< identifies type of face (e.g., triangular, quadrilateral, etc.)
	integer(kint) :: N !< Number of faces currently in use
	integer(kint) :: N_Active !< Number of undivided faces; these define the spatial discretization
	integer(kint) :: N_Max !< Maximum number of faces allowed in memory
	
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
	module procedure LogStatsPrivate
end interface

interface PrintDebugInfo
	module procedure PrintDebugPrivate
end interface

interface CountParents
	module procedure countParentFaces
end interface 

!
!----------------
! Logging
!----------------
!
logical(klog), save :: logInit = .FALSE.
type(Logger) :: log
character(len=28), save :: logKey = 'Faces'
integer(kint), parameter :: logLevel = DEBUG_LOGGING_LEVEL
character(len=MAX_STRING_LENGTH) :: logstring

contains

!> @brief Allocates memory for a new faces object.  New arrays are set to zero, and must be initialized separately.
!> @todo Move faceKind constants from ::numberkindsmodule to ::facesmodule? 
!>
!> @param self Target faces object
!> @param faceKind Face kind identifier (e.g., triangular, quadrilateral, etc.), as definded in NumberKindsModule
!> @param nMax Max amount of memory to allocate
subroutine newPrivate( self, faceKind, nMax )
	type(Faces), intent(out) :: self
	integer(kint), intent(in) :: faceKind
	integer(kint), intent(in) :: nMax

	if ( .NOT. logInit ) call InitLogger(log, procRank)

	if ( nMax <= 0 ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL, logKey, " NewFaces ERROR : invalid nMax.")
		return
	endif	

	if ( faceKind == QUAD_PANEL ) then
		allocate(self%vertices(4,nMax))
		allocate(self%edges(4,nMax))
		self%vertices = 0 
		self%edges = 0
	elseif ( faceKind == TRI_PANEL ) then
		allocate(self%vertices(3,nMax))
		allocate(self%edges(3,nMax))
		self%vertices = 0 
		self%edges = 0
	else
		call LogMessage(log, ERROR_LOGGING_LEVEL, logKey, " NewFaces ERROR : invalid face kind.")
		return
	endif
	
	allocate(self%centerParticle(nMax))
	self%centerParticle = 0
	allocate(self%children(4,nMax))
	self%children = 0
	allocate(self%hasChildren(nMax))
	self%hasChildren = .FALSE.
	allocate(self%parent(nMax))
	self%parent = 0
	self%N = 0
	self%N_Active = 0
	self%N_Max = nMax
	self%faceKind = faceKind
end subroutine

!> @brief Deletes and frees memory associated with a faces object
!> @param self Target faces object
subroutine deletePrivate(self)
	type(Faces), intent(inout) :: self
	if ( allocated(self%vertices)) deallocate(self%vertices)
	if ( allocated(self%edges)) deallocate(self%edges)
	if ( allocated(self%centerParticle)) deallocate(self%centerParticle)
	if ( allocated(self%children)) deallocate(self%children)
	if ( allocated(self%hasChildren)) deallocate(self%hasChildren)
	if ( allocated(self%parent)) deallocate(self%parent)
end subroutine

!> @brief Performs a deep copy of one faces object into another.  
!> Both objects must have been allocated correctly before calling this subroutine.
!> @param dest Target faces object
!> @param source Source faces object
subroutine copyPrivate( dest, source )
	type(Faces), intent(inout) :: dest
	type(Faces), intent(in) :: source
	!
	integer(kint) :: j
	
	if ( dest%N_Max < source%N ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL, logkey, " CopyFaces ERROR : not enough memory.")
		return
	endif
	if ( dest%faceKind /= source%faceKind) then
		call LogMessage(log, ERROR_LOGGING_LEVEL, logkey, " CopyFaces ERROR : face kind mismatch.")
		return
	endif
	dest%N = source%N
	dest%N_Active = source%N_Active
	do j = 1, source%N
		dest%vertices(:,j) = source%vertices(:,j)
		dest%edges(:,j) = source%vertices(:,j)
		dest%centerParticle(j) = source%centerParticle(j)
		dest%children(:,j) = source%children(:,j)
		dest%hasChildren(j) = source%hasChildren(j)
		dest%parent(j) = source%parent(j)
	enddo
end subroutine

!> @brief Prints detailed information about the contents of a faces object to the console.
!> Used for debugging.
!> @param self Target faces object
subroutine PrintDebugPrivate( self )
	type(Faces), intent(in) :: self
	integer(kint) :: i, j
	print *, " Faces DEBUG info : " 
	print *, "faces.N  = ", self%N
	print *, "faces.N_Max = ", self%N_Max
	print *, "faces.N_Active = ", self%N_Active
	print *, "faces.faceKind = ", self%faceKind
	print *, "n divided faces = ", count(self%hasChildren)
	print *, "faces.centerParticle = "
	do i = 1, self%N_Max
		print *, self%centerParticle(i)
	enddo
	print *, " faces.vertices = "
	do i = 1, self%N_Max
		do j = 1, size(self%vertices,1)
			write(6,'(I6)',advance='NO') self%vertices(j,i)
		enddo
		print *, " "
	enddo
	print *, " "
	print *, "faces.edges = "
	do i = 1, self%N_Max
		do j = 1, size(self%edges,1)
			write(6,'(I6)',advance='NO') self%edges(j,i)
		enddo
		print *, " "
	enddo
	print *, " "
	print *, "faces tree = "
	do i = 1, self%N_Max
		if ( self%hasChildren(i) ) then
			write(6,'(A,4X)',advance ='NO') 'T'
		else
			write(6,'(A,4X)',advance='NO' ) 'F'
		endif
		do j = 1, 4
			write(6, '(I6)', advance='NO') self%children(j,i)
		enddo
		print *, self%parent(i) 
	enddo	
end subroutine

!> @brief Writes basic info about a faces object to a specified loggermodule::Logger
!> @param self Target faces object
!> @param aLog Logger object to handle output
subroutine LogStatsPrivate(self, aLog)
	type(Faces), intent(in) :: self
	type(Logger), intent(inout) :: aLog 
	call LogMessage(aLog, TRACE_LOGGING_LEVEL, logkey, " Faces stats : ")
	call StartSection(aLog)
	call LogMessage(alog, TRACE_LOGGING_LEVEL, "faces.N = ", self%N)
	call LogMessage(aLog, TRACE_LOGGING_LEVEL, "faces.N_Max = ", self%N_Max )
	call LogMessage(aLog, TRACE_LOGGING_LEVEL, "faces.N_Active = ", self%N_Active)
	call LogMessage(aLog, TRACE_LOGGING_LEVEL, "n divided faces = ", count(self%hasChildren))
	if ( self%faceKind == TRI_PANEL ) then
		call LogMessage(aLog, TRACE_LOGGING_LEVEL, "faceKind = ", "TRI_PANEL")
	elseif ( self%faceKind == QUAD_PANEL ) then
		call LogMessage(aLog, TRACE_LOGGING_LEVEL, "faceKind = ", "QUAD_PANEL")
	else
		call LogMessage(aLog, TRACE_LOGGING_LEVEL, "faceKind = ", "invalid faceKind.")
	endif
	call EndSection(aLog)
end subroutine

!> @brief Returns the physical coordinates of a face's center particle
!> @param self Target faces object
!> @param faceIndex Index to target face
!> @param aParticles particlesmodule::particles object associated with this set of faces
!> @return Coordinate vector of face's center particle's position
function FaceCenterPhysCoord( self, faceIndex, aparticles)
	real(kreal) :: FaceCenterPhysCoord(3)
	type(Faces), intent(in) :: self
	integer(kint) :: faceIndex
	type(Particles), intent(in) :: aParticles
	
	FaceCenterPhysCoord(1) = aParticles%x(self%centerParticle(faceIndex))
	FaceCenterPhysCoord(2) = aParticles%y(self%centerParticle(faceIndex))
	if ( allocated(aParticles%z) ) then
		FaceCenterPhysCoord(3) = aParticles%z(self%centerParticle(faceIndex))
	else
		FaceCenterPhysCoord(3) = 0.0_kreal
	endif
end function

!> @brief Returns the Lagrangian coordinates of a face's center particle
!> @param self Target faces object
!> @param faceIndex Index to target face
!> @param aParticles particlesmodule::particles object associated with this set of faces
!> @return Coordinate vector of face's center particle's Lagrangian coordinate
function FaceCenterLagCoord( self, faceIndex, aparticles)
	real(kreal) :: FaceCenterLagCoord(3)
	type(Faces), intent(in) :: self
	integer(kint) :: faceIndex
	type(Particles), intent(in) :: aParticles
	
	FaceCenterLagCoord(1) = aParticles%x0(self%centerParticle(faceIndex))
	FaceCenterLagCoord(2) = aParticles%y0(self%centerParticle(faceIndex))
	if ( allocated(aParticles%z0) ) then
		FaceCenterLagCoord(3) = aParticles%z0(self%centerParticle(faceIndex))
	else
		FaceCenterLagCoord(3) = 0.0_kreal
	endif
end function

!> @brief Counts the number of levels in the faces quadtree above a particular face
!> @param self Target Faces object
!> @param index Target face
!> @return Number of levels in quadtree above face(i)(
function countParentFaces( self, index )
	integer(kint) :: countParentFaces
	type(Faces), intent(in) :: self
	integer(kint), intent(in) :: index
	!
	logical(klog) :: keepGoing
	integer(kint) :: parentIndex
	
	countParentFaces = 0
	keepGoing = ( self%parent(index) > 0 )
	parentIndex = self%parent(index)
	do while ( keepGoing )
		countParentFaces = countParentFaces + 1
		parentIndex = self%parent( parentIndex )
		keepGoing = ( parentIndex > 0 )
	enddo
end function 

!> @brief Counts the number of shared edges between two faces
!> Outputs an error of this number is > 1.
!> @param self Target faces object
!> @param face1 index to a face
!> @param face2 index to another face
!> @return SharedEdge = 0 if faces do not have a common edge, SharedEdge = 1 if faces have a common edge, 
!> SharedEdges > 1 is an error, probably with the mesh seed file.
function SharedEdge( self, face1, face2 )
	integer(kint) :: SharedEdge
	type(Faces), intent(in) :: self
	integer(kint), intent(in) :: face1, face2
	!
	integer(kint) :: i, j, shareCount
	
	SharedEdge = 0
	shareCount = 0
	do i = 1, 3
		do j = 1, 3
			if ( self%edges(i, face1) == self%edges(j, face2) ) then
				shareCount = shareCount + 1
				SharedEdge = self%edges(j,face2)
			endif
		enddo
	enddo
	if ( shareCount > 1 ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL, trim(logkey)//" SharedEdge ERROR at face ", face1)
		SharedEdge = 0
	endif
end function

!> @brief Returns the index to a vertex particle opposite to an edge in a triangular face.
!> @param self Target faces object
!> @param faceIndex target face
!> @param edgeIndex index to an edge from an edgesmodule::edges object
function ParticleOppositeTriEdge( self, faceIndex, edgeIndex )
	integer(kint) :: ParticleOppositeTriEdge
	type(Faces), intent(in) :: self
	integer(kint), intent(in) :: faceIndex
	integer(kint), intent(in) :: edgeIndex
	!
	integer(kint) :: i, foundEdge
	
	foundEdge = 0
	ParticleOppositeTriEdge = 0
	do i = 1, 3
		if ( self%edges(i, faceIndex) == edgeIndex ) then
			foundEdge = i
			exit
		endif
	enddo
	if ( foundEdge == 1 ) then
		ParticleOppositeTriEdge = self%vertices(3,faceIndex)
	elseif ( foundEdge == 2 ) then
		ParticleOppositeTriEdge = self%vertices(1,faceIndex)
	elseif ( foundEdge == 3 ) then
		ParticleOppositeTriEdge = self%vertices(2,faceIndex)
	else
		call LogMessage(log, ERROR_LOGGING_LEVEL, trim(logKey)//" ParticleOppositeTriEdge ERROR : ", "edge not found.")
	endif		
end function

!> @brief Outputs connectivity data to a legacy format .vtk ASCII file.
!> Note: This subroutine must be called in the correct order, after outputwritermodule::writevtkfileheader 
!> and particlesmodule::writevtkpoints.
!> @param self Target faces object
!> @param fileunit Integer unit associated with .vtk file
subroutine WriteFacesToVTKPolygons( self, fileunit )
	type(Faces), intent(in) :: self
	integer(kint), intent(in) :: fileunit
	!
	integer(kint) :: i, j, nCells, cellListSize, nVerts
	
	if ( self%faceKind == TRI_PANEL ) then
		nVerts = 3
	elseif ( self%faceKind == QUAD_PANEL ) then
		nVerts = 4
	endif
	nCells = nVerts * self%N_Active
	cellListSize = 4 * nCells
	
	write(fileunit, '(A,I8,A,I8)') "POLYGONS ", nCells, "   ", cellListSize
	do i = 1, self%N
		if ( .NOT. self%hasChildren(i) ) then
			do j = 1, nVerts
				write(fileunit,'(4I10)') 3, self%vertices(j,i)-1, self%vertices( mod(j,nVerts) + 1, i)-1, self%centerParticle(i)-1
			enddo
		endif
	enddo
end subroutine

!> @brief Outputs face areas data to a legacy format .vtk ASCII file.
!> Note: This subroutine must be called in the correct order, after all vtk point data has been written
!> by calls to fieldmodule::writefieldtovtkpointdata and before any other calls to fieldmodule::writefieldtovtkcelldata
!> @param self Target faces object
!> @param aParticles Particles object associated with this set of faces
!> @param fileunit Integer unit associated with .vtk file
subroutine WriteFaceAreaToVTKCellData(self, aParticles, fileunit )
	type(Faces), intent(in) :: self
	type(Particles), intent(in) :: aParticles
	integer(kint), intent(in) :: fileunit
	!
	integer(kint) :: i, j, nCells, nVerts
	
	if ( self%faceKind == TRI_PANEL ) then
		nVerts = 3
	elseif ( self%faceKind == QUAD_PANEL ) then
		nVerts = 4
	endif
	nCells = nVerts * self%N_Active
	
	write(fileunit,'(A,I8)') "CELL_DATA ", nCells
	write(fileunit,'(A)') "SCALARS faceArea double 1"
	write(fileunit,'(A)') "LOOKUP_TABLE default"
	do i = 1, self%N
		if ( .NOT. self%hasChildren(i) ) then
			do j = 1, nVerts
				write(fileunit,*) aParticles%area( self%centerParticle(i) )
			enddo
		endif
	enddo		
end subroutine

!> @brief Outputs data associated with a faces object, including all connectivity information, to a Matlab-readable .m file.
!> @param self Target faces object
!> @param fileunit Integer unit associated with a .m file
subroutine WriteFacesToMatlab( self, fileunit )
	type(Faces), intent(in) :: self
	integer(kint), intent(in) :: fileunit
	!
	integer(kint) :: i, j, nVerts
	
	if ( self%faceKind == TRI_PANEL ) then
		nVerts = 3
	elseif ( self%faceKind == QUAD_PANEL ) then
		nVerts = 4
	else
		nVerts = 0
	endif
	
	write(fileunit, '(A)', advance='NO') "faceVerts = ["
	do j = 1, nVerts - 1
		write(fileunit,'(I8,A)', advance = 'NO') self%vertices(j,1), ", "
	enddo
	write(fileunit,'(I8,A)') self%vertices(nVerts,1), "; ..."
	do i = 2, self%N-1
		do j = 1, nVerts - 1
			write(fileunit,'(I8,A)',advance = 'NO') self%vertices(j,i), ", "
		enddo
		write(fileunit,'(I8,A)' ) self%vertices(nVerts,i), "; ..."
	enddo
	do j = 1, nVerts - 1
		write(fileunit,'(I8,A)',advance='NO') self%vertices(j, self%N), ", "
	enddo
	write(fileunit,'(I8,A)' ) self%vertices(nVerts,self%N), "]; "

	write(fileunit, '(A)', advance='NO') "faceEdges = ["
	do j = 1, nVerts - 1
		write(fileunit,'(I8,A)', advance = 'NO') self%edges(j,1), ", "
	enddo
	write(fileunit,'(I8,A)') self%edges(nVerts,1), "; ..."
	do i = 2, self%N-1
		do j = 1, nVerts - 1
			write(fileunit,'(I8,A)',advance = 'NO') self%edges(j,i), ", "
		enddo
		write(fileunit,'(I8,A)' ) self%edges(nVerts,i), "; ..."
	enddo
	do j = 1, nVerts - 1
		write(fileunit,'(I8,A)',advance='NO') self%edges(j, self%N), ", "
	enddo
	write(fileunit,'(I8,A)' ) self%edges(nVerts,self%N), "]; "
	
	write(fileunit,'(A)',advance='NO') "faceHasChildren = ["
	if ( self%hasChildren(1) ) then
		write(fileunit,*) 1, ", ..."
	else
		write(fileunit,*) 0, ", ..."
	endif
	do i = 2, self%N-1
		if ( self%hasChildren(i) ) then
			write(fileunit,*) 1, ", ..."
		else
			write(fileunit,*) 0, ", ..."
		endif
	enddo
	if ( self%hasChildren(self%N) ) then
		write(fileunit,'(I4)',advance='NO') 1
	else
		write(fileunit,'(I4)',advance='NO') 0
	endif
	write(fileunit,'(A)') "];"

	write(fileunit,*) "faceCenterParticle = [", self%centerParticle(1), ", ..."
	do i = 2, self%N-1
		write(fileunit,*) self%centerParticle(i), ", ..."
	enddo
	write(fileunit,*) self%centerParticle(self%N), "]; "
end subroutine

!> @brief This is a primary subroutine used by any quadrilateral mesh (spherical or planar).
!> It divides a quadrilateral face into 4 children, adding a 4 leaves to the faces quadtree.
!> It creates new particles and new edges.  
!> It replaces divided edges with their appropriate children.
!> 
!> @param[inout] self Target faces object
!> @param[in] faceIndex index of face to be divided
!> @param[inout] aParticles Particles object associated with this set of faces.  
!> @param[inout] anEdges Edges object associated with this set of faces
subroutine DivideQuadFace( self, faceIndex, aParticles, anEdges )
	type(Faces), intent(inout) :: self
	integer(kint), intent(in) :: faceIndex
	type(Particles), intent(inout) :: aParticles
	type(Edges), intent(inout) :: anEdges
	!
	integer(kint) :: i, j, newFaceVerts(4,4), newFaceEdges(4,4)
	integer(kint) :: parentEdge, childEdge1, childEdge2
	real(kreal) :: quadCoords(3,4), lagQuadCoords(3,4), newFaceCenters(3,4), newFaceLagCenters(3,4)
	
	if ( self%faceKind /= QUAD_PANEL ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL, logkey, " DivideQuadFace ERROR : invalid faceKind.")
		return
	endif
	if ( self%N_Max < self%N + 4 ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL, logkey, " DivideQuadFace ERROR : not enough memory.")
		return
	endif
	
	newFaceVerts = 0
	newFaceEdges = 0
	
	!
	! connect parent vertices to child faces
	!
	do i = 1, 4
		newFaceVerts(i,i) = self%vertices(i,faceIndex)
	enddo
	!
	! loop over parent edges
	!
	do i = 1, 4
		parentEdge = self%edges(i, faceIndex)
		if ( anEdges%hasChildren(parentEdge) ) then
			!
			! parent edge was already divided by adjacent panel
			!
			childEdge1 = anEdges%child1(parentEdge)
			childEdge2 = anEdges%child2(parentEdge)
		else
			!
			! parent edge needs to be divided
			!
			childEdge1 = anEdges%N + 1
			childEdge2 = anEdges%N + 2
			call DivideEdge(anEdges, parentEdge, aParticles)
		endif
		
		!
		! connect child edges to new child faces
		!
		if ( positiveEdge(anEdges, faceIndex, parentEdge) ) then
			newFaceEdges(i, i ) = childEdge1
			anEdges%leftFace(childEdge1) = self%N + i
			
			newFaceEdges(i, mod(i,4) + 1) = childEdge2
			anEdges%leftFace(childEdge2) = self%N + mod(i,4) + 1
		else
			newFaceEdges(i,i) = childEdge2
			anEdges%rightFace(childEdge2) = self%N + i
			
			newFaceEdges(i, mod(i,4) + 1) = childEdge1
			anEdges%rightFace(childEdge1) = self%N + mod(i,4) + 1
		endif
		
		newFaceVerts( mod(i,4) + 1, i ) = anEdges%dest(childEdge1)
		newFaceVerts( i, mod(i,4) + 1 ) = anEdges%dest(childEdge1)
	enddo
	!
	! change parent center particle to vertex of new child panels
	!
	do i = 1, 4
		newFaceVerts( mod(i+1,4) + 1, i) = self%centerParticle(faceIndex)
	enddo
	
	!
	!	debugging : check vertex connectivity
	!
	do i = 1, 4
		do j = 1, 4
			if ( newFaceVerts(j,i) < 1 ) then
				write(logstring,*) " vertex connectivity ERROR at parent face ", faceIndex, ", child ", i, ", vertex ", j
				call LogMessage(log, ERROR_LOGGING_LEVEL,logkey//" DivideQuadFace :",logstring)
			endif
		enddo
	enddo
	
	!
	!	create new interior edges
	!
	call InsertEdge( anEdges, aParticles, newFaceVerts(2,1), newFaceVerts(3,1), self%N+1, self%N+2)
	newFaceEdges(2,1) = anEdges%N
	newFaceEdges(4,2) = anEdges%N
	
	call InsertEdge( anEdges, aParticles, newFaceVerts(1,3), newFaceVerts(4,3), self%N+4, self%N+3)
	newFaceEdges(4,3) = anEdges%N
	newFaceEdges(2,4) = anEdges%N
	
	call InsertEdge( anEdges, aParticles, newFaceVerts(3,2), newFaceVerts(4,2), self%N+2, self%N+3)
	newFaceEdges(3,2) = anEdges%N
	newFaceEdges(1,3) = anEdges%N
	
	call InsertEdge( anEdges, aParticles, newFaceVerts(2,4), newFaceVerts(1,4), self%N+1, self%N+4)
	newFaceEdges(1,4) = anEdges%N
	newFaceEdges(3,1) = anEdges%N
	
	!
	! debugging : check edge connectivity
	!
	do i = 1, 4
		do j = 1, 4
			if ( newFaceEdges(j,i) < 1 ) then
				write(logstring,*) " edge connectivity ERROR at parent face ", faceIndex, ", child ", i, ", vertex ", j
				call LogMessage(log, ERROR_LOGGING_LEVEL, logkey//" DivideQuadFace :",logstring)	
			endif
		enddo
	enddo
	
	!
	! create new center particles for child faces
	!
	quadCoords = 0.0_kreal
	lagQuadCoords = 0.0_kreal
	newFaceCenters = 0.0_kreal
	newFaceLagCenters = 0.0_kreal
	do i = 1, 4
		do j = 1, 4
			quadCoords(:,j) = PhysCoord(aParticles, newFaceVerts(j,i))
			lagQuadCoords(:,j) = LagCoord(aParticles, newFaceVerts(j,i))
		enddo
		if ( aParticles%geomKind == PLANAR_GEOM ) then
			newFaceCenters(1:2,i) = QuadCentroid( quadCoords(1:2,1), quadCoords(1:2,2), quadCoords(1:2,3), quadCoords(1:2,4))		
			newFaceLagCenters(1:2,i) = QuadCentroid( lagQuadCoords(1:2,1), lagQuadCoords(1:2,2), lagQuadCoords(1:2,3), lagQuadCoords(1:2,4))
		elseif ( aParticles%geomKind == SPHERE_GEOM ) then
			newFaceCenters(:,i) = SphereQuadCenter( quadCoords(:,1), quadCoords(:,2), quadCoords(:,3), quadCoords(:,4))
			newFaceLagCenters(:,i) = SphereQuadCenter( lagQuadCoords(:,1), lagQuadCoords(:,2), lagQuadCoords(:,3), lagQuadCoords(:,4))
		else
			call LogMessage(log,ERROR_LOGGING_LEVEL,logkey//" DivideQuadFace : ", "ERROR geomKind not implemented.")
			return
		endif
	enddo
	
	!
	!	create child faces
	!
	do i = 1, 4
		call InsertParticle( aParticles, newFaceCenters(:,i), newFaceLagCenters(:,i))
		call MakeParticleActive(aParticles, aParticles%N)
		self%centerParticle(self%N+i) = aParticles%N
		self%vertices(:,self%N+i) = newFaceVerts(:,i)
		self%edges(:,self%N+i) = newFaceEdges(:,i)
		self%children(i,faceIndex) = self%N + i
		self%parent(self%N+i) = faceIndex
		aParticles%area(self%centerParticle(self%N+i)) = QuadFaceArea(self, self%N+i, aParticles)
	enddo
	call MakeParticlePassive(aParticles, self%centerParticle(faceIndex) )
	aParticles%area(self%centerParticle(faceIndex)) = 0.0_kreal
	self%hasChildren(faceIndex) = .TRUE.
	self%N = self%N + 4
	self%N_Active = self%N_Active + 3
end subroutine


!> @brief This is a primary subroutine used by any triangular mesh (spherical or planar).
!> It divides a triangular face into 4 children, adding a 4 leaves to the faces quadtree.
!> It creates new particles and new edges.  
!> It replaces divided edges with their appropriate children.
!> 
!> @param[inout] self Target faces object
!> @param[in] faceIndex index of face to be divided
!> @param[inout] aParticles Particles object associated with this set of faces.  
!> @param[inout] anEdges Edges object associated with this set of faces
subroutine DivideTriFace( self, faceIndex, aParticles, anEdges )
	type(Faces), intent(inout) :: self
	integer(kint), intent(in) :: faceIndex
	type(Particles), intent(inout) :: aParticles
	type(Edges), intent(inout) :: anEdges
	!
	integer(kint) :: i, j, newFaceVerts(3,4), newFaceEdges(3,4)
	integer(kint) :: parentEdge, childEdge1, childEdge2
	real(kreal) :: triCoords(3,3), lagTriCoords(3,3), newFaceCenters(3,4), newFaceLagCenters(3,4)
	
	if ( self%faceKind /= TRI_PANEL ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL, logkey//" DivideTriFace ERROR : "," invalid faceKind.")
		return
	endif
	if ( self%N_Max < self%N + 4 ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL, logkey//" DivideTriFace ERROR : ", " not enough memory.")
		return
	endif
	
	newFaceVerts = 0
	newFaceEdges = 0
	!
	!	connect parent vertices to child faces
	!
	do i = 1, 3
		newFaceVerts(i,i) = self%vertices(i,faceIndex)
	enddo
	!
	!	loop over parent edges
	!
	do i = 1, 3
		parentEdge = self%edges(i,faceIndex)
		if ( anEdges%hasChildren(parentEdge) ) then
			!
			!	parent edge already divided by adjacent panel
			!
			childEdge1 = anEdges%child1(parentEdge)
			childEdge2 = anEdges%child2(parentEdge)
		else
			!
			!	divide parent edge
			!
			childEdge1 = anEdges%N + 1
			childEdge2 = anEdges%N + 2
			call DivideEdge(anEdges, parentEdge, aParticles)
		endif
		
		!
		!	connect child edges to new child faces
		!
		if ( positiveEdge( anEdges, faceIndex, parentEdge )  ) then
			newFaceEdges(i, i) = childEdge1
			anEdges%leftFace(childEdge1) = self%N + i
			
			newFaceEdges(i, mod(i,3)+1) = childEdge2
			anEdges%leftFace(childEdge2) = self%N + mod(i,3) + 1
		else
			newFaceEdges(i,i) = childEdge2
			anEdges%rightFace(childEdge2) = self%N + i
			
			newFaceEdges(i,mod(i,3)+1) = childEdge1
			anEdges%rightFace(childEdge1) = self%N + mod(i,3) + 1
		endif
		
		newFaceVerts( i, mod(i,3) + 1) = anEdges%dest(childEdge1)
		newFaceVerts( mod(i,3) + 1, i) = anEdges%dest(childEdge1)
	enddo
	
	newFaceVerts(:,4) = [ newFaceVerts(3,2), newFaceVerts(1,3), newFaceVerts(2,1) ]
	
	!
	! debugging : check vertex connectivity
	!
	do i = 1, 4
		do j = 1, 3
			if ( newFaceVerts(j,i) < 1 ) then
				write(logstring,*) " vertex connectivity ERROR at parent face ", faceIndex, ", child ", i, ", vertex ",j
				call LogMessage(log,ERROR_LOGGING_LEVEL,logkey//" DivideTriFace :", logString )
			endif
		enddo
	enddo
	
	!
	! create new interior edges
	!
	call InsertEdge( anEdges, aParticles, newFaceVerts(1,4), newFaceVerts(2,4), self%N+4, self%N+3)
	newFaceEdges(1,4) = anEdges%N
	newFaceEdges(1,3) = anEdges%N
	
	call InsertEdge( anEdges, aParticles, newFaceVerts(2,4), newFaceVerts(3,4), self%N+4, self%N+1)
	newFaceEdges(2,4) = anEdges%N
	newFaceEdges(2,1) = anEdges%N
		
	call InsertEdge( anEdges, aParticles, newFaceVerts(3,4), newFaceVerts(1,4), self%N+4, self%N+2)
	newFaceEdges(3,4) = anEdges%N
	newFaceEdges(3,2) = anEdges%N
	
	!
	! debugging : check edge connectivity
	!
	do i = 1, 4
		do j = 1, 3
			if ( newFaceEdges(j,i) < 1 ) then
				write(logstring,*) " edge connectivity ERROR at parent face ", faceIndex, ", child ", i, ", vertex ", j
				call LogMessage(log, ERROR_LOGGING_LEVEL, logkey//" DivideTriFace :", logString)
			endif
		enddo
	enddo
	
	!
	! create new center particles for child faces 1:3
	!
	triCoords = 0.0_kreal
	lagTriCoords = 0.0_kreal
	newFaceCenters = 0.0_kreal
	newFaceLagCenters = 0.0_kreal
	do i = 1, 4
		do j = 1, 3
			triCoords(:,j) = PhysCoord(aParticles, newFaceVerts(j,i))
			lagTriCoords(:,j)=LagCoord(aParticles, newFaceVerts(j,i))
		enddo
		if ( aParticles%geomKind == PLANAR_GEOM ) then
			newFaceCenters(1:2,i) = TriCentroid( triCoords(1:2,1), triCoords(1:2,2), triCoords(1:2,3) )
			newFaceLagCenters(1:2,i) = TriCentroid( lagTriCoords(1:2,1), lagTriCoords(1:2,2), lagTriCoords(1:2,3))
		elseif ( aParticles%geomKind == SPHERE_GEOM ) then
			newFaceCenters(:,i) = SphereTriCenter( triCoords(:,1), triCoords(:,2), triCoords(:,3))
			newFaceLagCenters(:,i) = SphereTriCenter(lagTriCoords(:,1),lagTriCoords(:,2), lagTriCoords(:,3))
		else
			call LogMessage(log, ERROR_LOGGING_LEVEL, logkey//" DivideTriFace ERROR : "," geomkind not implemented.")
			return
		endif
	enddo
	
	!
	! create child faces
	!
	do i = 1, 3
		call InsertParticle( aParticles, newFaceCenters(:,i), newFaceLagCenters(:,i))
		call MakeParticleActive( aParticles, aParticles%N)
		self%centerParticle( self%N + i) = aParticles%N
		self%vertices(:,self%N+i) = newFaceVerts(:,i)
		self%edges(:,self%N+i) = newFaceEdges(:,i)
		self%children(i,faceIndex) = self%N+i
		self%parent(self%N+i) = faceIndex
		aParticles%area(self%centerParticle(self%N+i)) = TriFaceArea(self, self%N+i, aParticles)
	enddo
	!	child 4 is special case : reuse parent face center particle
	self%centerParticle( self%N + 4 ) = self%centerParticle(faceIndex)
	self%vertices(:,self%N+4) = newFaceVerts(:,4)
	self%edges(:,self%N + 4) = newFaceEdges(:,4)
	self%children(4,faceIndex) = self%N + 4
	self%parent(self%N+4) = faceIndex
	aParticles%area(self%centerParticle(self%N+4)) = TriFaceArea(self, self%N+4, aParticles)
	
	self%hasChildren(faceIndex) = .TRUE.
	self%N = self%N + 4
	self%N_Active = self%N_Active + 3	
end subroutine


!> @brief Computes the area of a polyhedral face
!> @param self Target faces object
!> @param index Index of face whose area is needed
!> @param aParticles Particles object  associated with this set of faces
!> @param anEdges Edges object associated with this set of faces
function FaceArea( self, index, aParticles, anEdges )
	real(kreal) :: FaceArea
	type(Faces), intent(in) :: self
	integer(kint), intent(in) :: index
	type(Particles), intent(in) :: aParticles
	type(Edges), intent(in) :: anEdges
	!
	type(STDIntVector) :: leafEdges
	integer(kint) :: i, nEdges
	
	nEdges = self%faceKind
	FaceArea = 0.0_kreal
	do i = 1, nEdges
		call GetLeafEdgesFromParent( anEdges, self%edges(i, index), leafEdges )
		FaceArea = FaceArea + AreaFromLeafEdges( anEdges, aParticles, self%centerParticle(index), &
												 leafEdges%integers(1:leafEdges%N), leafEdges%N)
	enddo
end function

!> @brief Inserts a new face into the faces data structure.
!> Increases faces.N by 1. 
!> @param self Target faces object
!> @param centerParticle index to a the new face's center particle in a particlesmodule::particles object
!> @param vertIndices indices to the new face's vertices in a particlesmodule::particles object
!> @param edgeIndices indices to the new face's edges in an edgesmodule::edges object.
subroutine InsertFace( self, centerParticle, vertIndices, edgeIndices )
	type(Faces), intent(inout) :: self
	integer(kint), intent(in) :: centerParticle
	integer(kint), intent(in) :: vertIndices(:)
	integer(kint), intent(in) :: edgeIndices(:)
	
	if ( self%N+1 > self%N_Max ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL, logkey," InsertFace ERROR : not enough memory.")
		return
	endif
	if ( size(vertIndices) /= size(edgeIndices) .OR. size(vertIndices) /= size(self%vertices,1) ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL, logkey, " InsertFace ERROR : vertex and edge size mismatch.")
		return
	endif
	self%centerParticle(self%N+1) = centerParticle
	self%vertices(:,self%N+1) = vertIndices
	self%edges(:,self%N+1) = edgeIndices
	self%N = self%N + 1
end subroutine

!> @brief Determines the area of a quadrilateral face.  Assumes the face has exactly 4 edges 
!> (Does not account for adjacent face's having different levels of refinement).
!> 
!> @param self Target faces object
!> @param index Index to face whose area is needed
!> @param aParticles Particles object associated with this set of faces
function QuadFaceArea( self, index, aParticles )
	real(kreal) :: QuadFaceArea
	type(Faces), intent(in) :: self
	integer(kint), intent(in) :: index
	type(Particles), intent(in) :: aParticles
	!
	integer(kint) :: i
	real(kreal) :: centerVec(3), v1(3), v2(3)
	
	QuadFaceArea = 0.0_kreal
	centerVec = PhysCoord(aParticles, self%centerParticle(index))
	
	if ( aParticles%geomKind == PLANAR_GEOM ) then
		do i = 1, 4
			v1 = PhysCoord(aParticles, self%vertices(i, index))
			v2 = PhysCoord(aParticles, self%vertices(mod(i,4) + 1, index))
			QuadFaceArea = QuadFaceArea + TriArea(v1(1:2), centerVec(1:2), v2(1:2) )		
		enddo
	elseif ( aParticles%geomKind == SPHERE_GEOM ) then
		do i = 1, 4
			v1 = PhysCoord(aParticles, self%vertices(i, index))
			v2 = PhysCoord(aParticles, self%vertices(mod(i,4) + 1, index))
			QuadFaceArea = QuadFaceArea + SphereTriArea( v1, centerVec, v2 )
		enddo	
	endif
end function

!> @brief Determines the area of a triangular face.  Assumes the face has exactly 3 edges 
!> (Does not account for adjacent face's having different levels of refinement).
!> 
!> @param self Target faces object
!> @param index Index to face whose area is needed
!> @param aParticles Particles object associated with this set of faces
function TriFaceArea( self, index, aParticles )
	real(kreal) :: TriFaceArea
	type(Faces), intent(in) :: self
	integer(kint), intent(in) :: index
	type(Particles), intent(in) :: aParticles
	!
	integer(kint) :: i
	real(kreal) :: centerVec(3), v1(3), v2(3)
	
	TriFaceArea = 0.0_kreal
	centerVec = PhysCoord(aParticles, self%centerParticle(index))
	
	if ( aParticles%geomKind == PLANAR_GEOM ) then
		do i = 1, 3
			v1 = PhysCoord(aParticles, self%vertices(i, index))
			v2 = PhysCoord(aParticles, self%vertices(mod(i,3) + 1, index))
			TriFaceArea = TriFaceArea + TriArea(v1(1:2), centerVec(1:2), v2(1:2) )		
		enddo
	elseif ( aParticles%geomKind == SPHERE_GEOM ) then
		do i = 1, 3
			v1 = PhysCoord(aParticles, self%vertices(i, index))
			v2 = PhysCoord(aParticles, self%vertices(mod(i,3) + 1, index))
			TriFaceArea = TriFaceArea + SphereTriArea( v1, centerVec, v2 )
		enddo	
	endif
end function

!> @brief Computes the centroid (in physical space) of a face based on its vertices.
!> @param[in] self Target faces object
!> @param[in] index Index of face whose centroid is needed
!> @param[in] aParticles Particles object associated with this set of faces
!> @return Position vector of face centroid
function FaceCentroid(self, index, aParticles )
	real(kreal) :: FaceCentroid(3)
	type(Faces), intent(in) :: self
	integer(kint), intent(in) :: index
	type(Particles), intent(in) :: aParticles
	!
	integer(kint) :: nVerts, i
	
	!call LogMessage(log,DEBUG_LOGGING_LEVEL,trim(logKey)//" FaceCentroid : ", "entering.")
	
	FaceCentroid = 0.0_kreal
	nVerts = 0
	if ( self%faceKind == TRI_PANEL ) then
		nVerts = 3
	elseif( self%faceKind == QUAD_PANEL) then
		nVerts = 4
	endif
	do i = 1, nVerts
		FaceCentroid(1) = FaceCentroid(1) + aParticles%x( self%vertices(i,index) )
		FaceCentroid(2) = FaceCentroid(2) + aParticles%y( self%vertices(i,index) )
		if ( aParticles%geomKind /= PLANAR_GEOM ) &
			FaceCentroid(3) = FaceCentroid(3) + aParticles%z(self%vertices(i,index))
	enddo
	
	FaceCentroid = FaceCentroid / real(nVerts,kreal)
	
	if ( aParticles%geomKind == SPHERE_GEOM ) then
		FaceCentroid = FaceCentroid / sqrt(sum(FaceCentroid*FaceCentroid)) * SphereRadius
	endif
end function

!
!----------------
! Private methods
!----------------
!

!> @brief Initializes a logger for the Faces module
!> 
!> Output is controlled both by message priority and by MPI Rank
!> @param[out] aLog Target Logger object
!> @param[in] rank Rank of this processor
subroutine InitLogger(aLog,rank)
! Initialize a logger for this module and processor
	type(Logger), intent(out) :: aLog
	integer(kint), intent(in) :: rank
	write(logKey,'(A,A,I0.3,A)') trim(logKey),'_',rank,' : '
	call New(aLog,logLevel)
	logInit = .TRUE.
end subroutine

!> @}
end module