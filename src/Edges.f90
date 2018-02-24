module EdgesModule
!> @file Edges.f90
!> Provides a primitive data structure and methods for creating edges of polyhedral meshes.
!> @author Peter Bosler, Sandia National Laboratories Center for Computing Research
!>
!>
!> @defgroup Edges Edges
!> @brief Edges of polyhedral meshes connect to vertices and faces.
!>
!> Edges know the indices of their origin and destination particles (in a @ref Particles object), and the indices 
!> of their left face and right face (in a @ref Faces object).
!> 
!> The edges data structures is a "structure of arrays," so that all information about edge i is located at index i in the 
!> relevant array.  
!> For example, the index to the particle at the origin of edge i is `anEdges%%orig(i)` , and its destination particle is
!> `anEdges%%dest(i)`. 
!> An edge's left and right faces are accessed similarly, `anEdges%%leftFace(i)` gives the index of the face (in a @ref Faces object)
!> to the left of edge i.
!>
!> In addition the edges are organized into a binary tree, defined by the `hasChildren`, `children`, and `parent` arrays,
!> to facilitate faster searching through the data structure.
!> 
!> @{
use NumberKindsModule
use LoggerModule
use STDIntVectorModule
use ParticlesModule
use SphereGeomModule, only : crossProduct, SphereMidpoint, SphereTriArea, SphereDistance, ChordDistance
use PlaneGeomModule, only : Midpoint, TriArea
use ParticlesModule

implicit none

private
public Edges
public New, Delete, Copy
public InsertEdge, DivideEdge
public RecordIncidentEdgeAtParticles
public positiveEdge
public onBoundary
public ReplaceIncidentEdgeWithChild
public GetLeafEdgesFromParent, AreaFromLeafEdges
public EdgeLength, MaxEdgeLength, MinEdgeLength, AvgEdgeLength
public LogStats, PrintDebugInfo
public WriteEdgesToMatlab
public CountParents

!> @brief Edges know the indices (to @ref Particles) of their origin and destination, and the indices (to @ref Faces)
!> of their left face and right face.
type Edges
    integer(kint), allocatable :: orig(:) !< Integer array containing indices of particlesmodule::particles     
    integer(kint), allocatable :: dest(:) !< Integer array containing indices of particlesmodule::particles 
    integer(kint), allocatable :: rightFace(:) !< Integer array containing indices of facesmodule::faces
    integer(kint), allocatable :: leftFace(:) !< Integer array containing indices of facesmodule::faces
    integer(kint), allocatable :: child1(:) !< Integer array containing indices to child edges with same origin vertex as parent
    integer(kint), allocatable :: child2(:) !< Integer array containing indices to child edges with same destination vertex as parent
    logical(klog), allocatable :: hasChildren(:) !< hasChildren(i) is .TRUE. if edge i has been divided
    integer(kint), allocatable :: parent(:) !< Null for root edges.  Integer pointers into edges for divided edges
    integer(kint) :: N = 0 !< Number of edges currently in use
    integer(kint) :: N_Max = 0 !< Max number of edges allowed in memory
    
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
    module procedure countParentEdges
end interface 


interface DivideEdge
    module procedure divideLinearEdge
end interface

!
!----------------
! Logging
!----------------
!
logical(klog), save :: logInit = .FALSE.
type(Logger) :: log
character(len=28), save :: logKey = 'Edges'
integer(kint), parameter :: logLevel = DEBUG_LOGGING_LEVEL

contains
!
!----------------
! Public methods
!----------------
!
!> @brief Allocates memory for a new edges object.  All values are zeroed and must be initialized separately.
!> @param self Target edge object
!> @param nMax max number of edges allowed in memory
subroutine NewPrivate(self, nMax )
    type(Edges), intent(out) :: self
    integer(kint), intent(in) :: nMax
    
    if ( .NOT. logInit ) call InitLogger(log, procRank )
    
    if ( nMax <= 0 ) then
        call LogMessage(log, ERROR_LOGGING_LEVEL, logkey, "New Edges ERROR : invalid nMax.")
        return
    endif
    
    self%N_Max = nMax
    self%N = 0
    
    allocate(self%orig(nMax))
    self%orig = 0
    allocate(self%dest(nMax))
    self%dest = 0
    allocate(self%leftFace(nMax))
    self%leftFace = 0
    allocate(self%rightFace(nMax))
    self%rightFace = 0
    allocate(self%hasChildren(nMax))
    self%hasChildren = .FALSE.
    allocate(self%child1(nMax))
    self%child1 = 0
    allocate(self%child2(nMax))
    self%child2 = 0
    allocate(self%parent(nMax))
    self%parent = 0
end subroutine

!> @brief Deletes and frees memory associated with an edges object
!> @param self Target edge object
subroutine deletePrivate(self)
    type(Edges), intent(inout) :: self
    if ( allocated(self%orig)) deallocate(self%orig)
    if ( allocated(self%dest)) deallocate(self%dest)
    if ( allocated(self%leftFace)) deallocate(self%leftFace)
    if ( allocated(self%rightFace)) deallocate(self%rightFace)
    if ( allocated(self%hasChildren)) deallocate(self%hasChildren)
    if ( allocated(self%child1)) deallocate(self%child1)
    if ( allocated(self%child2)) deallocate(self%child2)
    if ( allocated(self%parent)) deallocate(self%parent)
end subroutine

!> @brief Performs a deep copy of one edges object into another.  They must have been allocated the same to avoid memory errors.
!> @param self Target edge object
!> @param other Source edge object
subroutine copyPrivate( self, other )
    type(Edges), intent(inout) :: self
    type(Edges), intent(in) :: other
    !
    integer(kint) :: j
    
    if ( self%N_Max < other%N ) then
        call LogMessage( log, ERROR_LOGGING_LEVEL, logkey, "CopyEdges ERROR : not enough memory.")
        return
    endif
    
    do j = 1, other%N
        self%orig(j) = other%orig(j)
        self%dest(j) = other%dest(j)
        self%leftFace(j) = other%leftFace(j)
        self%rightFace(j) = other%rightFace(j)
        self%hasChildren(j) = other%hasChildren(j)
        self%child1(j) = other%child1(j)
        self%child2(j) = other%child2(j)
        self%parent(j) = other%parent(j)
    enddo
    self%N = other%N
end subroutine

!> @brief Log summarizing statistics about an edges object
!> @param self Target edge object
!> @param aLog Target loggermodule::logger object for output
subroutine LogStatsPrivate( self, aLog )
    type(Edges), intent(in) :: self
    type(Logger), intent(inout) :: aLog
    call LogMessage(aLog, TRACE_LOGGING_LEVEL, logkey, " Edges Stats : ")
    call StartSection(aLog)
    call LogMessage(aLog, TRACE_LOGGING_LEVEL, "edges.N = ", self%N )
    call LogMessage(aLog, TRACE_LOGGING_LEVEL, "edges.N_Max = ", self%N_Max)
    call LogMessage(aLog, TRACE_LOGGING_LEVEL, "n divided edges = ", count(self%hasChildren) )
    call LogMessage(aLog, TRACE_LOGGING_LEVEL, "n leaf edges = ", self%N - count(self%hasChildren))
    call EndSection(aLog)
end subroutine

!> @brief Inserts an edge to the end of an edges object.  Increase edges.N by one.
!> @param self Target edge object
!> @param aParticles particlesmdoule::particles object associated with this set of edges
!> @param origIndex index to origin particle
!> @param destIndex index to destination particle
!> @param leftFace index to left face of new edge in faces object
!> @param rightFace index to right face of new edge in faces object
subroutine InsertEdge( self, aParticles, origIndex, destIndex, leftFace, rightFace,  interiorIndices )
    type(Edges), intent(inout) :: self
    type(Particles), intent(inout) :: aParticles
    integer(kint), intent(in) :: origIndex, destIndex, leftFace, rightFace
    integer(kint), dimension(:), intent(in), optional :: interiorIndices
    !
    integer(kint) :: n
    
    if ( self%N >= self%N_Max ) then
        call LogMessage(log, ERROR_LOGGING_LEVEL, logKey, " InsertEdge : out of memory. ")
        return
    endif
    
    n = self%N
    
    self%orig( n + 1 ) = origIndex
    self%dest( n + 1 ) = destIndex
    self%leftFace( n + 1) = leftFace
    self%rightFace(n + 1) = rightFace
        
    call RecordIncidentEdgeAtParticles( self, n + 1,    aParticles )
    
    self%N = n + 1
end subroutine

!> @brief Counts the number of levels in the edges binary tree above the edge at index
!> @param self Target edge object
!> @param index index of target edge
function countParentEdges( self, index )
    integer(kint) :: countParentEdges
    type(Edges), intent(in) :: self
    integer(kint), intent(in) :: index
    !
    logical(klog) :: keepGoing
    integer(kint) :: parentIndex
    
    countParentEdges = 0
    keepGoing = ( self%parent(index) > 0 )
    parentIndex = self%parent(index)
    do while ( keepGoing )
        countParentEdges = countParentEdges + 1
        parentIndex = self%parent( parentIndex )
        keepGoing = ( self%parent(parentIndex) > 0 )
    enddo
end function 

!> @brief Computes the maximum edge length amongst all leaf edges
!> @param self Target edge object
!> @param aParticles particles object associated with this set of edges
!> @return Maximum edge length
function MaxEdgeLength( self, aParticles )
    real(kreal) :: MaxEdgeLength
    type(Edges), intent(in) :: self
    type(Particles), intent(in) :: aParticles
    !
    integer(kint) :: i
    real(kreal) :: testLength
    
    MaxEdgeLength = 0.0_kreal
    do i = 1, self%N
        if ( .NOT. self%hasChildren(i) ) then
            testLength = EdgeLength(self, i, aParticles)
            if ( testLength > MaxEdgeLength ) MaxEdgeLength = testLength
        endif
    enddo
end function

!> @brief Computes the length of the edge at edgeIndex
!> @param self Target edge object
!> @param edgeIndex location of edge whose length is needed
!> @param aParticles particles object associated with this set of edges
!> @return Length of edge
function EdgeLength(self, edgeIndex, aParticles )
    real(kreal) :: EdgeLength
    integer(kint), intent(in) :: edgeIndex
    type(Edges), intent(in) :: self
    type(Particles), intent(in) :: aParticles
    !
    real(kreal) :: v0(3), v1(3)
    v0 = PhysCoord(aParticles, self%orig(edgeIndex))
    v1 = PhysCoord(aParticles, self%dest(edgeIndex))
    EdgeLength = 0.0_kreal
    if ( aParticles%geomKind == SPHERE_GEOM ) then
        EdgeLength = SphereDistance( v0, v1 )
    else
        EdgeLength = ChordDistance(v0, v1)
    endif
end function

!> @brief Computes the minimum edge length amongst all leaf edges
!> @param self Target edge object
!> @param aParticles particles object associated with this set of edges
!> @return Minimum edge length
function MinEdgeLength( self, aParticles)
    real(kreal) :: MinEdgeLength
    type(Edges), intent(in) :: self
    type(Particles), intent(in) :: aParticles
    !
    integer(kint) :: i
    real(kreal) :: testLen
    
    MinEdgeLength = EdgeLength(self, 1, aParticles)
    do i = 2, self%N
        if ( .NOT. self%hasChildren(i) ) then
            testLen = EdgeLength(self, i, aParticles)
            if ( testLen < MinEdgeLength) MinEdgeLength = testLen
        endif
    enddo
end function

!> @brief Computes the average edge length amongst all leaf edges
!> @param self Target edge object
!> @param aParticles particles object associated with this set of edges
!> @return average edge length
function AvgEdgeLength(self, aParticles)
    real(kreal) :: AvgEdgeLength
    type(Edges), intent(in) :: self
    type(Particles), intent(in) :: aParticles
    !
    integer(kint) :: i, nLeaves
    real(kreal) :: sum
    
    nLeaves = self%N - count(self%hasChildren(1:self%N))
    sum = 0.0_kreal
    do i = 1, self%N
        if ( .NOT. self%hasChildren(i) ) then
            sum = sum + EdgeLength(self, i, aParticles)
        endif
    enddo
    AvgEdgeLength = sum / real(nLeaves, kreal)
end function

!> @brief Writes the edges data structure to a .m file for later reading by Matlab.  For debugging.
!> @param self Target edge object
!> @param fileunit integer unit of .m file
subroutine WriteEdgesToMatlab( self, fileunit )
    type(Edges), intent(in) :: self
    integer(kint), intent(in) :: fileunit
    !
    integer(kint) :: i
    write(fileunit,*) "edgeVerts = [ ", self%orig(1), ", ", self%dest(1), "; ..."
    do i = 2, self%N-1
        write(fileunit, * ) self%orig(i), ", ", self%dest(i), "; ..."
    enddo
    write(fileunit, *) self%orig(self%N), ", ", self%dest(self%N), "]; "
    write(fileunit,'(A)',advance='NO') "edgeHasChildren = ["
    do i = 1, self%N - 1
        if ( self%hasChildren(i) ) then
            write(fileunit,*) 1, ", ..."
        else
            write(fileunit,*) 0, ", ..."
        endif
    enddo
    if ( self%hasChildren(self%N)) then
        write(fileunit,'(I4)', advance='NO') 1
    else
        write(fileunit,'(I4)', advance='NO') 0
    endif
    write(fileunit,'(A)') "];"
end subroutine

!> @brief Writes edges object values to the console.  For debugging.
!> @param self Target edge object
subroutine PrintDebugPrivate( self ) 
    type(Edges), intent(in) :: self
    integer(kint) :: i
    print *, "Edges DEBUG info : "
    print *, "edges.N = ", self%N
    print *, "edges.N_Max = ", self%N_Max
    print *, "edge records : "
    do i = 1, self%N_Max
        print *, self%orig(i), self%dest(i), self%leftFace(i), self%rightFace(i)
    enddo
    print *, "edge tree : "
    do i = 1, self%N_Max
        print *, self%hasChildren(i), self%child1(i), self%child2(i), self%parent(i)
    enddo
end subroutine

!> @brief Records the index of an incident edge at its origin and destination particles.  
!> Enables dual mesh functionality.
!> @todo This doesn't work.
!>
!> @param self Target edge object
!> @param edgeIndex target edge
!> @param aParticles particles object associated with this set of edges
subroutine RecordIncidentEdgeAtParticles( self, edgeIndex, aParticles )
    type(Edges), intent(in) :: self
    integer(kint), intent(in) :: edgeIndex
    type(Particles), intent(inout) :: aParticles
    !
    logical :: duplicateEdge
    integer(kint) :: j, origParticle, destParticle
    real(kreal) :: angleVal
    
    if ( self%hasChildren(edgeIndex) ) then
        return
    else
        origParticle = self%orig(edgeIndex)
        destParticle = self%dest(edgeIndex)
        !
        ! origin vertex
        !
        duplicateEdge = .FALSE.
        do j = 1, aParticles%nEdges( origParticle )
            if ( aParticles%incidentEdges( j, origParticle ) == edgeIndex ) duplicateEdge = .TRUE.
        enddo
        if ( .NOT. duplicateEdge ) then
            if ( aParticles%nEdges( origParticle ) >= MAX_VERTEX_DEGREE ) then
                call LogMessage(log,ERROR_LOGGING_LEVEL,logkey, " recordEdgeAtOrigin : out of memory.")
                return
            endif
            angleVal = edgeAngleAtOrig( self,  edgeIndex, aParticles)
            aParticles%incidentEdges( aParticles%nEdges( origParticle ) + 1, origParticle ) = edgeIndex
            aParticles%incidentAngles(aParticles%nEdges( origParticle ) + 1, origParticle ) = angleVal
            aParticles%nEdges(origParticle) = aParticles%nEdges(origParticle) + 1
        endif
        
        !
        ! destination vertex
        !
        duplicateEdge = .FALSE.
        do j = 1, aParticles%nEdges( self%dest( edgeIndex ))
            if ( aParticles%incidentEdges( j, self%dest(edgeIndex)) == edgeIndex ) duplicateEdge = .TRUE.
        enddo
        if ( .NOT. duplicateEdge ) then
            if ( aParticles%nEdges( destParticle) >= MAX_VERTEX_DEGREE ) then
                call LogMessage(log,ERROR_LOGGING_LEVEL,logkey, " recordEdgeAtDestination : out of memory.")
                return
            endif
            angleVal = edgeAngleAtDest( self, edgeIndex, aParticles )
            aParticles%incidentEdges( aParticles%nEdges(destParticle) + 1, destParticle ) = edgeIndex
            aParticles%incidentAngles(aParticles%nEdges(destParticle) + 1, destParticle ) = angleVal
            aParticles%nEdges(destParticle) = aParticles%nEdges(destParticle) + 1
        endif
    endif
end subroutine

!> @brief Determines the angle relative to the positive x-direction (for planar geometry) or the angle relative to the 
!> first edge listed at the origin particle (for spherical geometry).  
!> Used for dual-mesh functionality.
!> @todo This doesn't work.
!> 
!> @param self Target edge object
!> @param edgeIndex target edge
!> @param aParticles particles object associated with this set of edges
function edgeAngleAtOrig( self, edgeIndex, aParticles )
    real(kreal) :: edgeAngleAtOrig
    type(Edges), intent(in) :: self
    integer(kint), intent(in) :: edgeIndex
    type(Particles), intent(in) :: aParticles
    !
    real(kreal) :: edge1Vec(3), newEdgeVec(3), cp(3)
    
    edgeAngleAtOrig = 0.0_kreal
    if ( aParticles%geomKind == PLANAR_GEOM ) then
        !
        !   define angle relative to positive real axis
        !
        edgeAngleAtOrig = atan2( aParticles%y( self%dest(edgeIndex)) - aParticles%y( self%orig(edgeIndex)), &
                                 aParticles%x( self%dest(edgeIndex)) - aParticles%x( self%orig(edgeIndex)) )
    elseif ( aParticles%geomKind == SPHERE_GEOM ) then
        !
        !   define angle relative to first edge at particle
        !
        if ( aParticles%nEdges( self%orig(edgeIndex) ) > 1 ) then
            if ( self%orig(aParticles%incidentEdges(1,self%orig(edgeIndex))) == self%orig(edgeIndex) ) then
                edge1Vec = edgeVector( self, aParticles%incidentEdges(1,self%orig(edgeIndex)), aParticles )
            elseif ( self%dest(aParticles%incidentEdges(1,self%orig(edgeIndex))) == self%orig(edgeIndex) ) then
                edge1Vec = -edgeVector(self, aParticles%incidentEdges(1,self%orig(edgeIndex)), aParticles)
            else
                call LogMessage(log, ERROR_LOGGING_LEVEL, "edgeAngleAtOrig : ", "connectivity error.")
                return
            endif
            newEdgeVec = edgeVector( self, edgeIndex, aParticles)
            edge1Vec = edge1Vec/sqrt(sum( edge1Vec*edge1Vec))
            newEdgeVec = newEdgeVec/sqrt(sum( newEdgeVec*newEdgeVec))
            
            cp = crossProduct(edge1Vec, newEdgeVec)
            edgeAngleAtOrig = atan2( sqrt(sum(cp*cp)), sum(edge1Vec*newEdgeVec))
        endif
    else
        call LogMessage(log, WARNING_LOGGING_LEVEL, logkey, " geomKind not implemented yet.")
    endif
end function

!> @brief Determines the angle of edge at edgeIndex relative to the positive x-direction (for planar geometry) or the angle relative to the 
!> first edge listed at the destination particle (for spherical geometry).  
!> Used for dual-mesh functionality.
!> @todo This doesn't work.
!> 
!> @param self Target edge object
!> @param edgeIndex target edge
!> @param aParticles particles object associated with this set of edges
function edgeAngleAtDest( self, edgeIndex, aParticles )
    real(kreal) :: edgeAngleAtDest
    type(Edges), intent(in) :: self
    integer(kint), intent(in) :: edgeIndex
    type(Particles), intent(in) :: aParticles
    !
    real(kreal) :: edge1Vec(3), newEdgeVec(3), cp(3)
    
    edgeAngleAtDest = 0.0_kreal
    if ( aParticles%geomKind == PLANAR_GEOM ) then
        !
        !   define angle relative to positive real axis
        !
        edgeAngleAtDest = atan2( aParticles%y( self%orig(edgeIndex)) - aParticles%y( self%dest(edgeIndex)), &
                                 aParticles%x( self%orig(edgeIndex)) - aParticles%x( self%dest(edgeIndex)) )
    elseif ( aParticles%geomKind == SPHERE_GEOM ) then
        !
        !   define angle relative to first edge at particle
        !
        if ( aParticles%nEdges( self%dest(edgeIndex) ) > 1 ) then
            if ( self%orig(aParticles%incidentEdges(1,self%dest(edgeIndex))) == self%dest(edgeIndex) ) then
                edge1Vec = edgeVector(self, aParticles%incidentEdges(1,self%dest(edgeIndex)), aParticles)
            elseif ( self%dest(aParticles%incidentEdges(1,self%dest(edgeIndex))) == self%dest(edgeIndex)) then
                edge1Vec = -edgeVector(self, aParticles%incidentEdges(1,self%dest(edgeIndex)), aParticles)
            else
                call LogMessage(log, ERROR_LOGGING_LEVEL, "edgeAngleAtDest : ", "connectivity error.")
                return
            endif
        endif
        newEdgeVec = -edgeVector(self,edgeIndex,aParticles)
        
        edge1Vec = edge1Vec/sqrt(sum( edge1Vec*edge1Vec))
        newEdgeVec = newEdgeVec/sqrt(sum( newEdgeVec*newEdgeVec))
        
        cp = crossProduct(edge1Vec, newEdgeVec)
        edgeAngleAtDest = atan2( sqrt(sum(cp*cp)), sum(edge1Vec*newEdgeVec))
    else
        call LogMessage(log, WARNING_LOGGING_LEVEL, logkey, " geomKind not implemented yet.")
    endif
end function

!> @brief Returns a vector pointing from an edge's origin particle to its destination particle
!> @param self Target edge object
!> @param edgeIndex target edge
!> @param aParticles particles object associated with this set of edges
function edgeVector( self,  edgeIndex, aParticles )
    real(kreal), dimension(3) :: edgeVector
    type(Edges), intent(in) :: self
    integer(kint), intent(in) :: edgeIndex
    type(Particles), intent(in) :: aParticles
    edgeVector(1) = aParticles%x( self%dest(edgeIndex)) - aParticles%x( self%orig(edgeIndex))
    edgeVector(2) = aParticles%y( self%dest(edgeIndex)) - aParticles%y( self%orig(edgeIndex))
    if ( allocated( aParticles%z ) ) then
        edgeVector(3) = aParticles%z( self%dest(edgeIndex)) - aParticles%z( self%orig(edgeIndex))
    else
        edgeVector(3) = 0.0_kreal
    endif
end function

!> @brief Divides a parent edge, creating two child edges with the same orientation.
!> @param self Target edge object
!> @param edgeIndex target edge
!> @param aParticles particles object associated with this set of edges
subroutine divideLinearEdge( self, edgeIndex, aParticles )
    type(Edges), intent(inout) :: self
    integer(kint), intent(in) :: edgeIndex
    type(Particles), intent(inout) :: aParticles
    !
    real(kreal) :: midPt(3), lagMidPt(3), v0(3), v1(3), lV0(3), lV1(3)
    integer(kint) :: pInsertIndex
    
    if ( self%N + 2 > self%N_Max ) then
        call LogMessage(log, ERROR_LOGGING_LEVEL, logkey, " divideLinearEdge : out of memory.")
        return
    endif
    
    v0 = PhysCoord(aParticles, self%orig(edgeIndex))
    lV0 = LagCoord(aParticles, self%orig(edgeIndex))
    v1 = PhysCoord(aParticles, self%dest(edgeIndex))
    lV1 = LagCoord(aParticles, self%dest(edgeIndex))
        
    if ( aParticles%geomKind == SPHERE_GEOM ) then
        midPt = SphereMidpoint( v0, v1 )
        lagMidPt = SphereMidpoint( lv0, lv1 )
    else
        midPt = 0.5_kreal * (v0 + v1)
        lagMidPt = 0.5_kreal * ( lv0 + lv1)
    endif
    
    pInsertIndex = aParticles%N+1
    call InsertParticle( aParticles, midPt, lagMidPt )
    
    self%hasChildren(edgeIndex) = .TRUE.
    self%child1(edgeIndex) = self%N + 1
    self%child2(edgeIndex) = self%N + 2
    self%parent(self%N+1) = edgeIndex
    self%parent(self%N+2) = edgeIndex
    
    self%orig( self%N + 1 ) = self%orig(edgeIndex)
    self%dest( self%N + 1 ) = pInsertIndex
    self%leftFace( self%N + 1 ) = self%leftFace( edgeIndex ) 
    self%rightFace(self%N + 1 ) = self%rightFace( edgeIndex)
    
    self%orig( self%N + 2 ) = pInsertIndex
    self%dest( self%N + 2 ) = self%dest(edgeIndex)
    self%leftFace( self%N + 2 ) = self%leftFace( edgeIndex ) 
    self%rightFace(self%N + 1 ) = self%rightFace( edgeIndex)
    
    call replaceIncidentEdgeWithChild( self, edgeIndex, aParticles)
    call RecordIncidentEdgeAtParticles(self, self%N + 1, aParticles)
    call RecordIncidentEdgeAtParticles(self, self%N + 2, aParticles)
    
    self%N = self%N + 2
end subroutine



   

!> @brief Replaces an edge recorded at its incident particles with the appropriate child edge index.
!> @param self Target edge object
!> @param parentIndex index of edge that was divided
!> @param aParticles particles object associated with this set of edges
subroutine replaceIncidentEdgeWithChild( self, parentIndex, aParticles )
    type(Edges), intent(in) :: self
    integer(kint), intent(in) :: parentIndex
    type(Particles), intent(inout) :: aParticles
    !
    integer(kint) :: j, pEdgeIndex
    
    !
    ! replace parent edge at origin vertex
    !
    pEdgeIndex = 0
    do j = 1, aParticles%nEdges( self%orig(parentIndex) )
        if ( aParticles%incidentEdges( j, self%orig(parentIndex) ) == parentIndex ) then
            pEdgeIndex = j
            exit
        endif
    enddo
    
    if ( pEdgeIndex == 0 ) then
        call LogMessage(log, ERROR_LOGGING_LEVEL, logkey, "replaceIncidentEdgeWithChild ERROR : Parent not found.")
        return
    endif
    
    aParticles%incidentEdges( pEdgeIndex, self%orig(parentIndex)) = self%child1(parentIndex)
    
    !
    ! replace parent edge at destination vertex
    !
    pEdgeIndex = 0
    do j = 1, aParticles%nEdges( self%dest(parentIndex) )
        if ( aParticles%incidentEdges(j, self%dest(parentIndex)) == parentIndex ) then
            pEdgeIndex = j
            exit
        endif
    enddo
    
    if ( pEdgeIndex == 0 ) then
        call LogMessage(log, ERROR_LOGGING_LEVEL, logkey, "replaceIncidentEdgeWithChild ERROR : Parent not found.")
        return
    endif
    
    aParticles%incidentEdges( pEdgeIndex, self%dest(parentIndex)) = self%child2(parentIndex)
end subroutine

!> @brief Returns .TRUE. if an edge has positive orientation relative to the selected face
!> @todo This function should output an error if edgeIndex is not a member of faces.edges(:)
!>
!> @param anEdges Target edge object
!> @param faceIndex Integer index to a face in a facesmodule::faces object
!> @param edgeIndex Integer index to an edge
!> @return .TRUE. if edge has positive orientation relative to the chosen face
function positiveEdge( anEdges, faceIndex, edgeIndex )
    logical(klog) :: positiveEdge
    type(Edges), intent(in) :: anEdges
    integer(kint), intent(in) :: faceIndex
    integer(kint), intent(in) :: edgeIndex
    
    positiveEdge = ( faceIndex == anEdges%leftFace(edgeIndex))
end function

!> @brief Returns .TRUE. if an edge is on the boundary of a mesh (does not apply to spherical meshes, which have no boundaries)
!> @param anEdges Target edge object
!> @param edgeIndex target edge
function onBoundary( anEdges, edgeIndex )
    logical(klog) :: onBoundary
    type(Edges), intent(in) :: anEdges
    integer(kint), intent(in) :: edgeIndex
    onBoundary = ( anEdges%leftFace(edgeIndex) < 1 .OR. anEdges%rightFace(edgeIndex) < 1 )
end function

!> @brief Builds a STDIntVector containing the indices of all of an edge's children.
!>
!> @param anEdges Target edge object
!> @param parentIndex edge whose children are needed
!> @param leafEdges On output, the indices of parentEdge's children
subroutine GetLeafEdgesFromParent( anEdges, parentIndex, leafEdges )
    type(Edges), intent(in) :: anEdges
    integer(kint), intent(in) :: parentIndex
    type(STDIntVector), intent(out) :: leafEdges
    !
    integer(kint) :: i, nLeaves
    logical(klog) :: keepGoing
    
    call initialize(leafEdges)
    call leafEdges%pushBack(parentIndex)
    nLeaves = 1
    
    keepGoing = .FALSE.
    if ( anEdges%hasChildren(parentIndex) ) keepGoing = .TRUE.
    
    do while (keepGoing)
        do i = 1, nLeaves
            if ( anEdges%hasChildren( leafEdges%int(i) ) ) then
                call leafEdges%replace(i, anEdges%child1(i) )
                call leafEdges%insert(i+1, anEdges%child2(i))
            endif
        enddo
        nLeaves = leafEdges%N
        keepGoing = .FALSE.
        do i = 1, nLeaves
            if ( anEdges%hasChildren(leafEdges%int(i)) ) keepGoing = .TRUE.
        enddo
    enddo
end subroutine

!> @brief Computes the area represented by the triangles connecting the vertices of an edge's child particles and the center particle 
!> of their common face.  Used for adaptively refined meshes.
!>
!> @param self Target edge object
!> @param aParticles particles object associated with this set of edges
!> @param centerParticle index to the particle at the center of the edge's face
!> @param leafEdges indices of the child edges associated with a face whose neighbors may be more refined than itself
!> @param nLeaves number of leafEdges
!> @return Area represented by one side of a polyhedral face
function AreaFromLeafEdges( self, aParticles, centerParticle, leafEdges, nLeaves )
    real(kreal) :: AreaFromLeafEdges
    type(Edges), intent(in) :: self
    type(Particles), intent(in) :: aParticles
    integer(kint), intent(in) :: centerParticle
    integer(kint), intent(in) :: leafEdges(:)
    integer(kint), intent(in) :: nLeaves
    !
    integer(kint) :: i
    real(kreal) :: centerVec(3), v1Vec(3), v2Vec(3)
    
    AreaFromLeafEdges = 0.0_kreal
    centerVec = PhysCoord( aParticles, centerParticle )
    if ( aParticles%geomKind == PLANAR_GEOM ) then
        do i = 1, nLeaves
            v1Vec = PhysCoord( aParticles, self%orig(leafEdges(i)) )
            v2Vec = PhysCoord( aParticles, self%dest(leafEdges(i)) )
            AreaFromLeafEdges = AreaFromLeafEdges + TriArea( v1Vec(1:2), centerVec(1:2), v2Vec(1:2) )
        enddo
    elseif ( aParticles%geomKind == SPHERE_GEOM ) then
        do i = 1, nLeaves
            v1Vec = PhysCoord( aParticles, self%orig(leafEdges(i)) )
            v2Vec = PhysCoord( aParticles, self%dest(leafEdges(i)) )
            AreaFromLeafEdges = AreaFromLeafEdges + SphereTriArea( v1Vec, centerVec, v2Vec )
        enddo
    else
        call LogMessage(log,ERROR_LOGGING_LEVEL,logkey//" AreaFromLeafEdges ERROR : ", "geomKind not implemented.")
        return
    endif
end function


!
!----------------
! Private methods
!----------------
!

!> @brief Initializes a logger for the Edges module
!> 
!> Output is controlled both by message priority and by MPI Rank
!> @param aLog Target Logger object
!> @param rank Rank of this processor
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