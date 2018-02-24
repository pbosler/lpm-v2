module CubicQuadFacesModule

use NumberKindsModule
use LoggerModule
use ParticlesModule
use SphereGeomModule
use PlaneGeomModule
use EdgesModule
use CubicEdgesModule

implicit none
private

public CubicQuadFaces
public New, Delete, Copy
public InsertFace, DivideCubicQuadFace
public PhysicalCentroid, LagrangianCentroid !, FaceArea
public LogStats

type CubicQuadFaces
    integer(kint), allocatable :: interiorParticles(:,:)
    integer(kint), allocatable :: edges(:,:)
    logical(klog), allocatable :: hasChildren(:)
    integer(kint), allocatable :: children(:,:)
    integer(kint), allocatable :: parent(:)
    integer(kint) :: N = 0
    integer(kint) :: N_Active = 0
    integer(kint) :: N_Max = 0
    
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

!
!----------------
! Logging
!----------------
!
logical(klog), save :: logInit = .FALSE.
type(Logger) :: log
character(len=28), save :: logKey = 'CubicQuadFaces'
integer(kint), parameter :: logLevel = DEBUG_LOGGING_LEVEL
character(len=MAX_STRING_LENGTH) :: logstring

contains

subroutine newPrivate(self, nMax)
    type(CubicQuadFaces), intent(out) :: self
    integer(kint), intent(in) :: nMax
    
    if (.NOT. logInit) call InitLogger(log, procRank)
    
    if ( nMax <= 0 ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL, logKey, " NewCubicQuadFaces ERROR : invalid nMax.")
		return
	endif
	
	self%N_Max = nMax
	
	allocate(self%interiorParticles(4,nMax))
	allocate(self%edges(4,nMax))
	allocate(self%hasChildren(nMax))
	allocate(self%children(4,nMax))
	allocate(self%parent(nMax))
	
	self%interiorParticles = 0
	self%edges = 0
	self%hasChildren = .FALSE.
	self%children = 0
	self%parent = 0
end subroutine

subroutine deletePrivate(self)
    type(CubicQuadFaces), intent(inout) :: self
    if (allocated(self%interiorParticles)) then
        deallocate(self%interiorParticles)
        deallocate(self%edges)
        deallocate(self%hasChildren)
        deallocate(self%children)
        deallocate(self%parent)
    endif
end subroutine

subroutine copyPrivate(self, other)
    type(CubicQuadFaces), intent(inout) :: self
    type(CubicQuadFaces), intent(in) :: other
    !
    integer(kint) :: j
    
    if ( self%N_Max < other%N ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL, logkey, " CopyCubicQuadFaces ERROR : not enough memory.")
		return
	endif
    
    self%N = other%N
    self%N_Active = other%N_Active
    do j = 1, other%N
        self%edges(:,j) = other%edges(:,j)
        self%interiorParticles(:,j) = other%interiorParticles(:,j)
        self%hasChildren(j) = other%hasChildren(j)
        self%children(:,j) = other%children(:,j)
        self%parent(j) = other%parent(j)
    enddo
end subroutine

subroutine LogStatsPrivate(self, aLog)
    type(CubicQuadFaces), intent(in) :: self
    type(Logger), intent(inout) :: aLog
    call LogMessage(aLog, TRACE_LOGGING_LEVEL, logkey, " CubicQuadFaces stats : ")
	call StartSection(aLog)
	call LogMessage(alog, TRACE_LOGGING_LEVEL, "cubicQuadFaces.N = ", self%N)
	call LogMessage(aLog, TRACE_LOGGING_LEVEL, "cubicQuadFaces.N_Max = ", self%N_Max )
	call LogMessage(aLog, TRACE_LOGGING_LEVEL, "cubicQuadFaces.N_Active = ", self%N_Active)
	call LogMessage(aLog, TRACE_LOGGING_LEVEL, "n divided faces = ", count(self%hasChildren))
	call EndSection(aLog)
end subroutine

subroutine InsertFace(self, interiorIndices, edgeIndices)
    type(CubicQuadFaces), intent(inout) :: self
    integer(kint), dimension(4) :: interiorIndices, edgeIndices
    
    if ( self%N+1 > self%N_Max ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL, logkey," InsertFace ERROR : not enough memory.")
		return
	endif
	
	self%interiorParticles(:,self%N+1) = interiorIndices
	self%edges(:,self%N+1) = edgeIndices
	self%N = self%N + 1
end subroutine

subroutine DivideCubicQuadFace(self, index, aParticles, anEdges)
    type(CubicQuadFaces), intent(inout) :: self
    integer(kint), intent(in) :: index
    type(Particles), intent(inout) :: aParticles
    type(CubicEdges), intent(inout) :: anEdges
    !
    integer(kint) :: i, j, parentEdge, edgeKid1, edgeKid2, nParticles, ctr
    integer(kint), dimension(4) :: parentVerts
    integer(kint), dimension(4,4) :: newFaceVerts, newFaceEdges, newIntInds
    real(kreal), parameter :: denomSqrt5 = 1.0_kreal / sqrt(5.0_kreal)
    real(kreal), dimension(3) :: center, lagCenter, v0, v1, v2, v3, l0, l1, l2, l3
    real(kreal), dimension(3,4) :: interiorXyz, lagInterior
    real(kreal), dimension(3,8) :: edgeMidpoints, lagEdgeMidpoints
    
    if ( self%N_Max < self%N + 4 ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL, logkey, " DivideQuadFace ERROR : not enough memory.")
		return
	endif
	
	newFaceVerts = 0
	newFaceEdges = 0
	
	parentVerts = getVertices(self, index, anEdges)
	do i=1,4
	    newFaceVerts(i,i) = parentVerts(i)
	enddo
	
	
	!
	!   loop over parent edges
	!
	do i=1,4
	    parentEdge = self%edges(i, index)
	    if (anEdges%hasChildren(parentEdge)) then
	        edgeKid1 = anEdges%child1(parentEdge)
	        edgeKid2 = anEdges%child2(parentEdge)
	    else
	        edgeKid1 = anEdges%N+1
	        edgeKid2 = anEdges%N+2
	        call DivideEdge(anEdges, parentEdge, aParticles)
	    endif
	    
	    if (positiveEdge(anEdges%edges, index, parentEdge)) then
	        newFaceEdges(i,i) = edgeKid1
	        anEdges%leftFace(edgeKid1) = self%N + i
	        
	        newFaceEdges(i, mod(i,4)+1) = edgeKid2
	        anEdges%leftFace(edgeKid2) = self%N + mod(i,4) + 1
	    else
	        newFaceEdges(i,i) = edgeKid2
	        anEdges%rightFace(edgeKid2) = self%N + i
	        
	        newFaceEdges(i, mod(i,4)+1) = edgeKid1
	        anEdges%rightFace(edgeKid1) = self%N + mod(i,4) + 1
	    endif
	    
	    newFaceVerts(mod(i,4)+1, i) = anEdges%dest(edgeKid1)
	    newFaceVerts(i, mod(i,4)+1) = anEdges%dest(edgeKid1)
	enddo
	
	!
	!   insert parent panel centroid
	!
	nParticles = aParticles%N
	center = PhysicalCentroid(self, index, anEdges, aParticles)
	lagCenter = LagrangianCentroid(self, index, anEdges, aParticles)
	call InsertParticle(aParticles, center, lagCenter) ! nParticles + 1
	do i=1,4
	    newFaceVerts(mod(i+1,4)+1,i) = nParticles+1
	enddo
	
	do j=1,4
	    do i=1,4
	        if (newFaceVerts(i,j) == 0) then
	            write(logstring,*) " vertex connectivity ERROR at parent face ", index, ", child ", i, ", vertex ", j
				call LogMessage(log, ERROR_LOGGING_LEVEL,logkey//" DivideCubicQuadFace :",logstring)
	        endif
	    enddo
	enddo
	
	!
	!   re-assign parent interior particle indices to child panel 1 
	!
	newIntInds(:,1) = self%interiorParticles(:,index)
	v0 = PhysCoord(aParticles, newFaceVerts(1,1))
	v1 = PhysCoord(aParticles, newFaceVerts(2,1))
	v2 = PhysCoord(aParticles, newFaceVerts(3,1))
	v3 = PhysCoord(aParticles, newFaceVerts(4,1))
	l0 = LagCoord(aParticles, newFaceVerts(1,1))
	l1 = LagCoord(aParticles, newFaceVerts(2,1))
	l2 = LagCoord(aParticles, newFaceVerts(3,1))
	l3 = LagCoord(aParticles, newFaceVerts(4,1))
	interiorXyz = calcInteriorPoints(v0, v1, v2, v3)
	lagInterior = calcInteriorPoints(l0, l1, l2, l3)
	edgeMidpoints(:,1) = bilinearMap(v0, v1, v2, v3, -denomSqrt5, -1.0_kreal)
	edgeMidpoints(:,2) = bilinearMap(v0, v1, v2, v3, denomSqrt5, -1.0_kreal)
	edgeMidpoints(:,7) = bilinearMap(v0, v1, v2, v3, 1.0_kreal, -denomSqrt5)
	edgeMidpoints(:,8) = bilinearMap(v0, v1, v2, v3, 1.0_kreal, denomSqrt5)
	lagEdgeMidpoints(:,1) = bilinearMap(l0, l1, l2, l3, -denomSqrt5, -1.0_kreal)
	lagEdgeMidpoints(:,2) = bilinearMap(l0, l1, l2, l3, denomSqrt5, -1.0_kreal)
	lagEdgeMidpoints(:,7) = bilinearMap(l0, l1, l2, l3, 1.0_kreal, -denomSqrt5)
	lagEdgeMidpoints(:,8) = bilinearMap(l0, l1, l2, l3, 1.0_kreal, denomSqrt5)
	
	do i=1,4
	    call ReplaceParticle(aParticles, newIntInds(i,1), interiorXyz(:,i), lagInterior(:,i))
	enddo
	
	!
	!   create interior particles for remaining children
	!
	ctr = 1
	do j=2,4
	    v0 = PhysCoord(aParticles, newFaceVerts(1,j))
	    v1 = PhysCoord(aParticles, newFaceVerts(2,j))
	    v2 = PhysCoord(aParticles, newFaceVerts(3,j))
	    v3 = PhysCoord(aParticles, newFaceVerts(4,j))
	    l0 = LagCoord(aParticles, newFaceVerts(1,j))
	    l1 = LagCoord(aParticles, newFaceVerts(2,j))
	    l2 = LagCoord(aParticles, newFaceVerts(3,j))
	    l3 = LagCoord(aParticles, newFaceVerts(4,j))
	    interiorXyz = calcInteriorPoints(v0, v1, v2, v3)
        lagInterior = calcInteriorPoints(l0, l1, l2, l3)
        if (j==3) then
            edgeMidpoints(:,3) = bilinearMap(v0, v1, v2, v3, -denomSqrt5, 1.0_kreal)
            edgeMidpoints(:,4) = bilinearMap(v0, v1, v2, v3, denomSqrt5, 1.0_kreal)
            edgeMidpoints(:,5) = bilinearMap(v0, v1, v2, v3, -1.0_kreal, -denomSqrt5)
            edgeMidpoints(:,6) = bilinearMap(v0, v1, v2, v3, -1.0_kreal, denomSqrt5)
            lagEdgeMidpoints(:,3) = bilinearMap(l0, l1, l2, l3, -denomSqrt5, 1.0_kreal)
            lagEdgeMidpoints(:,4) = bilinearMap(l0, l1, l2, l3, denomSqrt5, 1.0_kreal)
            lagEdgeMidpoints(:,5) = bilinearMap(l0, l1, l2, l3, -1.0_kreal, -denomSqrt5)
            lagEdgeMidpoints(:,6) = bilinearMap(l0, l1, l2, l3, -1.0_kreal, denomSqrt5)
        endif
        do i=1,4
            call InsertParticle(aParticles, interiorXyz(:,i), lagInterior(:,i))
            newIntInds(i,j) = nParticles + 1 + ctr
            ctr = ctr + 1
        enddo
	enddo
	
	!
	!   create new interior edges
	! 
	do i=1,8
	    call InsertParticle(aParticles, edgeMidpoints(:,i), lagEdgeMidpoints(:,i)) ! nParticles + 13 + i
	enddo
	
	call InsertEdge(anEdges, newFaceVerts(2,1), newFaceVerts(3,1), self%N+1, self%N+2, (/nParticles+14, nParticles+15/))
	newFaceEdges(2,1) = anEdges%N
	newFaceEdges(4,2) = anEdges%N
	
	call InsertEdge(anEdges, newFaceVerts(1,3), newFaceVerts(4,3), self%N+4, self%N+3, (/nParticles+16, nParticles+17/))
	newFaceEdges(4,3) = anEdges%N
	newFaceEdges(2,4) = anEdges%N
	
	call InsertEdge(anEdges, newFaceVerts(3,2), newFaceVerts(4,2), self%N+2, self%N+3, (/nParticles+18, nParticles+19/))
	newFaceEdges(3,2) = anEdges%N
	newFaceEdges(1,3) = anEdges%N
	
	call InsertEdge(anEdges, newFaceVerts(2,4), newFaceVerts(1,4), self%N+1, self%N+4, (/nParticles+20, nParticles+21/))
	newFaceEdges(1,4) = anEdges%N
	newFaceEdges(3,1) = anEdges%N
	
	do i = 1, 4
		do j = 1, 4
			if ( newFaceEdges(j,i) < 1 ) then
				write(logstring,*) " edge connectivity ERROR at parent face ", index, ", child ", i, ", vertex ", j
				call LogMessage(log, ERROR_LOGGING_LEVEL, logkey//" DivideCubicQuadFace :",logstring)	
			endif
		enddo
	enddo
	
	!
	!   create child faces
	!
	do i=1,4
	    self%interiorParticles(:,i) = newIntInds(:,i)
	    self%edges(:,self%N+i) = newFaceEdges(:,i)
	    self%children(i,index) = self%N + i
	    self%parent(self%N+i) = index
	enddo
	self%hasChildren(index) = .TRUE.
	self%N = self%N + 4
	self%N_Active = self%N_Active + 3
end subroutine

function PhysicalCentroid(self, index, anEdges, aParticles)
    real(kreal), dimension(3) :: PhysicalCentroid
    type(CubicQuadFaces), intent(in) :: self
    integer(kint), intent(in) :: index
    type(CubicEdges), intent(in) :: anEdges
    type(Particles), intent(in) :: aParticles
    !
    integer(kint) :: i
    integer(kint), dimension(4) :: vertInds
    
    PhysicalCentroid = 0.0_kreal
    vertInds = getVertices(self, index, anEdges)
    do i=1,4
        PhysicalCentroid = PhysicalCentroid + PhysCoord(aParticles, vertInds(i))
    enddo 
    PhysicalCentroid = 0.25_kreal * PhysicalCentroid
end function

function LagrangianCentroid(self, index, anEdges, aParticles)
    real(kreal), dimension(3) :: LagrangianCentroid
    type(CubicQuadFaces), intent(in) :: self
    integer(kint), intent(in) :: index
    type(CubicEdges), intent(in) :: anEdges
    type(Particles), intent(in) :: aParticles
    !
    integer(kint) :: i
    integer(kint), dimension(4) :: vertInds
    
    LagrangianCentroid = 0.0_kreal
    vertInds = getVertices(self, index, anEdges)
    do i=1,4
        LagrangianCentroid = LagrangianCentroid + LagCoord(aParticles, vertInds(i))
    enddo 
    LagrangianCentroid = 0.25_kreal * LagrangianCentroid
end function


function getVertices(self, index, anEdges)
    integer(kint), dimension(4) :: getVertices
    type(CubicQuadFaces), intent(in) :: self
    integer(kint), intent(in) :: index
    type(CubicEdges), intent(in) :: anEdges
    !
    integer(kint) :: j
    
    do j=1,4
        if (positiveEdge(anEdges%edges, index, self%edges(j, index))) then
            getVertices(j) = anEdges%orig(self%edges(j, index))
        else
            getVertices(j) = anEdges%dest(self%edges(j, index))
        endif
    enddo
end function

pure function calcInteriorPoints( v0, v1, v2, v3)
    real(kreal), dimension(3,4) :: calcInteriorPoints
    real(kreal), dimension(3), intent(in) :: v0, v1, v2, v3
    !
    real(kreal), parameter :: denomSqrt5 = 1.0_kreal / sqrt(5.0_kreal)
    
    calcInteriorPoints(:,1) = bilinearMap(v0, v1, v2, v3, -denomSqrt5, denomSqrt5)
    calcInteriorPoints(:,2) = bilinearMap(v0, v1, v2, v3, -denomSqrt5, -denomSqrt5)
    calcInteriorPoints(:,3) = bilinearMap(v0, v1, v2, v3, denomSqrt5, -denomSqrt5)
    calcInteriorPoints(:,4) = bilinearMap(v0, v1, v2, v3, denomSqrt5, denomSqrt5)
end function

pure function bilinearMap(v0, v1, v2, v3, s1, s2)
    real(kreal), dimension(3) :: bilinearMap
    real(kreal), dimension(3), intent(in) :: v0, v1, v2, v3
    real(kreal), intent(in) :: s1, s2
    
    bilinearMap = 0.25_kreal * ( (1.0_kreal-s1)*(1.0_kreal-s2)*v0 + (1.0_kreal+s1)*(1.0_kreal-s2)*v1 + &
        (1.0_kreal+s1)*(1.0_kreal+s2)*v3 + (1.0_kreal-s1)*(1.0_kreal+s2)*v3)
end function

end module