module EdgesOOModule

use NumberKindsModule
use LoggerModule
use ParticlesOOModule
use SphereGeomModule
use PlaneGeomModule
use STDIntVectorModule

implicit none
private

public Edges, CubicEdges, LinearEdges

type, abstract :: Edges
    integer(kint), allocatable :: orig(:) !< Integer array containing indices of particlesmodule::particles     
    integer(kint), allocatable :: interiorParticles(:,:)
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
        procedure :: init 
        procedure :: copy => copyLinear
        procedure :: insert => insertLinear
        procedure(divide), deferred :: divide
        procedure :: onBoundary
        procedure :: positiveOrientation
        procedure :: getLeafEdges
        procedure :: length
        procedure :: maxLength
        procedure :: minLength
        procedure :: avgLength
        procedure :: logStats
        procedure :: countParents
        procedure :: areaFromLeaves
end type

interface
    subroutine divide(self, index, aParticles)
        import :: Edges
        import :: kint
        import :: Particles
        implicit none
        class(Edges), intent(inout) :: self
        integer(kint), intent(in) :: index
        class(Particles), intent(inout) :: aParticles
    end subroutine
end interface

type, extends(Edges) :: LinearEdges
    contains
        final :: deleteLinear
        procedure :: divide => divideLinear
end type

type, extends(Edges) :: CubicEdges
    contains
        procedure :: copy => copyCubic
        procedure :: insert => insertCubic
        procedure :: divide => divideCubic
        final :: deleteCubic
end type

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

subroutine init(self, nMax)
    class(Edges), intent(inout) :: self
    integer(kint), intent(in) :: nMax
    
    if (.NOT. logInit) call InitLogger(log, procRank)
    
    if ( nMax <= 0 ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL, logKey, " invalid nMax.")
		return
	endif
	
	self%N_Max = nMax
	self%N = 0
	
	allocate(self%orig(nmax))
	allocate(self%dest(nMax))
	allocate(self%leftFace(nMax))
	allocate(self%rightFace(nMax))
	allocate(self%hasChildren(nMax))
	allocate(self%child1(nMax))
	allocate(self%child2(nMax))
	allocate(self%parent(nMax))
	self%orig = 0
	self%dest = 0
	self%leftFace = 0
	self%rightFace = 0
	self%hasChildren = .FALSE.
	self%child1 = 0
	self%child2 = 0
	self%parent = 0
	
	select type (self)
	class is (LinearEdges)
	    ! initialize base class (already done)
	class is (CubicEdges)
	    allocate(self%interiorParticles(2,nMax))
	    self%interiorParticles = 0
	class default
	    call LogMessage(log, ERROR_LOGGING_LEVEL, logKey, " init error: invalid type.")
	    return 
	end select
end subroutine

subroutine deleteLinear(self)
    type(LinearEdges), intent(inout) :: self
    if (allocated(self%orig)) deallocate(self%orig)
    if (allocated(self%dest)) deallocate(self%dest)
    if (allocated(self%leftFace)) deallocate(self%leftFace)
    if (allocated(self%rightFace)) deallocate(self%rightFace)
    if (allocated(self%child1)) deallocate(self%child1)
    if (allocated(self%child2)) deallocate(self%child2)
    if (allocated(self%hasChildren)) deallocate(self%hasChildren)
    if (allocated(self%parent)) deallocate(self%parent)
end subroutine

subroutine deleteCubic(self)
    type(CubicEdges), intent(inout) :: self
    if (allocated(self%interiorParticles)) deallocate(self%interiorParticles)
    if (allocated(self%orig)) deallocate(self%orig)
    if (allocated(self%dest)) deallocate(self%dest)
    if (allocated(self%leftFace)) deallocate(self%leftFace)
    if (allocated(self%rightFace)) deallocate(self%rightFace)
    if (allocated(self%child1)) deallocate(self%child1)
    if (allocated(self%child2)) deallocate(self%child2)
    if (allocated(self%hasChildren)) deallocate(self%hasChildren)
    if (allocated(self%parent)) deallocate(self%parent)
end subroutine

subroutine copyLinear(self, other)
    class(Edges), intent(inout) :: self
    class(Edges), intent(in) :: other

    if (self%N < other%N) then
        call LogMessage(log, ERROR_LOGGING_LEVEL, logKey, " copy error: not enough memory.")
    endif
    
    self%orig(1:other%N) = other%orig(1:other%N)
    self%dest(1:other%N) = other%dest(1:other%N)
    self%leftFace(1:other%N) = other%leftFace(1:other%N)
    self%rightFace(1:other%N) = other%rightFace(1:other%N)
    self%hasChildren(1:other%N) = other%hasChildren(1:other%N)
    self%parent(1:other%N) = other%parent(1:other%N)
    self%child1(1:other%N) = other%child1(1:other%N)
    self%child2(1:other%N) = other%child2(1:other%N)
end subroutine

subroutine copyCubic(self, other)
    class(CubicEdges), intent(inout) :: self
    class(Edges), intent(in) :: other
    
    call copyLinear(self%edges, other)
    select type(other)
        class is (CubicEdges)
            self%interiorParticles(:,1:other%N) = other%interiorParticles(:,1:other%N)
    end select
end subroutine

!> @brief Inserts a new edge; particle allocation/defintion takes place elsewhere
subroutine insertLinear(self, origIndex, destIndex, leftFace, rightFace, intrinds)
    class(Edges), intent(inout) :: self
    integer(kint), intent(in) :: origIndex, destIndex, leftFace, rightFace
    integer(kint), dimension(:), intent(in), optional :: intrinds ! ignored for linear edges
    !
    integer(kint) :: nn
    
    nn = self%N
    
    if ( self%N >= self%N_Max ) then
        call LogMessage(log, ERROR_LOGGING_LEVEL, logKey, " InsertEdge : out of memory. ")
        return
    endif
    
    self%orig(nn+1) = origIndex
    self%dest(nn+1) = destIndex
    self%leftFace(nn+1) = leftFace
    self%rightFace(nn+1) = rightFace
    
    self%N = nn+1
end subroutine

subroutine insertCubic(self, origIndex, destIndex, leftFace, rightFace, intrInds)
    class(CubicEdges), intent(inout) :: self
    integer(kint), intent(in) :: origIndex, destIndex, leftFace, rightFace
    integer(kint), dimension(:), intent(in), optional :: intrInds
    !
    integer(kint) :: nn
    
    nn = self%N
    
    if ( self%N >= self%N_Max ) then
        call LogMessage(log, ERROR_LOGGING_LEVEL, logKey, " InsertEdge : out of memory. ")
        return
    endif
    
    self%orig(nn+1) = origIndex
    self%dest(nn+1) = destIndex
    self%leftFace(nn+1) = leftFace
    self%rightFace(nn+1) = rightFace
    if (present(intrinds)) then
        self%interiorParticles(:,nn+1) = intrInds
    endif
    self%N = nn+1
end subroutine

pure function onBoundary(self, index)
    logical(klog) :: onBoundary
    class(Edges), intent(in) :: self
    integer(kint), intent(in) :: index
    onBoundary = (self%leftFace(index) < 1 .or. self%rightFace(index) < 1)
end function

pure function positiveOrientation(self, edgeIndex, faceIndex)
    logical(klog) :: positiveOrientation
    class(Edges), intent(in) :: self
    integer(kint), intent(in) :: edgeIndex, faceIndex
    positiveOrientation = (self%leftFace(edgeIndex) == faceIndex)
end function

function getLeafEdges(self, index)
    type(STDIntVector) :: getLeafEdges
    class(Edges), intent(in) :: Self
    integer(kint), intent(in) :: index
    !
    integer(kint) :: i, nLeaves, workingEdge
    logical(klog) :: keepGoing
    
    call initialize(getLeafEdges)
    call getLeafEdges%pushBack(index)
    nLeaves = 1
    
    keepGoing = self%hasChildren(index)
    
    do while (keepGoing)
        do i = 1, nLeaves
            workingEdge = getLeafEdges%int(i)
            if (self%hasChildren(workingEdge)) then
                call getLeafEdges%replace(i, self%child1(workingEdge))
                call getLeafEdges%insert(i+1, self%child2(workingEdge))
            endif
        enddo
        nLeaves = getLeafEdges%N
        keepGoing = .FALSE.
        do i = 1, nLeaves
            if (self%hasChildren(getLeafEdges%int(i))) keepGoing = .TRUE.
        enddo
    enddo
end function

subroutine divideLinear(self, index, aParticles)
    class(LinearEdges), intent(inout) :: self
    integer(kint), intent(in) :: index
    class(Particles), intent(inout) :: aParticles
    !
    real(kreal), dimension(3) :: midPt, lagMidPt, physOrig, physDest, lagOrig, lagDest
    integer(kint) :: nParticles
    
    if (self%N + 2 > self%N_Max) then
        call LogMessage(log, ERROR_LOGGING_LEVEL, logkey, " divideEdge error : not enough memory.")
        return
    endif
    
    physOrig = aParticles%physCoord(self%orig(index))
    physDest = aParticles%physCoord(self%dest(index))
    
    lagOrig = aParticles%lagCoord(self%orig(index))
    lagDest = aParticles%lagCoord(self%dest(index))
    
    if (aParticles%geomKind == SPHERE_GEOM) then
        midPt = SphereMidpoint(physOrig, physDest)
        lagMidPt = SphereMidpoint(lagOrig, lagDest)
    else
        midPt = 0.5_kreal * (physOrig + physDest)
        lagMidPt = 0.5_kreal * (lagOrig + lagDest)
    endif
    
    nParticles = aParticles%N
    call aParticles%insert(midPt, lagMidPt)
    
    self%hasChildren(index) = .TRUE.
    self%child1(index) = self%N+1
    self%child2(index) = self%N+2
    self%parent(self%N+1) = index
    self%parent(self%N+2) = index
    
    self%orig(self%N+1) = self%orig(index)
    self%dest(self%N+1) = nParticles+1
    self%orig(self%N+2) = nParticles+1
    self%dest(self%N+2) = self%dest(index)
    
    self%leftFace(self%N+1) = self%leftFace(index)
    self%leftFace(self%N+2) = self%leftFace(index)
    self%rightFace(self%N+1) = self%rightFace(index)
    self%rightFace(self%N+2) = self%rightFace(index)
    
    self%N = self%N + 2
end subroutine

subroutine divideCubic(self, index, aParticles)
    class(CubicEdges), intent(inout) :: self
    integer(kint), intent(in) :: index
    class(Particles), intent(inout) :: aParticles
    !
    real(kreal), dimension(3) :: midPt, lagMidPt, physOrig, physDest, lagOrig, lagDest
    integer(kint) :: nParticles
    real(kreal), dimension(3,4) :: newPts, newLagPts
    
    if (self%N + 2 > self%N_Max) then
        call LogMessage(log, ERROR_LOGGING_LEVEL, logkey, " divideEdge error : not enough memory.")
        return
    endif
    
    physOrig = aParticles%physCoord(self%orig(index))
    physDest = aParticles%physCoord(self%dest(index))
    
    lagOrig = aParticles%lagCoord(self%orig(index))
    lagDest = aParticles%lagCoord(self%dest(index))
    
    if (aParticles%geomKind == SPHERE_GEOM) then
        midPt = SphereMidpoint(physOrig, physDest)
        lagMidPt = SphereMidpoint(lagOrig, lagDest)
        
        newPts(:,1) = pointAlongSphereVector(physOrig, midPt, -oosqrt5)
        newPts(:,2) = pointAlongSphereVector(physOrig, midPt, oosqrt5)
        newPts(:,3) = pointAlongSphereVector(midPt, physDest, -oosqrt5)
        newPts(:,4) = pointAlongSphereVector(midPt, physDest, oosqrt5)
        
        newLagPts(:,1) = pointAlongSphereVector(lagOrig, lagMidPt, -oosqrt5)
        newLagPts(:,2) = pointAlongSphereVector(lagOrig, lagMidPt, oosqrt5)
        newLagPts(:,3) = pointAlongSphereVector(lagMidPt, lagDest, -oosqrt5)
        newLagPts(:,4) = pointAlongSphereVector(lagMidPt, lagDest, oosqrt5)
    else
        midPt = 0.5_kreal * (physOrig + physDest)
        lagMidPt = 0.5_kreal * (lagOrig + lagDest)
        
        newPts(:,1) = pointAlongChordVector(physOrig, midPt, -oosqrt5)
        newPts(:,2) = pointAlongChordVector(physOrig, midPt, oosqrt5)
        newPts(:,3) = pointAlongChordVector(midPt, physDest, -oosqrt5)
        newPts(:,4) = pointAlongChordVector(midPt, physDest, oosqrt5)
        
        newLagPts(:,1) = pointAlongChordVector(lagOrig, lagMidPt, -oosqrt5)
        newLagPts(:,2) = pointAlongChordVector(lagOrig, lagMidPt, oosqrt5)
        newLagPts(:,3) = pointAlongChordVector(lagMidPt, lagDest, -oosqrt5)
        newLagPts(:,4) = pointAlongChordVector(lagMidPt, lagDest, oosqrt5)
    endif
    
    nParticles = aParticles%N
    
    call aParticles%replace(self%interiorParticles(1,index), newPts(:,1), newLagPts(:,1))
    call aParticles%insert(newPts(:,2), newLagPts(:,2)) ! nParticles + 1
    call aParticles%insert(midPt, lagMidPt) ! nParticles + 2 = midpoint
    call aParticles%insert(newPts(:,3), newLagPts(:,3)) ! nParticles + 3
    call aParticles%replace(self%interiorParticles(2, index), newPts(:,4), newLagPts(:,4))
    
    self%hasChildren(index) = .TRUE.
    self%child1(index) = self%N+1
    self%child2(index) = self%N+2
    self%parent(self%N+1) = index
    self%parent(self%N+2) = index
    
    self%orig(self%N+1) = self%orig(index)
    self%interiorParticles(:,self%N+1) = (/ self%interiorParticles(1,index), nParticles+1 /)
    self%dest(self%N+1) = nParticles+2
    
    self%orig(self%N+2) = nParticles+2
    self%interiorParticles(:,self%N+2) = (/ nParticles+3, self%interiorParticles(2, index) /)
    self%dest(self%N+2) = self%dest(index)
    
    self%leftFace(self%N+1) = self%leftFace(index)
    self%leftFace(self%N+2) = self%leftFace(index)
    self%rightFace(self%N+1) = self%rightFace(index)
    self%rightFace(self%N+2) = self%rightFace(index)
    
    self%N = self%N + 2
end subroutine

pure function length(self, index, aParticles)
    real(kreal) :: length
    class(Edges), intent(in) :: self
    integer(kint), intent(in) :: index
    class(Particles), intent(in) :: aParticles
    !
    real(kreal) :: orig(3), dest(3)

    orig = aParticles%physCoord(self%orig(index))
    dest = aParticles%physCoord(self%dest(index))
    if (aParticles%geomKind == SPHERE_GEOM) then
        length = SphereDistance(orig, dest)
    else
        length = ChordDistance(orig, dest)
    endif
end function

pure function maxLength(self, aParticles)
    real(kreal) :: maxLength
    class(Edges), intent(in) :: self
    class(Particles), intent(in) :: aParticles
    !
    integer(kint) :: i
    real(kreal) :: testlen
    
    maxLength = dzero
    do i=1, self%N
        if (.not. self%hasChildren(i)) then
            testlen = self%length(i, aParticles)
            if (testlen > maxLength) maxLength = testlen
        endif
    enddo
end function

pure function minLength(self, aParticles)
    real(kreal) :: minLength
    class(Edges), intent(in) :: self
    class(Particles), intent(in) :: aParticles
    !
    integer(kint) :: i
    real(kreal) :: testlen
    
    minLength = self%length(1, aParticles)
    do i=1, self%N
        if (.not. self%hasChildren(i)) then
            testlen = self%length(i, aParticles)
            if (testlen < minLength) minLength = testlen
        endif
    enddo
end function

pure function avgLength(self, aParticles)
    real(kreal) :: avgLength
    class(Edges), intent(in) :: Self
    class(Particles), intent(in) :: aparticles
    !
    integer(kint) :: i, nleaves
    
    nleaves = 0
    avgLength = dzero
    do i=1, self%N
        if (.not. self%hasChildren(i)) then
            avgLength = avgLength + self%length(i, aParticles)
            nLeaves = nLeaves + 1
        endif
    enddo
    avgLength = avgLength / real(nLeaves, kreal)
end function

subroutine logStats(self, alog)
    class(Edges), intent(in) :: self
    type(Logger), intent(inout) :: alog
    call LogMessage(aLog, TRACE_LOGGING_LEVEL, logkey, " Edges Stats : ")
    call StartSection(aLog)
    call LogMessage(aLog, TRACE_LOGGING_LEVEL, "edges.N = ", self%N )
    call LogMessage(aLog, TRACE_LOGGING_LEVEL, "edges.N_Max = ", self%N_Max)
    call LogMessage(aLog, TRACE_LOGGING_LEVEL, "n divided edges = ", count(self%hasChildren) )
    call LogMessage(aLog, TRACE_LOGGING_LEVEL, "n leaf edges = ", self%N - count(self%hasChildren))
    call EndSection(aLog)
end subroutine

pure function countParents(self, index)
    integer(kint) :: countParents
    class(Edges), intent(in) :: self
    integer(kint), intent(in) :: index
    !
    logical(klog) :: keepGoing
    integer(kint) :: parentIndex
    
    countParents = 0
    keepGoing = (self%parent(index) > 0)
    parentIndex = self%parent(index)
    do while (keepGoing)
        countParents = countParents + 1
        parentIndex = self%parent(parentIndex)
        keepGoing = (self%parent(index) > 0)
    enddo
end function

pure function areaFromLeaves(self, aParticles, centroid, leaves)
    real(kreal) :: areaFromLeaves
    class(Edges), intent(in) :: self
    class(Particles), intent(in) :: aParticles
    real(kreal), dimension(3), intent(in) :: centroid
    integer(kint), dimension(:), intent(in) :: leaves
    !
    integer(kint) :: i, nLeaves
    real(kreal), dimension(3) :: orig, dest
    
    areaFromLeaves = dzero
    nLeaves = size(leaves)
    if (aParticles%geomKind == PLANAR_GEOM) then
        do i=1, nLeaves
            orig = aParticles%physCoord(self%orig(leaves(i)))
            dest = aParticles%physCoord(self%dest(leaves(i)))
            areaFromLeaves = areaFromLeaves + TriArea(orig, centroid, dest)
        enddo
    elseif (aParticles%geomKind == SPHERE_GEOM) then
        do i=1, nLeaves
            orig = aParticles%physCoord(self%orig(leaves(i)))
            dest = aParticles%physCoord(self%dest(leaves(i)))
            areaFromLeaves = areaFromLeaves + SphereTriArea(orig, centroid, dest)
        enddo
    endif
end function

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

end module