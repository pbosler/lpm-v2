module EdgesOOModule
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
use UtilitiesModule
use LoggerModule
use ParticlesOOModule
use SphereGeomModule
use PlaneGeomModule
use STDIntVectorModule
use CubicGLLModule

implicit none
private

public Edges, CubicEdges, LinearEdges

!> @brief Edges know the indices (to @ref Particles) of their origin and destination, and the indices (to @ref Faces)
!> of their left face and right face.
type, abstract :: Edges
    integer(kint), allocatable :: orig(:) !< Integer array containing indices of particlesmodule::particles
    integer(kint), allocatable :: interiorParticles(:,:) !< Used for high-order edges only
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
        procedure :: writeMatlab
        procedure :: vector
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

!> @brief Allocates memory for a new edges object.  All values are zeroed and must be initialized separately.
!> @param self Target edge object
!> @param nMax max number of edges allowed in memory
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

!> @brief Deletes and frees memory associated with an edges object
!> @param self Target edge object
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

!> @brief Deletes and frees memory associated with an edges object
!> @param self Target edge object
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

!> @brief Writes the edges data structure to a .m file for later reading by Matlab.  For debugging.
!> @param self Target edge object
!> @param fileunit integer unit of .m file
subroutine writeMatlab(self, fileunit)
    class(Edges), intent(in) :: self
    integer(kint), intent(in) :: fileunit
    !
    integer(kint) :: i, j, nint

    write(fileunit,'(A)',advance='no') 'edge_origs = ['
    do i=1, self%N-1
        write(fileunit,'(I8,A)',advance='no') self%orig(i), ','
    enddo
    write(fileunit, *) self%orig(self%N), '];'

    write(fileunit,'(A)',advance='no') 'edge_dests = ['
    do i=1, self%N-1
        write(fileunit,'(I8,A)',advance='no') self%dest(i), ','
    enddo
    write(fileunit, *) self%dest(self%N), '];'

    write(fileunit,'(A)',advance='no') 'edge_lefts = ['
    do i=1, self%N-1
        write(fileunit,'(I8,A)',advance='no') self%leftFace(i), ','
    enddo
    write(fileunit, *) self%leftFace(self%N), '];'

    write(fileunit,'(A)',advance='no') 'edge_rights = ['
    do i=1, self%N-1
        write(fileunit,'(I8,A)',advance='no') self%rightFace(i), ','
    enddo
    write(fileunit, *) self%rightFace(self%N), '];'

    write(fileunit,'(A)',advance='no') 'edge_kids = ['
    do i=1, self%N-1
        write(fileunit,'(2(I8,A))',advance='no') self%child1(i), ',', self%child2(i), ';'
    enddo
    write(fileunit, *) self%child1(self%N), ', ', self%child2(self%N), '];'

    if (allocated(self%interiorParticles)) then
        nint = size(self%interiorParticles,1)
        write(fileunit,'(A)',advance='no') 'edge_interiors = ['
        do i=1,self%N-1
            do j=1, nint-1
                write(fileunit,'(I8,A)',advance='no') self%interiorParticles(j,i), ','
            enddo
            write(fileunit,*) self%interiorParticles(nint,i), ';'
        enddo
        do j=1, nint-1
            write(fileunit,'(I8,A)',advance='no') self%interiorParticles(j,self%N), ','
        enddo
        write(fileunit,*) self%interiorParticles(nint,self%N), '];'
    endif
end subroutine

!> @brief Performs a deep copy of one edges object into another.  They must have been allocated the same to avoid memory errors.
!> @param self Target edge object
!> @param other Source edge object
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

!> @brief Performs a deep copy of one edges object into another.  They must have been allocated the same to avoid memory errors.
!> @param self Target edge object
!> @param other Source edge object
subroutine copyCubic(self, other)
    class(CubicEdges), intent(inout) :: self
    class(Edges), intent(in) :: other
    select type(other)
        class is (CubicEdges)
            self%orig(1:other%N) = other%orig(1:other%N)
            self%dest(1:other%N) = other%dest(1:other%N)
            self%leftFace(1:other%N) = other%leftFace(1:other%N)
            self%rightFace(1:other%N) = other%rightFace(1:other%N)
            self%hasChildren(1:other%N) = other%hasChildren(1:other%N)
            self%parent(1:other%N) = other%parent(1:other%N)
            self%child1(1:other%N) = other%child1(1:other%N)
            self%child2(1:other%N) = other%child2(1:other%N)
            self%interiorParticles(:,1:other%N) = other%interiorParticles(:,1:other%N)
            self%N = other%N
        class default
            call LogMessage(log, ERROR_LOGGING_LEVEL, trim(logkey)//" copy cubic edges error: ", "type mismatch.")
    end select
end subroutine

!> @brief Inserts an edge to the end of an edges object.  Increase edges.N by one.
!> Particle and Face insertion/definition must take place elsewhere.
!> @param self Target edge object
!> @param origIndex index to origin particle
!> @param destIndex index to destination particle
!> @param leftFace index to left face of new edge in faces object
!> @param rightFace index to right face of new edge in faces object
!> @param unused (low-order)
subroutine insertLinear(self, origIndex, destIndex, leftFace, rightFace, intrinds)
    class(Edges), intent(inout) :: self
    integer(kint), intent(in) :: origIndex, destIndex, leftFace, rightFace
    integer(kint), dimension(:), intent(in), optional :: intrinds ! ignored for linear edges
    !
    integer(kint) :: nn

    nn = self%N

    if ( self%N > self%N_Max ) then
        call LogMessage(log, ERROR_LOGGING_LEVEL, logKey, " InsertEdge : out of memory. ")
        return
    endif

    self%orig(nn+1) = origIndex
    self%dest(nn+1) = destIndex
    self%leftFace(nn+1) = leftFace
    self%rightFace(nn+1) = rightFace

    self%N = nn+1
end subroutine

!> @brief Inserts an edge to the end of an edges object.  Increase edges.N by one.
!> Particle and Face insertion/definition must take place elsewhere.
!> @param self Target edge object
!> @param origIndex index to origin particle
!> @param destIndex index to destination particle
!> @param leftFace index to left face of new edge in faces object
!> @param rightFace index to right face of new edge in faces object
!> @param intrinds indices of the 2 interior particles of new cubic edge
subroutine insertCubic(self, origIndex, destIndex, leftFace, rightFace, intrInds)
    class(CubicEdges), intent(inout) :: self
    integer(kint), intent(in) :: origIndex, destIndex, leftFace, rightFace
    integer(kint), dimension(:), intent(in), optional :: intrInds
    !
    integer(kint) :: nn

    nn = self%N

    if ( self%N > self%N_Max ) then
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

!> @brief Returns .TRUE. if an edge is on the boundary of a mesh (does not apply to spherical meshes, which have no boundaries)
!> @param anEdges Target edge object
!> @param edgeIndex target edge
pure function onBoundary(self, index)
    logical(klog) :: onBoundary
    class(Edges), intent(in) :: self
    integer(kint), intent(in) :: index
    onBoundary = (self%leftFace(index) < 1 .or. self%rightFace(index) < 1)
end function

!> @brief Returns .TRUE. if an edge has positive orientation relative to the selected face
!> @todo This function should output an error if edgeIndex is not a member of faces.edges(:)
!>
!> @param self Target edge object
!> @param faceIndex Integer index to a face in a facesmodule::faces object
!> @param edgeIndex Integer index to an edge
!> @return .TRUE. if edge has positive orientation relative to the chosen face
pure function positiveOrientation(self, edgeIndex, faceIndex)
    logical(klog) :: positiveOrientation
    class(Edges), intent(in) :: self
    integer(kint), intent(in) :: edgeIndex, faceIndex
    positiveOrientation = (self%leftFace(edgeIndex) == faceIndex)
end function

!> @brief Builds a STDIntVector containing the indices of all of an edge's children.
!>
!> @param anEdges Target edge object
!> @param parentIndex edge whose children are needed
!> @param leafEdges On output, the indices of parentEdge's children
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

!> @brief Divides a parent edge, creating two child edges with the same orientation and adds particles to a Particles object.
!> @param self Target edge object
!> @param edgeIndex target edge
!> @param aParticles particles object associated with this set of edges
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

!> @brief Divides a parent edge, creating two child edges with the same orientation and adds particles to a Particles object.
!> @param self Target edge object
!> @param edgeIndex target edge
!> @param aParticles particles object associated with this set of edges
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

        newPts(:,1) = pointAlongSphereVector(physOrig, midPt, -oor5)
        newPts(:,2) = pointAlongSphereVector(physOrig, midPt, oor5)
        newPts(:,3) = pointAlongSphereVector(midPt, physDest, -oor5)
        newPts(:,4) = pointAlongSphereVector(midPt, physDest, oor5)

        newLagPts(:,1) = pointAlongSphereVector(lagOrig, lagMidPt, -oor5)
        newLagPts(:,2) = pointAlongSphereVector(lagOrig, lagMidPt, oor5)
        newLagPts(:,3) = pointAlongSphereVector(lagMidPt, lagDest, -oor5)
        newLagPts(:,4) = pointAlongSphereVector(lagMidPt, lagDest, oor5)
    else
        midPt = 0.5_kreal * (physOrig + physDest)
        lagMidPt = 0.5_kreal * (lagOrig + lagDest)

        newPts(:,1) = pointAlongChordVector(physOrig, midPt, -oor5)
        newPts(:,2) = pointAlongChordVector(physOrig, midPt, oor5)
        newPts(:,3) = pointAlongChordVector(midPt, physDest, -oor5)
        newPts(:,4) = pointAlongChordVector(midPt, physDest, oor5)

        newLagPts(:,1) = pointAlongChordVector(lagOrig, lagMidPt, -oor5)
        newLagPts(:,2) = pointAlongChordVector(lagOrig, lagMidPt, oor5)
        newLagPts(:,3) = pointAlongChordVector(lagMidPt, lagDest, -oor5)
        newLagPts(:,4) = pointAlongChordVector(lagMidPt, lagDest, oor5)
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

!>@brief Returns a vector pointing from edge orig to edge dest
!>@param self Target edges object
!>@param edgeIndex
!>@param aParticles particles object associated with this set of edges
pure function vector(self, index, aParticles)
    real(kreal), dimension(3) :: vector
    class(Edges), intent(in) :: self
    integer(kint), intent(in) :: index
    class(Particles), intent(in) :: aParticles
    !
    real(kreal), dimension(3) :: v1, v2
    
    v1 = aParticles%physCoord(self%orig(index))
    v2 = aParticles%physCoord(self%dest(index))
    vector = v2-v1
end function

!> @brief Computes the length of the edge at edgeIndex
!> @param self Target edge object
!> @param edgeIndex location of edge whose length is needed
!> @param aParticles particles object associated with this set of edges
!> @return Length of edge
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

!> @brief Computes the maximum edge length amongst all leaf edges
!> @param self Target edge object
!> @param aParticles particles object associated with this set of edges
!> @return Minimum edge length
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

!> @brief Computes the minimum edge length amongst all leaf edges
!> @param self Target edge object
!> @param aParticles particles object associated with this set of edges
!> @return Minimum edge length
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

!> @brief Computes the average edge length amongst all leaf edges
!> @param self Target edge object
!> @param aParticles particles object associated with this set of edges
!> @return Minimum edge length
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

!> @brief Log summarizing statistics about an edges object
!> @param self Target edge object
!> @param aLog Target loggermodule::logger object for output
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

!> @brief Counts the number of levels in the edges binary tree above the edge at index
!> @param self Target edge object
!> @param index index of target edge
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

!> @brief Computes the area represented by the triangles connecting the vertices of an edge's child particles and the center particle 
!> of their common face.  Used for adaptively refined meshes.
!>
!> @param self Target edge object
!> @param aParticles particles object associated with this set of edges
!> @param centroid coordinates of the face centroid
!> @param leaves indices of the child edges associated with a face whose neighbors may be more refined than itself
!> @param nLeaves number of leafEdges
!> @return areaFromLeaves represented by one side of a polyhedral face
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

!>@}
end module
