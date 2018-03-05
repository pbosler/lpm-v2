module FacesOOModule

use NumberKindsModule
use UtilitiesModule
use STDIntVectorModule
use LoggerModule
use ParticlesOOModule
use EdgesOOModule
use PlaneGeomModule
use SphereGeomModule

implicit none
private

public Faces, TriLinearFaces, QuadLinearFaces, QuadCubicFaces

type, abstract :: Faces
    integer(kint), allocatable :: centerParticles(:,:)
    integer(kint), allocatable :: vertices(:,:)  !< indices to faces' vertex particles in a particlesmodule::particles object
	integer(kint), allocatable :: edges(:,:) !< indices to faces' edges in an edgesmodule::edges object
	logical(klog), allocatable :: hasChildren(:) !< hasChildren(i) is .TRUE. if face i has been divided
	integer(kint), allocatable :: children(:,:) !< indices to child faces in a faces object
	integer(kint), allocatable :: parent(:) !< indices to parent faces in a faces object
	integer(kint) :: N !< Number of faces currently in use
	integer(kint) :: N_Active !< Number of undivided faces; these define the spatial discretization
	integer(kint) :: N_Max !< Maximum number of faces allowed in memory
	real(kreal), allocatable :: area(:)

	contains
	    procedure :: init
	    procedure :: copy
        procedure :: insert
        procedure(divide), deferred :: divide
        procedure(setArea), deferred :: setArea
        procedure :: physCentroid
        procedure :: lagCentroid
        procedure :: countParents
        procedure :: sharedEdge
        procedure :: logStats
        procedure :: writeMatlab
end type

interface
    subroutine divide(self, index, aParticles, anEdges)
        import :: Faces
        import :: Particles
        import :: Edges
        import :: kint
        implicit none
        class(Faces), intent(inout) :: self
        integer(kint), intent(in) :: index
        class(Particles), intent(inout) :: aParticles
        class(Edges), intent(inout) :: anEdges
    end subroutine
end interface



interface
    pure function setArea(self, index, aParticles)
        import :: Faces
        import :: Particles
        import :: kint
        import :: kreal
        implicit none
        real(kreal) :: setArea
        class(Faces), intent(in) :: self
        integer(kint), intent(in) :: index
        class(Particles), intent(in) :: aParticles
    end function
end interface


type, extends(Faces) :: TriLinearFaces
    contains
        procedure :: divide => divideTri
        final :: deleteTri
        procedure :: particleOppositeTriEdge
        procedure :: setarea => triFaceArea
end type

type, extends(Faces) :: QuadLinearFaces
    contains
        procedure :: divide => divideQuadLinear
        final :: deleteQuad
        procedure :: setarea => quadFaceArea
end type

type, extends(Faces) :: QuadCubicFaces
    contains
        procedure :: divide => divideQuadCubic
        procedure :: setarea => quadCubicArea
        final :: deleteQuadCubic
        procedure, private :: calcInteriorPts
        procedure, private :: getVerticesFromEdge
end type

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

subroutine init(self, faceKind, nMax)
    class(Faces), intent(inout) :: self
    integer(kint), intent(in) :: faceKind, nMax

    if (.not. logInit) call InitLogger(log, procRank)

	if ( nMax <= 0 ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL, logKey, " NewFaces ERROR : invalid nMax.")
		return
	endif

	if ( faceKind == QUAD_PANEL) then
		allocate(self%vertices(4,nMax))
		allocate(self%edges(4,nMax))
        allocate(self%centerParticles(1,nMax))
	elseif ( faceKind == TRI_PANEL ) then
		allocate(self%vertices(3,nMax))
		allocate(self%edges(3,nMax))
		allocate(self%centerParticles(1,nMax))
	elseif (faceKind == QUAD_CUBIC_PANEL) then
	    allocate(self%vertices(12, nMax))
	    allocate(self%edges(4,nMax))
	    allocate(self%centerParticles(4,nMax))
	else
		call LogMessage(log, ERROR_LOGGING_LEVEL, logKey, " NewFaces ERROR : invalid face kind.")
		return
	endif
	self%centerParticles = 0
    self%vertices = 0
    self%edges = 0

    allocate(self%children(4,nMax))
    allocate(self%hasChildren(nMax))
    allocate(self%parent(nMax))
    allocate(self%area(nMax))
    self%children = 0
    self%hasChildren = .FALSE.
    self%parent = 0
    self%area = dzero
    self%N = 0
    self%N_Active = 0
    self%N_Max = nMax

!    self%faceKind = faceKind
end subroutine

subroutine deleteTri(self)
    type(TriLinearFaces), intent(inout) :: self
    if (allocated(self%vertices)) deallocate(self%vertices)
    if (allocated(self%edges)) deallocate(self%edges)
    if (allocated(self%centerParticles)) deallocate(self%centerParticles)
    if (allocated(self%children)) deallocate(self%children)
    if (allocated(self%hasChildren)) deallocate(self%hasChildren)
    if (allocated(self%parent)) deallocate(self%parent)
    if (allocated(self%area)) deallocate(self%area)
end subroutine

subroutine deleteQuad(self)
    type(QuadLinearFaces), intent(inout) :: self
    if (allocated(self%vertices)) deallocate(self%vertices)
    if (allocated(self%edges)) deallocate(self%edges)
    if (allocated(self%centerParticles)) deallocate(self%centerParticles)
    if (allocated(self%children)) deallocate(self%children)
    if (allocated(self%hasChildren)) deallocate(self%hasChildren)
    if (allocated(self%parent)) deallocate(self%parent)
    if (allocated(self%area)) deallocate(self%area)
end subroutine

subroutine deleteQuadCubic(self)
    type(QuadCubicFaces), intent(inout) :: self
    if (allocated(self%vertices)) deallocate(self%vertices)
    if (allocated(self%edges)) deallocate(self%edges)
    if (allocated(self%centerParticles)) deallocate(self%centerParticles)
    if (allocated(self%children)) deallocate(self%children)
    if (allocated(self%hasChildren)) deallocate(self%hasChildren)
    if (allocated(self%parent)) deallocate(self%parent)
    if (allocated(self%area)) deallocate(self%area)
end subroutine

subroutine copy(self, other)
    class(Faces), intent(inout) :: self
    class(Faces), intent(in) :: other

    if (self%N < other%N) then
        call LogMessage(log, ERROR_LOGGING_LEVEL, logKey, " copy error : not enough memory.")
        return
    endif
!    if (self%faceKind /= other%faceKind) then
!        call LogMessage(log, ERROR_LOGGING_LEVEL, logkey, " copy error : facekind mismatch.")
!        return
!    endif
    self%N = other%N
    self%N_Active = other%N_Active
    self%vertices(:,1:other%N) = other%vertices(:, 1:other%N)
    self%edges(:,1:other%N) = other%edges(:,1:other%N)
    self%centerParticles(:,1:other%N) = other%centerParticles(:,1:other%N)
    self%children(:,1:other%N) = other%children(:,1:other%N)
    self%hasChildren(1:other%N) = other%hasChildren(1:other%N)
    self%parent(1:other%N) = other%parent(1:other%N)
    self%area(1:other%N) = other%area(1:other%N)
end subroutine

subroutine logStats(self, alog)
    class(Faces), intent(in) :: self
    type(Logger), intent(inout) :: alog
    !
    call LogMessage(alog, TRACE_LOGGING_LEVEL, logkey, " Faces stats:")
    call StartSection(aLog)
    call LogMessage(alog, TRACE_LOGGING_LEVEL, "N = ", self%N)
    call LogMessage(alog, TRACE_LOGGING_LEVEL, "N_Active = ", self%N_Active)
    call LogMessage(alog, TRACE_LOGGING_LEVEL, "N_Max = ", self%N_Max)
    call LogMessage(aLog, TRACE_LOGGING_LEVEL, "total area = ", sum(self%area(1:self%N)))
    call EndSection(aLog)
end subroutine

subroutine insert(self, centerInds, vertInds, edgeInds, area)
    class(Faces), intent(inout) :: self
    integer(kint), dimension(:), intent(in) :: centerInds, vertInds, edgeInds
    real(kreal), intent(in), optional :: area

    if (self%N+1 > self%N_Max) then
        call LogMessage(log, ERROR_LOGGING_LEVEL, logkey, " insert error : not enough memory.")
        return
    endif
    if (size(centerInds) /= size(self%centerParticles,1) .or. size(vertInds) /= size(self%vertices,1) ) then
        call LogMessage(log, ERROR_LOGGING_LEVEL, logkey, " insert error : size mismatch.")
        return
    endif

    self%centerParticles(:,self%N+1) = centerInds
    self%vertices(:,self%N+1) = vertInds
    self%edges(:,self%N+1) = edgeInds
    if (present(area)) then
        self%area(self%N+1) = area
    endif
    self%N = self%N + 1
end subroutine

!pure function triCorners(self, index, aParticles)
!    real(kreal), dimension(:,:) :: triCorners
!    class(TriLinearFaces), intent(in) :: self
!    class(Particles), intent(in) :: aParticles
!    integer(kint), intent(in) :: index
!    !
!    integer(kint) :: i
!    do i=1,3
!        triCorners(:,i) = aParticles%physCoord(self%vertices(i,index))
!    enddo
!end function
!
!pure function quadCorners(self, index, aParticles)
!    real(kreal), dimension(:,:) :: quadCorners
!    class(quadLinearFaces), intent(in) :: self
!    class(Particles), intent(in) :: aParticles
!    integer(kint), intent(in) :: index
!    !
!    integer(kint) :: i
!    do i=1,4
!        quadCorners(:,i) = aParticles%physCoord(self%vertices(i,index))
!    enddo
!end function
!
!pure function quadCubicCorners(self, index, aParticles)
!    real(kreal), dimension(:,:) :: quadCubicCorners
!    class(QuadCubicFaces), intent(in) :: self
!    integer(kint), intent(in) :: index
!    class(Particles), intent(in) :: aParticles
!    !
!    integer(kint) :: i
!
!    do i=1,4
!        quadCubicCorners(:,i) = aParticles%physCoord(self%vertices(mod(3*i+9,12),index))
!    enddo
!end function

subroutine divideQuadLinear(self, index, aParticles, anEdges)
    class(QuadLinearFaces), intent(inout) :: self
    integer(kint), intent(in) :: index
    class(Particles), intent(inout) :: aParticles
    class(Edges), intent(inout) :: anEdges
    !
    integer(kint) :: i, j, newFaceVerts(4,4), newFaceEdges(4,4)
    integer(kint) :: parentEdge, childEdge1, childEdge2
    real(kreal) :: quadCoords(3,4), lagQuadCoords(3,4), newFaceCenters(3,4), newFaceLagCenters(3,4)

    if (self%N + 4 > self%N_Max) then
        call LogMessage(log, ERROR_LOGGING_LEVEL, " divideQuadLinear error : ", "not enough memory.")
        return
    endif

    newFaceVerts = 0
    newFaceEdges = 0
    !
    !   connect parent vertices to child faces
    !
    do i=1, 4
        newFaceVerts(i,i) = self%vertices(i,index)
    enddo

    !
    !   loop over parent edges
    !
    do i=1, 4
        parentEdge = self%edges(i,index)
        if (anEdges%hasChildren(parentEdge)) then
            !
            !   edges has already been divided by adjacent panel
            !
            childEdge1 = anEdges%child1(parentEdge)
            childEdge2 = anEdges%child2(parentEdge)
        else
            !
            !   divide parent edge
            !
            childEdge1 = anEdges%N+1
            childEdge2 = anEdges%N+2
            call anEdges%divide(parentEdge, aParticles)
        endif

        !
        !   connect child edges to new child faces
        !
        if (anEdges%positiveOrientation(parentEdge, index)) then
            newFaceEdges(i,i) = childEdge1
            anEdges%leftFace(childEdge1) = self%N + i

            newFaceEdges(i, mod(i,4)+1) = childEdge2
            anEdges%leftFace(childedge2) = self%N + mod(i,4) + 1
        else
            newFaceEdges(i,i) = childEdge2
            anEdges%rightFace(childEdge2) = self%N+i

            newFaceEdges(i, mod(i,4)+1) = childEdge1
            anEdges%rightFace(childEdge1) = self%N + mod(i,4) + 1
        endif

        newFaceVerts(mod(i,4)+1, i) = anEdges%dest(childEdge1)
        newFaceVerts(i, mod(i,4)+1) = anEdges%dest(childEdge1)
    enddo

    !
    !   change parent center particle to vertex of new child panels
    !
    do i=1,4
        newFaceVerts(mod(i+1,4)+1, i) = self%centerParticles(1,index)
    enddo

    !
	!	debugging : check vertex connectivity
	!
	do i = 1, 4
		do j = 1, 4
			if ( newFaceVerts(j,i) < 1 ) then
				write(logstring,*) " vertex connectivity ERROR at parent face ", index, ", child ", i, ", vertex ", j
				call LogMessage(log, ERROR_LOGGING_LEVEL,logkey//" DivideQuadFace :",logstring)
			endif
		enddo
	enddo

	!
	!   create new interior edges
	!
	call anEdges%insert(newFaceVerts(2,1), newFaceVerts(3,1), self%N+1, self%N+2)
	newFaceEdges(2,1) = anEdges%N
	newFaceEdges(4,2) = anEdges%N

	call anEdges%insert(newFaceVerts(1,3), newFaceVerts(4,3), self%N+4, self%N+3)
	newFaceEdges(4,3) = anEdges%N
	newFaceEdges(2,4) = anEdges%N

	call anEdges%insert(newFaceVerts(3,2), newFaceVerts(4,2), self%N+2, self%N+3)
	newFaceEdges(3,2) = anEdges%N
	newFaceEdges(1,3) = anEdges%N

	call anEdges%insert(newFaceVerts(2,4), newFaceVerts(1,4), self%N+1, self%N+4)
	newFaceEdges(1,4) = anEdges%N
	newFaceEdges(3,1) = anEdges%N

	!
	! debugging : check edge connectivity
	!
	do i = 1, 4
		do j = 1, 4
			if ( newFaceEdges(j,i) < 1 ) then
				write(logstring,*) " edge connectivity ERROR at parent face ", index, ", child ", i, ", vertex ", j
				call LogMessage(log, ERROR_LOGGING_LEVEL, logkey//" DivideQuadFace :",logstring)
			endif
		enddo
	enddo

	!
	!   create new particles for child face centers
	!
	quadCoords = dzero
	lagQuadCoords = dzero
	newFaceCenters = dzero
	newFaceLagCenters = dzero
	do i=1, 4
	    do j=1,4
	        quadCoords(:,j) = aParticles%physCoord(newFaceVerts(j,i))
	        lagQuadCoords(:,j) = aParticles%lagCoord(newFaceVerts(j,i))
	    enddo
	    if (aParticles%geomKind == PLANAR_GEOM) then
	        newFaceCenters(1:2,i) = QuadCentroid(quadCoords(1:2,1), quadCoords(1:2,2), quadCoords(1:2,3), quadCoords(1:2,4))
	        newFaceLagCenters(1:2,i) = QuadCentroid(lagQuadCoords(1:2,1), lagQuadCoords(1:2,2), lagQuadCoords(1:2,3), &
	                 lagQuadCoords(1:2,4))
	    elseif (aParticles%geomKind == SPHERE_GEOM) then
	        newFaceCenters(:,i) = SphereQuadCenter(quadCoords(:,1), quadCoords(:,2), quadCoords(:,3), quadCoords(:,4))
	        newFaceLagCenters(:,i) = SphereQuadCenter(lagQuadCoords(:,1), lagQuadCoords(:,1), lagQuadCoords(:,1), &
	             lagQuadCoords(:,4))
	    endif
	enddo

	!
	!   create child faces
	!
	do i=1, 4
	    call aParticles%insert(newFaceCenters(:,i), newFaceLagCenters(:,i))
	    self%centerParticles(1,self%N+i) = aParticles%N
	    self%vertices(:,self%N+i) = newFaceVerts(:,i)
	    self%edges(:,self%N+i) = newFaceEdges(:,i)
	    self%children(i, index) = self%N+i
	    self%parent(self%N+i) = index
	    aParticles%weight(self%centerParticles(1,self%N+i)) = self%setArea(self%N+i, aParticles)
	    self%area(self%N+i) = self%setArea(self%N+i,aParticles)
	enddo
	self%area(index) = dzero
	aParticles%weight(self%centerParticles(1,index)) = dzero
	self%hasChildren(index) = .TRUE.
	self%N = self%N+4
	self%N_Active = self%N_Active + 3
end subroutine



function getVerticesFromEdge(self, faceIndex, edgeIndex, anEdges)
    integer(kint), dimension(4) :: getVerticesFromEdge
    class(QuadCubicFaces), intent(in) :: self
    integer(kint), intent(in) :: faceIndex, edgeIndex
    class(CubicEdges), intent(in) :: anEdges

    if (anEdges%positiveOrientation(edgeIndex, faceIndex)) then
        getVerticesFromEdge(1) = anEdges%orig(edgeIndex)
        getVerticesFromEdge(2:3) = anEdges%interiorParticles(:,edgeIndex)
        getVerticesFromEdge(4) = anEdges%dest(edgeIndex)
    else
        getVerticesFromEdge(1) = anEdges%dest(edgeIndex)
        getVerticesFromEdge(2) = anEdges%interiorParticles(2,edgeIndex)
        getVerticesFromEdge(3) = anEdges%interiorParticles(1,edgeIndex)
        getVerticesFromEdge(4) = anEdges%orig(edgeIndex)
    endif
end function

subroutine divideQuadCubic(self, index, aParticles, anEdges)
    class(QuadCubicFaces), intent(inout) :: self
    integer(kint), intent(in) :: index
    class(Particles), intent(inout) :: aParticles
    class(Edges), intent(inout) :: anEdges
    !
    integer(kint) :: i, j, parentEdge, childEdge1, childEdge2
    integer(kint) :: nParticles, nEdges, nFaces
    integer(kint), dimension(12,4) :: newFaceVerts
    integer(kint), dimension(4,4) :: newFaceEdges, newFaceCenters
    real(kreal), dimension(3,4) :: physCorners, lagCorners, childCorners, lagChildCorners
    real(kreal), dimension(3,4) :: edgeMidpts, lagEdgeMidpts
    real(kreal), dimension(3) :: physCenter, lagCenter
    real(kreal), dimension(3,8) :: newIntEdgePts, lagNewIntEdgePts
    integer(kint), dimension(4) :: edgeMidptInds
    real(kreal), dimension(3,16) :: newFaceIntPts, lagNewFaceIntPts
    real(kreal), dimension(16) :: jac

    if ( self%N+4 > self%N_Max ) then
        write(logstring,'(2(A,I8))') 'faces.n = ', self%N, ', faces.N_Max = ', self%N_Max
		call LogMessage(log, ERROR_LOGGING_LEVEL, trim(logkey)//" DivideQuadFace ERROR : not enough memory.", logstring)
		return
	endif
	select type(anEdges)
	class is (CubicEdges)
	    newFaceEdges = 0
	    newFaceVerts = 0
	    newFaceCenters = 0
	    !
	    !   find corners and centroid of parent
	    !

        do i=1,4
            physCorners(:,i) = aParticles%physCoord(self%vertices(mod(3*i+9,12)+1,index))
            lagCorners(:,i) = aParticles%lagCoord(self%vertices(mod(3*i+9,12)+1, index))
        enddo
        if (aParticles%geomKind == SPHERE_GEOM) then
            physCenter = SphereQuadCenter(physCorners(:,1), physCorners(:,2), physCorners(:,3), physCorners(:,4))
            lagCenter = SphereQuadCenter(lagCorners(:,1), lagCorners(:,2), lagCorners(:,3), lagCorners(:,4))
        else
            physCenter = sum(physCorners,2) * 0.25_kreal
            lagCenter = sum(lagCorners,2) * 0.25_kreal
        endif
        !
        !   loop over parent edges
        !
        do i=1,4
            parentEdge = self%edges(i,index)
            if (anEdges%hasChildren(parentEdge)) then
                childEdge1 = anEdges%child1(parentEdge)
                childEdge2 = anEdges%child2(parentEdge)
            else
                childEdge1 = anEdges%N+1
                childEdge2 = anEdges%N+2
                call anEdges%divide(parentEdge, aParticles)
            endif
            edgeMidpts(:,i) = aParticles%physCoord(anEdges%dest(childEdge1))
            lagEdgeMidpts(:,i) = aParticles%lagCoord(anEdges%dest(childEdge1))
            edgeMidptInds(i) = anEdges%dest(childEdge1)

            if (anEdges%positiveOrientation(parentEdge, index)) then
                newFaceEdges(i,i) = childEdge1
                anEdges%leftFace(childEdge1) = self%N+i

                newFaceEdges(i, mod(i,4)+1) = childEdge2
                anEdges%leftFace(childedge2) = self%N + mod(i,4)+1
            else
                newFaceEdges(i,i) = childEdge2
                anEdges%rightFace(childedge2) = self%N+i

                newFaceEdges(i,mod(i,4)+1) = childEdge1
                anEdges%rightFace(childEdge1) = self%N+ mod(i,4)+1
            endif
        enddo
        !
        !   create interior edges
        !
        nEdges = anEdges%N
        nParticles = aParticles%N
        call aParticles%insert(physCenter, lagCenter) ! nParticles + 1
        newFaceVerts(7,1) = nParticles+1
        newFaceVerts(10,2)= nParticles+1
        newFaceVerts(1,3) = nParticles+1
        newFaceVerts(4,4) = nParticles+1
        if (aParticles%geomKind == SPHERE_GEOM) then
            newIntEdgePts(:,1) = pointAlongSphereVector(edgeMidpts(:,1), physCenter, gll_cubic_qp(2))
            newIntEdgePts(:,2) = pointAlongSphereVector(edgeMidpts(:,1), physCenter, gll_cubic_qp(3))
            newIntEdgePts(:,3) = pointAlongSphereVector(edgeMidpts(:,2), physCenter, gll_cubic_qp(2))
            newIntEdgePts(:,4) = pointAlongSphereVector(edgeMidpts(:,2), physCenter, gll_cubic_qp(3))
            newIntEdgePts(:,5) = pointAlongSphereVector(physCenter, edgeMidpts(:,3), gll_cubic_qp(2))
            newIntEdgePts(:,6) = pointAlongSphereVector(physCenter, edgeMidpts(:,3), gll_cubic_qp(3))
            newIntEdgePts(:,7) = pointAlongSphereVector(physCenter, edgeMidpts(:,4), gll_cubic_qp(2))
            newIntEdgePts(:,8) = pointAlongSphereVector(physCenter, edgeMidpts(:,4), gll_cubic_qp(3))

            lagNewIntEdgePts(:,1) = pointAlongSphereVector(lagEdgeMidPts(:,1), lagCenter, gll_cubic_qp(2))
            lagNewIntEdgePts(:,2) = pointAlongSphereVector(lagEdgeMidPts(:,1), lagCenter, gll_cubic_qp(3))
            lagNewIntEdgePts(:,3) = pointAlongSphereVector(lagEdgeMidPts(:,2), lagCenter, gll_cubic_qp(2))
            lagNewIntEdgePts(:,4) = pointAlongSphereVector(lagEdgeMidPts(:,2), lagCenter, gll_cubic_qp(3))
            lagNewIntEdgePts(:,5) = pointAlongSphereVector(lagCenter, lagEdgeMidPts(:,3), gll_cubic_qp(2))
            lagNewIntEdgePts(:,6) = pointAlongSphereVector(lagCenter, lagEdgeMidPts(:,3), gll_cubic_qp(3))
            lagNewIntEdgePts(:,7) = pointAlongSphereVector(lagCenter, lagEdgeMidPts(:,4), gll_cubic_qp(2))
            lagNewIntEdgePts(:,8) = pointAlongSphereVector(lagCenter, lagEdgeMidPts(:,4), gll_cubic_qp(3))
        else
            newIntEdgePts(:,1) = pointAlongChordVector(edgeMidpts(:,1), physCenter, gll_cubic_qp(2))
            newIntEdgePts(:,2) = pointAlongChordVector(edgeMidpts(:,1), physCenter, gll_cubic_qp(3))
            newIntEdgePts(:,3) = pointAlongChordVector(edgeMidpts(:,2), physCenter, gll_cubic_qp(2))
            newIntEdgePts(:,4) = pointAlongChordVector(edgeMidpts(:,2), physCenter, gll_cubic_qp(3))
            newIntEdgePts(:,5) = pointAlongChordVector(physCenter, edgeMidpts(:,3), gll_cubic_qp(2))
            newIntEdgePts(:,6) = pointAlongChordVector(physCenter, edgeMidpts(:,3), gll_cubic_qp(3))
            newIntEdgePts(:,7) = pointAlongChordVector(physCenter, edgeMidpts(:,4), gll_cubic_qp(2))
            newIntEdgePts(:,8) = pointAlongChordVector(physCenter, edgeMidpts(:,4), gll_cubic_qp(3))

            lagNewIntEdgePts(:,1) = pointAlongChordVector(lagEdgeMidPts(:,1), lagCenter, gll_cubic_qp(2))
            lagNewIntEdgePts(:,2) = pointAlongChordVector(lagEdgeMidPts(:,1), lagCenter, gll_cubic_qp(3))
            lagNewIntEdgePts(:,3) = pointAlongChordVector(lagEdgeMidPts(:,2), lagCenter, gll_cubic_qp(2))
            lagNewIntEdgePts(:,4) = pointAlongChordVector(lagEdgeMidPts(:,2), lagCenter, gll_cubic_qp(3))
            lagNewIntEdgePts(:,5) = pointAlongChordVector(lagCenter, lagEdgeMidPts(:,3), gll_cubic_qp(2))
            lagNewIntEdgePts(:,6) = pointAlongChordVector(lagCenter, lagEdgeMidPts(:,3), gll_cubic_qp(3))
            lagNewIntEdgePts(:,7) = pointAlongChordVector(lagCenter, lagEdgeMidPts(:,4), gll_cubic_qp(2))
            lagNewIntEdgePts(:,8) = pointAlongChordVector(lagCenter, lagEdgeMidPts(:,4), gll_cubic_qp(3))
        endif
        do i=1,8
            call aParticles%insert(newIntEdgePts(:,i), lagNewIntEdgePts(:,i)) ! nParticles + 1 + i
        enddo
        call anEdges%insert(edgeMidptInds(1), nParticles+1, self%N+1, self%N+2, [nParticles+2, nParticles+3]) ! nedges +1
        newFaceEdges(2,1) = nEdges+1
        newFaceEdges(4,2) = nEdges+1
        newFaceVerts(5,1) = nParticles+2
        newFaceVerts(12,2)= nParticles+2
        newFaceVerts(6,1) = nParticles+3
        newFaceVerts(11,2)= nParticles+3
        call anEdges%insert(edgeMidptInds(2), nParticles+1, self%N+2, self%N+3, [nParticles+4, nParticles+5])! nedges +2
        newFaceEdges(3,2) = nEdges+2
        newFaceEdges(1,3) = nEdges+2
        newFaceVerts(8,2) = nParticles+4
        newFaceVerts(3,3) = nParticles+4
        newFaceVerts(9,2) = nParticles+5
        newFaceVerts(2,3) = nParticles+5
        call anEdges%insert(nParticles+1, edgeMidptInds(3), self%N+4, self%N+3, [nParticles+6, nParticles+7])! nedges +3
        newFaceEdges(2,4) = nedges+3
        newFaceEdges(4,3) = nEdges+3
        newFaceVerts(12,3) = nParticles+6
        newFaceVerts(5,4) = nParticles+6
        newFaceVerts(11,3) = nParticles+7
        newFaceVerts(6,4) = nParticles+7
        call anEdges%insert(nParticles+1, edgeMidptInds(4), self%N+1, self%N+4, [nParticles+8, nParticles+9])! nedges +4
        newFaceEdges(3,1) = nedges+4
        newFaceEdges(1,4) = nedges+4
        newFaceVerts(3,4) = nParticles + 8
        newFaceVerts(8,1) = nParticles+8
        newFaceVerts(2,4) = nParticles+9
        newFaceVerts(9,1) = nParticles+9

        !
        !   locate face interior particles
        !
        childCorners(:,1) = physCorners(:,1)
        childCorners(:,2) = edgeMidpts(:,1)
        childCorners(:,3) = physCenter
        childCorners(:,4) = edgeMidpts(:,4)
        lagChildCorners(:,1) = lagCorners(:,1)
        lagChildCorners(:,2) = lagEdgeMidpts(:,1)
        lagChildCorners(:,3) = lagCenter
        lagChildCorners(:,4) = lagEdgeMidpts(:,4)
        newFaceIntPts(:,1:4) = quad16InteriorPts(childCorners)
        lagNewFaceIntPts(:,1:4) = quad16InteriorPts(lagChildCorners)

        childCorners(:,1) = edgeMidpts(:,1)
        childCorners(:,2) = physCorners(:,2)
        childCorners(:,3) = edgeMidpts(:,2)
        childCorners(:,4) = physCenter
        lagChildCorners(:,1) = lagEdgeMidpts(:,1)
        lagChildCorners(:,2) = lagCorners(:,2)
        lagChildCorners(:,3) = lagEdgeMidpts(:,2)
        lagChildCorners(:,4) = lagCenter
        newFaceIntPts(:,5:8) = quad16InteriorPts(childCorners)
        lagNewFaceIntPts(:,5:8) = quad16InteriorPts(lagChildCorners)

        childCorners(:,1) = physCenter
        childCorners(:,2) = edgeMidpts(:,2)
        childCorners(:,3) = physCorners(:,3)
        childCorners(:,4) = edgeMidpts(:,3)
        lagChildCorners(:,1) = lagCenter
        lagChildCorners(:,2) = lagEdgeMidpts(:,2)
        lagChildCorners(:,3) = lagCorners(:,3)
        lagChildCorners(:,4) = lagEdgeMidpts(:,3)
        newFaceIntPts(:,9:12) = quad16InteriorPts(childCorners)
        lagNewFaceIntPts(:,9:12) = quad16InteriorPts(lagChildCorners)

        childCorners(:,1) = edgeMidpts(:,4)
        childCorners(:,2) = physCenter
        childCorners(:,3) = edgeMidpts(:,3)
        childCorners(:,4) = physCorners(:,4)
        lagChildCorners(:,1) = lagEdgeMidpts(:,4)
        lagChildCorners(:,2) = lagCenter
        lagChildCorners(:,3) = lagEdgeMidpts(:,3)
        lagChildCorners(:,4) = lagCorners(:,4)
        newFaceIntPts(:,13:16) = quad16InteriorPts(childCorners)
        lagNewFaceIntPts(:,13:16) = quad16InteriorPts(lagChildCorners)

        ! move parent particles to new locations in child faces
        call aParticles%replace(self%centerParticles(1,index), newFaceIntPts(:,1), lagNewFaceIntPts(:,1))
        call aParticles%replace(self%centerParticles(2,index), newFaceIntPts(:,6), lagNewFaceIntPts(:,6))
        call aParticles%replace(self%centerParticles(3,index), newFaceIntPts(:,11), lagNewFaceIntPts(:,11))
        call aParticles%replace(self%centerParticles(4,index), newFaceIntPts(:,16), lagNewFaceIntPts(:,16))
        do i=1,4
            newFaceCenters(i,i) = self%centerParticles(i,index)
        enddo
        ! insert new interior particles for child faces
        nParticles = aParticles%N
        do i=1,4
            if (i /= 1 ) then
                call aParticles%insert(newFaceIntPts(:,i), lagNewFaceIntPts(:,i))
                newFaceCenters(i,1) = aParticles%N
            endif
        enddo
        do i=1,4
            if (i/=2) then
                call aParticles%insert(newFaceIntPts(:,i+4), lagNewFaceIntPts(:,i+4))
                newFaceCenters(i,2) = aParticles%N
            endif
        enddo
        do i=1,4
            if (i/=3) then
                call aParticles%insert(newFaceIntPts(:,i+8), lagNewFaceIntPts(:,i+8))
                newFaceCenters(i,3) = aParticles%N
            endif
        enddo
        do i=1,4
            if (i/=4) then
                call aParticles%insert(newFaceIntPts(:,i+12), lagNewFaceIntPts(:,i+12))
                newFaceCenters(i,4) = aParticles%N
            endif
        enddo
        !
        !   set face vertices
        !
        if (anEdges%positiveOrientation(self%edges(1,index), index)) then
            newFaceVerts(1,1) = anEdges%orig(newFaceEdges(1,1))
            newFaceVerts(2,1) = anEdges%interiorParticles(1,newFaceEdges(1,1))
            newFaceVerts(3,1) = anEdges%interiorParticles(2,newFaceCenters(1,1))
            newFaceVerts(4,1) = anEdges%dest(newFaceEdges(1,1))
            newFaceVerts(1,2) = anEdges%dest(newFaceEdges(1,1))

            newFaceVerts(2,2) = anEdges%interiorParticles(1,newFaceEdges(1,2))
            newFaceVerts(3,2) = anEdges%interiorParticles(2,newFaceEdges(1,2))
            newFaceVerts(4,2) = anEdges%dest(newFaceEdges(1,2))
        else
            newFaceVerts(1,1) = anEdges%dest(newFaceEdges(1,1))
            newFaceVerts(2,1) = anEdges%interiorParticles(2,newFaceEdges(1,1))
            newFaceVerts(3,1) = anEdges%interiorParticles(1,newFaceEdges(1,1))
            newFaceVerts(4,1) = anEdges%orig(newFaceEdges(1,1))
            newFaceVerts(1,2) = anEdges%orig(newFaceEdges(1,1))

            newFaceVerts(2,2) = anEdges%interiorParticles(2,newFaceEdges(1,2))
            newFaceVerts(3,2) = anEdges%interiorParticles(1,newFaceEdges(1,2))
            newFaceVerts(4,2) = anEdges%orig(newFaceEdges(1,2))
        endif
        if (anEdges%positiveOrientation(self%edges(2,index), index)) then
            newFaceVerts(5,2) = anEdges%interiorParticles(1,newFaceEdges(2,2))
            newFaceVerts(6,2) = anEdges%interiorParticles(2,newFaceEdges(2,2))
            newFaceVerts(7,2) = anEdges%dest(newFaceEdges(2,2))
            newFaceVerts(4,3) = anEdges%dest(newFaceEdges(2,2))

            newFaceVerts(5,3) = anEdges%interiorParticles(1,newFaceEdges(2,3))
            newFaceVerts(6,3) = anEdges%interiorParticles(2,newFaceEdges(2,3))
            newFaceVerts(7,3) = anEdges%dest(newFaceEdges(2,3))
        else
            newFaceVerts(5,2) = anEdges%interiorParticles(2,newFaceEdges(2,2))
            newFaceVerts(6,2) = anEdges%interiorParticles(1,newFaceEdges(2,2))
            newFaceVerts(7,2) = anEdges%orig(newFaceEdges(2,2))
            newFaceVerts(4,3) = anEdges%orig(newFaceEdges(2,2))

            newFaceVerts(5,3) = anEdges%interiorParticles(2,newFaceEdges(2,3))
            newFaceVerts(6,3) = anEdges%interiorParticles(2,newFaceEdges(2,3))
            newFaceVerts(7,3) = anEdges%orig(newFaceEdges(2,3))
        endif
        if (anEdges%positiveOrientation(self%edges(3,index), index)) then
            newFaceVerts(8,3) = anEdges%interiorParticles(1,newFaceEdges(3,3))
            newFaceVerts(9,3) = anEdges%interiorParticles(2,newFaceEdges(3,3))
            newFaceVerts(10,3)= anEdges%dest(newFaceEdges(3,3))
            newFaceVerts(7,4) = anEdges%dest(newFaceEdges(3,3))

            newFaceVerts(8,4) = anEdges%interiorParticles(1,newFaceEdges(3,4))
            newFaceVerts(9,4) = anEdges%interiorParticles(2,newFaceEdges(3,4))
            newFaceVerts(10,4)= anEdges%dest(newFaceEdges(3,4))
        else
            newFaceVerts(8,3) = anEdges%interiorParticles(2,newFaceEdges(3,3))
            newFaceVerts(9,3) = anEdges%interiorParticles(1,newFaceEdges(3,3))
            newFaceVerts(10,3)= anEdges%orig(newFaceEdges(3,3))
            newFaceVerts(7,4) = anEdges%orig(newFaceEdges(3,3))

            newFaceVerts(8,4) = anEdges%interiorParticles(2,newFaceEdges(3,4))
            newFaceVerts(9,4) = anEdges%interiorParticles(1,newFaceEdges(3,4))
            newFaceVerts(10,4)= anEdges%orig(newFaceEdges(3,4))
        endif
        if (anEdges%positiveOrientation(self%edges(4,index), index)) then
            newFaceVerts(11,4) = anEdges%interiorParticles(1,newFaceEdges(4,4))
            newFaceVerts(12,4) = anEdges%interiorParticles(2,newFaceEdges(4,4))
            newFaceVerts(1,4) = anEdges%dest(newFaceEdges(4,4))
            newFaceVerts(10,1) = anEdges%dest(newFaceEdges(4,4))

            newFaceVerts(11,1) = anEdges%interiorParticles(1,newFaceEdges(4,1))
            newFaceVerts(12,1) = anEdges%interiorParticles(2,newFaceEdges(4,1))
        else
            newFaceVerts(11,4) = anEdges%interiorParticles(2,newFaceEdges(4,4))
            newFaceVerts(12,4) = anEdges%interiorParticles(1,newFaceEdges(4,4))
            newFaceVerts(1,4) = anEdges%orig(newFaceEdges(4,4))
            newFaceVerts(10,1) = anEdges%orig(newFaceEdges(4,4))

            newFaceVerts(11,1) = anEdges%interiorParticles(2,newFaceEdges(4,1))
            newFaceVerts(12,1) = anEdges%interiorParticles(1,newFaceEdges(4,1))
        endif

        !check connectivity
        do i=1,4
            do j=1,4
                if (newFaceEdges(j,i) == 0) then
                    write(logstring,'(3(A,I4))') 'divideCubicFace connectivity error at parent face ', &
                        index, '; child face ', i, ' edge ', j
                    call LogMessage(log, ERROR_LOGGING_LEVEL, trim(logkey)//" ", trim(logstring))
                endif
                if (newFaceCenters(j,i) == 0) then
                    write(logstring,'(3(A,I4))') 'divideCubicFace connectivity error at parent face ',&
                         index, '; child face ', i, ' centerParticle', j
                    call LogMessage(log, ERROR_LOGGING_LEVEL, trim(logkey)//" ", trim(logstring))
                endif
            enddo
            do j=1,12
                if (newFaceVerts(j,i) == 0) then
                    write(logString,'(3(A,I4))') 'divideCubicFace connectivity error at parent face ', &
                        index, '; child face ', i, ' vertex ', j
                    call LogMessage(log, ERROR_LOGGING_LEVEL, trim(logkey)//" ", trim(logstring))
                endif
            enddo
        enddo

!        print *, 'inserting faces'
!        do i=1,4
!            write(6, '(A,I2,A,I2,A,4(I4))') 'parent ', index, ', child ', i, ' centerInds ', newFaceCenters
!            write(6, '(A, 12(I4))') '...face verts', newFaceVerts(:,i)
!            write(6,'(A,4(I4))') '...face edges', newFaceEdges
!            do j=1,4
!                write(6,'(A,I2, A, 4I4)') '.../...', j, ': ', anEdges%orig(newFaceEdges(j,i)), &
!                    anEdges%interiorParticles(:,newFaceEdges(j,i)), &
!                    anEdges%dest(newFaceEdges(j,i))
!            enddo
!        enddo


        nFaces = self%n
        do i=1,4
            call self%insert(newFaceCenters(:,i), newFaceVerts(:,i), newFaceEdges(:,i))
            self%children(i,index) = self%N
            self%parent(self%N) = index
        enddo
        self%hasChildren(index) = .TRUE.
        self%area(index) = dzero


!        call self%writeMatlab(6)
!        call anEdges%writeMatlab(6)

        do i=1,4
            self%area(nFaces+i) = self%setArea(nFaces+i, aParticles)
        enddo

        do i=1, 4
            do j=1,4
                physCorners(:,j) = aParticles%physCoord(newFaceVerts(mod(3*j+9,12)+1,i))
            enddo
            do j=1,12
                jac(j) = bilinearPlaneJacobian(physCorners, quad16_vertex_qp(1,j), quad16_vertex_qp(2,j))
            enddo
            do j=1,4
                jac(12+j) = bilinearPlaneJacobian(physCorners, quad16_center_qp(1,j), quad16_center_qp(2,j))
            enddo
            do j=1,12
                aparticles%weight(newFaceVerts(j,i)) = quad16_vertex_qw(j) * jac(j)
            enddo
            do j=1,4
                aparticles%weight(newFaceCenters(j,i)) = quad16_center_qw(j)*jac(12+j)
            enddo
        enddo
        self%N_Active = self%N_Active +3

    class default
        call logMessage(log, ERROR_LOGGING_LEVEL, trim(logkey)//" divideCubicFace error : ", " cubic edges required.")
        return
    end select
end subroutine

pure function calcInteriorPts(self, vertXyz)
    class(QuadCubicFaces), intent(in) :: self
    real(kreal), dimension(3,4) :: calcInteriorPts
    real(kreal), dimension(3,4), intent(in) :: vertXyz
    calcInteriorPts(:,1) = bilinearMap(vertXyz, -oor5, oor5)
    calcInteriorPts(:,2) = bilinearMap(vertXyz, -oor5, -oor5)
    calcInteriorPts(:,3) = bilinearMap(vertXyz, oor5, -oor5)
    calcInteriorPts(:,4) = bilinearMap(vertXyz, oor5, oor5)
end function

pure function physCentroid(self, index, aParticles)
    real(kreal), dimension(3) :: physCentroid
    class(Faces), intent(in) :: self
    integer(kint), intent(in) :: index
    class(Particles), intent(in) :: aParticles
    !
    integer(kint) :: i
    real(kreal) :: norm

    physCentroid = dzero
    select type(self)
        class is (TriLinearFaces)
            do i=1, 3
                physCentroid = physCentroid + aParticles%physCoord(self%vertices(i, index))
            enddo
            physCentroid = physCentroid / real(3, kreal)
        class is (QuadLinearFaces)
            do i=1, 4
                physCentroid = physCentroid + aParticles%physCoord(self%vertices(i, index))
            enddo
            physCentroid = physCentroid / real(4, kreal)
        class is (QuadCubicFaces)
            do i=1,4
                physCentroid = physCentroid + aParticles%physCoord(self%vertices(mod(3*i+9,12)+1,index))
            enddo
            physCentroid = 0.25_kreal * physCentroid
    end select

    if (aParticles%geomKind == SPHERE_GEOM) then
        norm = sqrt(sum(physCentroid*physCentroid))
        physCentroid = physCentroid / norm
    endif
end function

pure function lagCentroid(self, index, aParticles)
    real(kreal), dimension(3) :: lagCentroid
    class(Faces), intent(in) :: self
    integer(kint), intent(in) :: index
    class(Particles), intent(in) :: aParticles
    !
    integer(kint) :: i, nv, nc
    real(kreal) :: norm

    lagCentroid = dzero
    nv = size(self%vertices,1)
    nc = size(self%centerParticles,1)
    do i=1, nv
        lagCentroid = lagCentroid + aParticles%lagCoord(self%vertices(i, index))
    enddo
    do i=1,nc
        lagCentroid = lagCentroid + aParticles%lagCoord(self%centerParticles(i, index))
    enddo
    lagCentroid = lagCentroid / real(nv + nc, kreal)
    if (aParticles%geomKind == SPHERE_GEOM) then
        norm = sqrt(sum(lagCentroid*lagCentroid))
        lagCentroid = lagCentroid / norm
    endif
end function

function countParents(self, index)
    integer(kint) :: countParents
    class(Faces), intent(in) :: self
    integer(kint), intent(in) :: index
    !
    logical(klog) :: keepGoing
    integer(kint) :: parentIndex

    countParents = 0
    keepGoing = (self%parent(index) > 0)
    parentIndex = self%parent(index)
    do while (keepGOing)
        countParents = countParents + 1
        parentIndex = self%parent(parentIndex)
        keepGoing = (parentIndex > 0)
    enddo
end function

function area(self, index, aParticles, anEdges)
    real(kreal) :: area
    class(Faces), intent(in) :: self
    integer(kint), intent(in) :: index
    class(Particles), intent(in) :: aParticles
    class(Edges), intent(in) :: anEdges
    !
    integer(kint) :: i, nEdges
    type(STDIntVector) :: leafEdges
    real(kreal), dimension(3) :: cntd

    nEdges = size(self%edges,1)
    area = dzero
    cntd = self%physCentroid(index, aParticles)
    do i=1, nEdges
        leafEdges = anEdges%getLeafEdges(self%edges(i, index))
        area = area + anEdges%areaFromLeaves(aParticles, cntd, leafEdges%integers)
    enddo
end function

function sharedEdge(self, face1, face2)
    integer(kint) :: sharedEdge
    class(Faces), intent(in) :: self
    integer(kint), intent(in) :: face1, face2
    !
    integer(kint) :: i, j, shareCount, nEdges

    sharedEdge = 0
    shareCount = 0
    nEdges = size(self%edges,1)

    do i=1, nEdges
        do j=1, nEdges
            if (self%edges(i, face1) == self%edges(j, face2)) then
                shareCount = shareCount + 1
                sharedEdge = self%edges(j,face2)
            endif
        enddo
    enddo
    if (shareCount>1) then
        call LogMessage(log, ERROR_LOGGING_LEVEL, trim(logKey)//" sharedEdge ERROR at face ", face1)
        sharedEdge = 0
    endif
end function

subroutine divideTri(self, index, aParticles, anEdges)
    class(TriLinearFaces), intent(inout) :: self
    integer(kint), intent(in) :: index
    class(Particles), intent(inout) :: aParticles
    class(Edges), intent(inout) :: anEdges
    !
    integer(kint) :: i, j, newFaceVerts(3,4), newFaceEdges(3,4)
    integer(kint) :: parentEdge, childEdge1, childEdge2
    real(kreal) :: triCoords(3,3), lagTriCoords(3,3), newFaceCenters(3,4), newFaceLagCenters(3,4)

!    if (self%faceKind /= TRI_PANEL) then
!         call LogMessage(log, ERROR_LOGGING_LEVEL, logkey, " divide tri face error : invalid facekind.")
!         return
!    endif
    if (self%N + 4 > self%N_Max) then
        call LogMessage(log, ERROR_LOGGING_LEVEL, logkey, " divide tri face error : not enough memory.")
        return
    endif

    newFaceVerts = 0
    newFaceEdges = 0
    !
    !   connect parent vertices to child faces
    !
    do i=1,3
        newFaceVerts(i,i) = self%vertices(i, index)
    enddo

    !
    !   loop over parent edges
    !
    do i=1,3
        parentEdge = self%edges(i, index)
        if (anEdges%hasChildren(parentEdge)) then
            !
            !   parent edge was already divided by adjacent panel
            !
            childEdge1 = anEdges%child1(parentEdge)
            childedge2 = anEdges%child2(parentEdge)
        else
            !
            !   divide parent edge
            !
            childEdge1 = anEdges%N + 1
            childEdge2 = anEdges%N + 2
            call anEdges%divide(parentEdge, aParticles)
        endif

        !
        !   connect child edges to child faces
        !
        if (anEdges%positiveOrientation(parentEdge, index)) then
            newFaceEdges(i, i) = childEdge1
            anEdges%leftFace(childEdge1) = self%N + i

            newFaceEdges(i, mod(i,3)+1) = childEdge2
            anEdges%leftFace(childedge2) = self%N + mod(i,3) + 1
        else
            newFaceEdges(i,i) = childEdge2
            anEdges%rightFace(childEdge2) = self%N + i

            newFaceEdges(i, mod(i,3) + 1) = childEdge1
            anEdges%rightFace(childEdge1) = self%N + mod(i,3) + 1
        endif

        newFaceVerts(i, mod(i,3) + 1) = anEdges%dest(childEdge1)
        newFaceVerts(mod(i,3)+1, i) = anEdges%dest(childEdge1)
    enddo

    newFaceVerts(:,4) = (/ newFaceVerts(3,2), newFaceVerts(1,3), newFaceVerts(2,1) /)

    !
	! debugging : check vertex connectivity
	!
	do i = 1, 4
		do j = 1, 3
			if ( newFaceVerts(j,i) < 1 ) then
				write(logstring,*) " vertex connectivity ERROR at parent face ", index, ", child ", i, ", vertex ",j
				call LogMessage(log,ERROR_LOGGING_LEVEL,logkey//" DivideTriFace :", logString )
			endif
		enddo
	enddo

	!
	!   create new interior edges
	!
	call anEdges%insert(newFaceVerts(1,4), newFaceVerts(2,4), self%N + 4, self%N+3)
	newFaceEdges(1,4) = anEdges%N
	newFaceEdges(1,3) = anEdges%N

	call anEdges%insert(newFaceVerts(2,4), newFaceVerts(3,4), self%N+4, self%N+1)
	newFaceEdges(2,4) = anEdges%N
	newFaceEdges(2,1) = anEdges%N

	call anEdges%insert(newFaceVerts(3,4), newFaceVerts(1,4), self%N+4, self%N+2)
	newFaceEdges(3,4) = anEdges%N
	newFaceEdges(3,2) = anEdges%N

	!
	! debugging : check edge connectivity
	!
	do i = 1, 4
		do j = 1, 3
			if ( newFaceEdges(j,i) < 1 ) then
				write(logstring,*) " edge connectivity ERROR at parent face ", index, ", child ", i, ", vertex ", j
				call LogMessage(log, ERROR_LOGGING_LEVEL, logkey//" DivideTriFace :", logString)
			endif
		enddo
	enddo

	!
	!   create new center particles for child faces 1:3
	!
	triCoords = dzero
	lagTriCoords = dzero
	newFaceCenters = dzero
	newFaceLagCenters = dzero
	do i=1,4
	    do j=1,3
	        triCoords(:,j) = aParticles%physCoord(newFaceVerts(j,i))
	        lagTriCoords(:,j) = aParticles%lagCoord(newFaceVerts(j,i))
	    enddo
	    if (aParticles%geomKind == PLANAR_GEOM) then
	        newFaceCenters(1:2,i) = TriCentroid(triCoords(1:2,1), triCoords(1:2,2), triCoords(1:2,3))
	        newFaceLagCenters(1:2,i) = TriCentroid(lagTriCoords(1:2,1), lagTriCoords(1:2,2), lagTriCoords(1:2,3))
	    elseif (aParticles%geomKind == SPHERE_GEOM) then
	        newFaceCenters(:,i) = SphereTriCenter(triCoords(:,1), triCoords(:,2), triCoords(:,3))
	        newFaceLagCenters(:,i) = SphereTriCenter(lagTriCoords(:,1), lagTriCoords(:,2), lagTriCoords(:,3))
	    endif
	enddo

	!
	!   create child faces
	!
	do i=1,3
	    call aParticles%insert(newFaceCenters(:,i), newFaceLagCenters(:,i))
	    self%centerParticles(1,self%N+i) = aParticles%N
	    self%vertices(:,self%N+i) = newFaceVerts(:,i)
	    self%edges(:, self%N+i) = newFaceEdges(:,i)
	    self%children(i, index) = self%N+i
	    self%parent(self%N+i) = index
	    self%area(self%N+i) = self%setArea(self%N+i,aParticles)
	    aParticles%weight(self%centerParticles(1,self%N+i)) = self%setArea(self%N+i, aParticles)
	enddo
	! special case for child 4: re-use parent face centerParticle
	self%centerParticles(1,self%N+4) = self%centerParticles(1,index)
	self%vertices(:,self%N+4) = newFaceVerts(:,4)
	self%edges(:,self%N+4) = newFaceEdges(:,4)
	self%children(4,index) = self%N+4
	self%parent(self%N+4) = index
	self%area(self%N+4) = self%setArea(self%N+4, aParticles)
	aParticles%weight(self%centerParticles(1,self%N+4)) = self%setArea(self%N+4, aParticles)

    self%area(index) = dzero
	self%hasChildren(index) = .TRUE.
	self%N = self%N + 4
	self%N_Active = self%N_Active + 3
end subroutine

pure function triFaceArea(self, index, aParticles)
    real(kreal) :: triFaceArea
    class(TriLinearFaces), intent(in) :: self
    integer(kint), intent(in) :: index
    class(Particles), intent(in) ::  aParticles
    !
    integer(kint) :: i
    real(kreal) :: center(3), v1(3), v2(3)

    triFaceArea = dzero
    center = aParticles%physCoord(self%centerParticles(1,index))
    if (aParticles%geomKind == PLANAR_GEOM) then
        do i=1, 3
            v1 = aParticles%physCoord(self%vertices(i, index))
            v2 = aParticles%physCoord(self%vertices(mod(i,3)+1, index))
            triFaceArea = triFaceArea + TriArea(v1(1:2), center(1:2), v2(1:2))
        enddo
    elseif (aParticles%geomKind == SPHERE_GEOM) then
        do i=1,3
            v1 = aParticles%physCoord(self%vertices(i, index))
            v2 = aParticles%physCoord(self%vertices(mod(i,3)+1, index))
            triFaceArea = triFaceArea + SphereTriArea(v1, center, v2)
        enddo
    endif
end function

pure function quadCubicArea(self, index, aParticles)
    real(kreal) :: quadCubicArea
    class(QuadCubicFaces), intent(in) :: self
    integer(kint), intent(in) :: index
    class(Particles), intent(in) :: aParticles
    !
    integer(kint) :: i
    real(kreal) :: center(3), verts(3,4), v1(3), v2(3)

    quadCubicArea = dzero
    center = self%physCentroid(index, aParticles)

    verts(:,1) = aParticles%physCoord(self%vertices(1,index))
    verts(:,2) = aParticles%physCoord(self%vertices(4,index))
    verts(:,3) = aParticles%physCoord(self%vertices(7,index))
    verts(:,4) = aParticles%physCoord(self%vertices(10,index))
    if (aParticles%geomKind == PLANAR_GEOM) then
        do i=1,4
            v1 = verts(:,i)
            v2 = verts(:,mod(i,4)+ 1)
            quadCubicArea = quadCubicArea + TriArea(v1(1:2), center(1:2), v2(1:2))
        enddo
    elseif (aParticles%geomKind == SPHERE_GEOM) then
        do i=1,4
            v1 = verts(:,i)
            v2 = verts(:,mod(i,4)+ 1)
            quadCubicArea = quadCubicArea + SphereTriArea(v1, center, v2)
        enddo
    endif
end function

pure function quadFaceArea(self, index, aParticles)
    real(kreal) :: quadFaceArea
    class(QuadLinearFaces), intent(in) :: self
    integer(kint), intent(in) :: index
    class(Particles), intent(in) :: aParticles
    !
    integer(kint) :: i
    real(kreal) :: center(3), v1(3), v2(3)

    quadFaceArea = dzero
    center = aParticles%physCoord(self%centerParticles(1,index))
    if (aParticles%geomKind == PLANAR_GEOM) then
        do i=1,4
            v1 = aParticles%physCoord(self%vertices(i,index))
            v2 = aParticles%physCoord(self%vertices(mod(i,4)+1, index))
            quadFaceArea = quadFaceArea + TriArea(v1(1:2), center(1:2), v2(1:2))
        enddo
    elseif (aParticles%geomKind == SPHERE_GEOM) then
        do i=1,4
            v1 = aParticles%physCoord(self%vertices(i,index))
            v2 = aParticles%physCoord(self%vertices(mod(i,4)+1, index))
            quadFaceArea = quadFaceArea + SphereTriArea(v1, center, v2)
        enddo
    endif
end function

function particleOppositeTriEdge(self, faceIndex, edgeIndex)
    integer(kint) :: particleOppositeTriEdge
    class(TriLinearFaces), intent(in) :: self
    integer(kint), intent(in) :: faceIndex, edgeIndex
    !
    integer(kint) :: i, foundEdge

    foundEdge = 0
    particleOppositeTriEdge = 0
    do i=1,3
        if (self%edges(i,faceIndex) == edgeIndex) then
            foundEdge = i
            exit
        endif
    enddo
    if (foundEdge > 0) then
        particleOppositeTriEdge = self%vertices(mod(foundEdge+1,3) + 1, faceIndex)
    else
        call LogMessage(log, ERROR_LOGGING_LEVEL, trim(logKey)//" particleOppositeTriEdge error : ", " edge not found.")
    endif
end function

subroutine writeMatlab(self, fileunit)
    class(Faces), intent(in) :: self
    integer(kint), intent(in) :: fileunit
    !
    integer(kint) :: i, j, nv, nc, ne

    nv = size(self%vertices,1)
    nc = size(self%centerParticles,1)
    ne = size(self%edges,1)

    write(fileunit,'(A)',advance='no') 'face_verts = ['
    do i=1, self%N-1
        do j=1, nv-1
            write(fileunit, '(I8,A)', advance='no') self%vertices(j,i), ','
        enddo
        write(fileunit,*) self%vertices(nv,i), ';'
    enddo
    do j=1, nv-1
        write(fileunit,'(I8,A)',advance='no') self%vertices(j,self%N),','
    enddo
    write(fileunit, *) self%vertices(nv,self%N), '];'

    write(fileunit,'(A)', advance='no') 'face_edges = ['
    do i=1, self%N-1
        do j=1, ne-1
            write(fileunit, '(I8,A)', advance='no') self%edges(j,i), ','
        enddo
        write(fileunit,*) self%edges(ne,i), ';'
    enddo
    do j=1, ne-1
        write(fileunit,'(I8,A)',advance='no') self%edges(j,self%N),','
    enddo
    write(fileunit, *) self%edges(ne,self%N), '];'

    write(fileunit,'(A)', advance='no') 'face_centers = ['
    do i=1, self%N-1
        do j=1, nc-1
            write(fileunit, '(I8,A)', advance='no') self%centerParticles(j,i), ','
        enddo
        write(fileunit,*) self%centerParticles(nc,i), ';'
    enddo
    do j=1, nc-1
        write(fileunit,'(I8,A)',advance='no') self%centerParticles(j,self%N),','
    enddo
    write(fileunit, *) self%centerParticles(nc,self%N), '];'

    write(fileunit,'(A)', advance='no') 'face_children = ['
    do i=1, self%N-1
        do j=1,3
            write(fileunit,'(I8,A)',advance='no') self%children(j,i), ','
        enddo
        write(fileunit,*) self%children(4,i), ';'
    enddo
    do j=1,3
        write(fileunit,'(I8,A)',advance='no') self%children(j,self%N), ','
    enddo
    write(fileunit,*) self%children(4,self%N), '];'

    write(fileunit, '(A)',advance='no') 'face_parent = ['
    do i=1, self%N-1
        write(fileunit, '(I8,A)',advance='no') self%parent(i), ','
    enddo
    write(fileunit, *) self%parent(self%N), '];'

    write(fileunit,'(A)',advance='no') 'face_area = ['
    do i=1, self%N-1
        write(fileunit,'(F24.12,A)',advance='no') self%area(i), ','
    enddo
    write(fileunit,*) self%area(self%N), '];'

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
	write(logKey,'(A,A,I0.3,A)') trim(logKey),'_',rank,' : '
	call New(aLog,logLevel)
	logInit = .TRUE.
end subroutine

end module
