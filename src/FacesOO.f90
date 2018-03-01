module FacesOOModule

use NumberKindsModule
use STDIntVectorModule
use LoggerModule
use ParticlesOOModule
use EdgesOOModule
use PlaneGeomModule
use SphereGeomModule

implicit none
private

public Faces, bilinearMap, TriLinearFaces, QuadLinearFaces, QuadCubicFaces

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
	
	contains
	    procedure :: init
	    procedure :: copy
        procedure :: insert
        procedure(divide), deferred :: divide
        procedure :: physCentroid
        procedure :: lagCentroid
        procedure :: countParents
        procedure :: area
        procedure :: sharedEdge
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


type, extends(Faces) :: TriLinearFaces
    contains
        procedure :: divide => divideTri
        final :: deleteTri
        procedure :: particleOppositeTriEdge
        procedure, private :: triFaceArea
end type

type, extends(Faces) :: QuadLinearFaces
    contains    
        procedure :: divide => divideQuadLinear
        final :: deleteQuad
        procedure, private :: quadFaceArea
end type

type, extends(QuadLinearFaces) :: QuadCubicFaces
    contains
        procedure :: divide => divideQuadCubic
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
    self%children = 0
    self%hasChildren = .FALSE.
    self%parent = 0
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
end subroutine

subroutine deleteQuad(self)
    type(QuadLinearFaces), intent(inout) :: self
    if (allocated(self%vertices)) deallocate(self%vertices)
    if (allocated(self%edges)) deallocate(self%edges)
    if (allocated(self%centerParticles)) deallocate(self%centerParticles)
    if (allocated(self%children)) deallocate(self%children)
    if (allocated(self%hasChildren)) deallocate(self%hasChildren)
    if (allocated(self%parent)) deallocate(self%parent)
end subroutine

subroutine deleteQuadCubic(self)
    type(QuadCubicFaces), intent(inout) :: self
    if (allocated(self%vertices)) deallocate(self%vertices)
    if (allocated(self%edges)) deallocate(self%edges)
    if (allocated(self%centerParticles)) deallocate(self%centerParticles)
    if (allocated(self%children)) deallocate(self%children)
    if (allocated(self%hasChildren)) deallocate(self%hasChildren)
    if (allocated(self%parent)) deallocate(self%parent)
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
end subroutine

subroutine insert(self, centerInds, vertInds, edgeInds)
    class(Faces), intent(inout) :: self
    integer(kint), dimension(:), intent(in) :: centerInds, vertInds, edgeInds
    
    if (self%N+1 >= self%N_Max) then
        call LogMessage(log, ERROR_LOGGING_LEVEL, logkey, " insert error : not enough memory.")
        return
    endif
    if (size(centerInds) /= size(self%centerParticles,1) .or. size(vertInds) /= size(self%vertices,1) .or. &
        size(edgeInds) /= size(vertInds) ) then
        call LogMessage(log, ERROR_LOGGING_LEVEL, logkey, " insert error : size mismatch.")
        return
    endif
    
    self%centerParticles(:,self%N+1) = centerInds
    self%vertices(:,self%N+1) = vertInds
    self%edges(:,self%N+1) = edgeInds
    self%N = self%N + 1
end subroutine

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
	    aParticles%weight(self%centerParticles(1,self%N+i)) = self%quadFaceArea(self%N+i, aParticles)
	enddo
	aParticles%weight(self%centerParticles(1,index)) = dzero
	self%hasChildren(index) = .TRUE.
	self%N = self%N+4
	self%N_Active = self%N_Active + 3
end subroutine

pure function bilinearMap(vertXyz, s1, s2)
    real(kreal), dimension(3) :: bilinearMap
    real(kreal), dimension(3,4), intent(in) :: vertXyz
    real(kreal), intent(in) :: s1, s2
    
    bilinearMap = 0.25_kreal * ( (1.0_kreal-s1)*(1.0_kreal-s2)*vertXyz(:,1) + (1.0_kreal+s1)*(1.0_kreal-s2)*vertXyz(:,2) + &
        (1.0_kreal+s1)*(1.0_kreal+s2)*vertXyz(:,3) + (1.0_kreal-s1)*(1.0_kreal+s2)*vertXyz(:,4))
end function

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
    integer(kint) :: nParticles, nEdges
    integer(kint), dimension(12,4) :: newFaceVerts
    integer(kint), dimension(4,4) :: newFaceEdges, newFaceCenters
    integer(kint), dimension(4) :: parentEdgeMidpoint
    real(kreal), dimension(3) :: quadCtr, lagQuadCtr
    real(kreal), dimension(3,8) :: edgePts, lagEdgePts
    integer(kint), dimension(8) :: edgePtInds
    real(kreal), dimension(3,4) :: quadCoords, lagQuadCoords
    real(kreal), dimension(3,4) :: newPhysCenters, newLagCenters
    
    if ( self%N_Max < self%N + 4 ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL, logkey, " DivideQuadFace ERROR : not enough memory.")
		return
	endif
	select type(anEdges)
	class is (CubicEdges)
        newFaceVerts = 0
        do i=1,12
            newFaceVerts(i, mod(i/3,4) + 1) = self%vertices(i,index)
        enddo
        do i=1,4
            newFaceCenters(i,i) = self%centerParticles(i,index)
        enddo
    
        !
        !   loop over parent edges to divide face boundaries
        !
        do i=1, 4
            parentEdge = self%edges(i,index)
            if (anEdges%hasChildren(parentEdge)) then
                !
                !   edge has already been divided by adjacent panel
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
                anEdges%leftFace(childEdge1) = self%N+i
            
                newFaceEdges(i, mod(i,4)+1) = childEdge2
                anEdges%leftFace(childEdge2) = self%N + mod(i,4) + 1
            else
                newFaceEdges(i,i) = childEdge2
                anEdges%rightFace(childEdge2) = self%N+i
            
                newFaceEdges(i, mod(i,4)+1) = childEdge1
                anEdges%rightFace(childEdge1) = self%N + mod(i,4) + 1
            endif
            
            parentEdgeMidpoint(i) = anEdges%dest(childEdge1)
        enddo
        !
        !   connect boundary edges and vertices to child faces
        !
        do i=1,4
            parentEdge = self%edges(i,index)
            childEdge1 = anEdges%child1(parentEdge)
            childEdge2 = anEdges%child2(parentEdge)
            if (anEdges%dest(childEdge1) /= anEdges%orig(childEdge2)) then
                call LogMessage(log, ERROR_LOGGING_LEVEL, &
                    trim(logkey)//" connectivity error, dividequadpanel, edge children, subpanel ", i)
            endif
            if (anEdges%positiveOrientation(parentEdge, index)) then
                select case (i)
                case(1)
                    if (anEdges%orig(childEdge1) /= self%vertices(1,index)) then
                        call LogMessage(log, ERROR_LOGGING_LEVEL, &
                            trim(logkey)//" connectivity error, dividequadpanel, subpanel 1 ", "edge 1, vertex 1")
                    endif
                    if (anEdges%interiorParticles(1,childEdge1) /= self%vertices(2,index)) then
                        call LogMessage(log, ERROR_LOGGING_LEVEL, &
                            trim(logkey)//" connectivity error, dividequadpanel, subpanel 1 ", "edge 1, int. p1")
                    endif
            
                    newFaceVerts(1,1) = anEdges%orig(childEdge1)
                    newFaceVerts(2:3,1) = anEdges%interiorParticles(:,childEdge1)
                    newFaceVerts(4,1) = anEdges%dest(childEdge1)
                    
                    newFaceVerts(1,2) = anEdges%orig(childEdge2)
                    newFaceVerts(2:3,2) = anEdges%interiorParticles(:, childEdge2)
                    newFaceVerts(4,2) = anEdges%dest(childEdge2)
                    
                    if (newFaceVerts(3,2) /= anEdges%interiorParticles(2,childEdge2)) then
                        call LogMessage(log, ERROR_LOGGING_LEVEL, &
                            trim(logkey)//" connectivity error, dividequadpanel, subpanel 2 ", "vertex 3, child2 int. p2")
                    endif
                case(2)
                    if (newFaceVerts(4,2) /= anEdges%orig(childEdge1)) then
                        call LogMessage(log, ERROR_LOGGING_LEVEL, &
                            trim(logkey)//" connectivity error, dividequadpanel, subpanels ", "1 and 2")
                    endif
                    if (anEdges%interiorParticles(1,childEdge1) /= self%vertices(5,index)) then
                        call LogMessage(log, ERROR_LOGGING_LEVEL, &
                            trim(logkey)//" connectivity error, dividequadpanel, subpanel 2 ", "vertex 5, child1 int. p1")
                    endif
                    newFaceVerts(5:6,2) = anEdges%interiorParticles(:,childEdge1)
                    newFaceVerts(7,2) = anEdges%dest(childEdge1)
                    
                    newFaceVerts(4,3) = anEdges%orig(childEdge2)
                    newFaceVerts(5:6,3) = anEdges%interiorParticles(:,childEdge2)
                    newFaceVerts(7,3) = anEdges%dest(childEdge2)
                case(3)
                    if (newFaceVerts(7,3) /= anEdges%orig(childEdge1)) then
                        call LogMessage(log, ERROR_LOGGING_LEVEL, &
                            trim(logkey)//" connectivity error, dividequadpanel, subpanels ", "2 and 3")
                    endif
                    if (anEdges%interiorParticles(1,childEdge1) /= self%vertices(8,index)) then
                        call LogMessage(log, ERROR_LOGGING_LEVEL, &
                            trim(logkey)//" connectivity error, dividequadpanel, subpanel 3 ", "vertex 8, child1 int. p1")
                    endif
                    newFaceVerts(8:9,3) = anEdges%interiorParticles(:, childEdge1)
                    newFaceVerts(10,3) = anEdges%dest(childEdge1)
                    
                    newFaceVerts(7,4) = anEdges%orig(childEdge2)
                    newFaceVerts(8:9,4) = anEdges%interiorParticles(:,childedge2)
                    newFaceVerts(10,4) = anEdges%dest(childEdge2)
                case(4)
                    if (newFaceVerts(10,4) /= anEdges%orig(childEdge1)) then
                        call LogMessage(log, ERROR_LOGGING_LEVEL, &
                            trim(logkey)//" connectivity error, dividequadpanel, subpanels ", "3 and 4")
                    endif
                    if (anEdges%interiorParticles(1,childEdge1) /= self%vertices(11,index)) then
                        call LogMessage(log, ERROR_LOGGING_LEVEL, &
                            trim(logkey)//" connectivity error, dividequadpanel, subpanel 4 ", "vertex 11, child1 int. p1")
                    endif
                    newFaceVerts(11:12,4) = anEdges%interiorParticles(:,childEdge1)
                    
                    newFaceVerts(11:12,1) = anEdges%interiorParticles(:,childEdge2)
                    if (anEdges%dest(childEdge2) /= newFaceVerts(1,1)) then
                        call LogMessage(log, ERROR_LOGGING_LEVEL, &
                            trim(logkey)//" connectivity error, dividequadpanel, subpanel 4 ", "vertex 12, child2 dest")
                    endif
                end select
            else
                select case (i)
                case(1)
                    if (anEdges%dest(childEdge2) /= self%vertices(1,index)) then
                        call LogMessage(log, ERROR_LOGGING_LEVEL, &
                            trim(logkey)//" connectivity error, dividequadpanel, subpanel 1 ", "edge 1, vertex 1")
                    endif
                    if (anEdges%interiorParticles(2,childEdge2) /= self%vertices(2,index)) then
                        call LogMessage(log, ERROR_LOGGING_LEVEL, &
                            trim(logkey)//" connectivity error, dividequadpanel, subpanel 1 ", "edge 1, int. p2")
                    endif
                    
                    newFaceVerts(1,1) = anEdges%dest(childEdge2)
                    newFaceVerts(2,1) = anEdges%interiorParticles(2,childEdge2)
                    newFaceVerts(3,1) = anEdges%interiorParticles(1,childEdge2)
                    newFaceVerts(4,1) = anEdges%orig(childEdge2)
                    
                    newFaceVerts(1,2) = anEdges%dest(childEdge1)
                    newFaceVerts(2,2) = anEdges%interiorParticles(2, childEdge1)
                    newFaceVerts(3,2) = anEdges%interiorParticles(1, childEdge1)
                    newFaceVerts(4,2) = anEdges%orig(childEdge1)
                    
                    if (newFaceVerts(3,2) /= anEdges%interiorParticles(1,childEdge1)) then
                        call LogMessage(log, ERROR_LOGGING_LEVEL, &
                            trim(logkey)//" connectivity error, dividequadpanel, subpanel 2 ", "vertex 3, child1 int. p1")
                    endif
                case(2)
                    if (newFaceVerts(4,2) /= anEdges%dest(childEdge2)) then
                        call LogMessage(log, ERROR_LOGGING_LEVEL, &
                            trim(logkey)//" connectivity error, dividequadpanel, subpanels ", "1 and 2")
                    endif
                    if (anEdges%interiorParticles(2,childEdge2) /= self%vertices(5,index)) then
                        call LogMessage(log, ERROR_LOGGING_LEVEL, &
                            trim(logkey)//" connectivity error, dividequadpanel, subpanel 2 ", "vertex 5, child2 int. p2")
                    endif
                    newFaceVerts(5,2) = anEdges%interiorParticles(2,childEdge2)
                    newFaceVerts(6,2) = anEdges%interiorParticles(1,childEdge2)
                    newFaceVerts(7,2) = anEdges%orig(childEdge2)
                    
                    newFaceVerts(4,3) = anEdges%dest(childEdge1)
                    newFaceVerts(5,3) = anEdges%interiorParticles(2,childEdge1)
                    newFaceVerts(6,3) = anEdges%interiorParticles(1,childEdge1)
                    newFaceVerts(7,3) = anEdges%orig(childEdge1)
                case(3)
                    if (newFaceVerts(7,3) /= anEdges%Dest(childEdge2)) then
                        call LogMessage(log, ERROR_LOGGING_LEVEL, &
                            trim(logkey)//" connectivity error, dividequadpanel, subpanels ", "2 and 3")
                    endif
                    if (anEdges%interiorParticles(2,childEdge2) /= self%vertices(8,index)) then
                        call LogMessage(log, ERROR_LOGGING_LEVEL, &
                            trim(logkey)//" connectivity error, dividequadpanel, subpanel 3 ", "vertex 8, child2 int. p2")
                    endif
                    newFaceVerts(8,3) = anEdges%interiorParticles(2, childEdge2)
                    newFaceVerts(9,3) = anEdges%interiorParticles(1, childEdge2)
                    newFaceVerts(10,3) = anEdges%orig(childEdge2)
                    
                    newFaceVerts(7,4) = anEdges%dest(childEdge1)
                    newFaceVerts(8,4) = anEdges%interiorParticles(2,childedge1)
                    newFaceVerts(9,4) = anEdges%interiorParticles(1,childedge1)
                    newFaceVerts(10,4) = anEdges%orig(childEdge1)
                case(4)
                    if (newFaceVerts(10,4) /= anEdges%dest(childEdge2)) then
                        call LogMessage(log, ERROR_LOGGING_LEVEL, &
                            trim(logkey)//" connectivity error, dividequadpanel, subpanels ", "3 and 4")
                    endif
                    if (anEdges%interiorParticles(2,childEdge2) /= self%vertices(11,index)) then
                        call LogMessage(log, ERROR_LOGGING_LEVEL, &
                            trim(logkey)//" connectivity error, dividequadpanel, subpanel 4 ", "vertex 11, child2 int. p2")
                    endif
                    newFaceVerts(11,4) = anEdges%interiorParticles(2,childEdge2)
                    newFaceVerts(12,4) = anEdges%interiorParticles(1,childEdge2)
                    
                    newFaceVerts(11,1) = anEdges%interiorParticles(2,childEdge1)
                    newFaceVerts(12,1) = anEdges%interiorParticles(1,childEdge1)
                    if (anEdges%orig(childEdge1) /= newFaceVerts(1,1)) then
                        call LogMessage(log, ERROR_LOGGING_LEVEL, &
                            trim(logkey)//" connectivity error, dividequadpanel, subpanel 4 ", "vertex 12, child1 orig")
                    endif
                end select
            endif
        enddo
        
        !
        !   create new interior edges
        !
        do i=1,4
            quadCoords(:,i) = aparticles%physCoord(self%vertices(mod(3*i+9,12)+1,index))
            lagQuadCoords(:,i) = aParticles%lagCoord(self%vertices(mod(3*i+9,12)+1, index))
        enddo
        if(aParticles%geomKind == SPHERE_GEOM) then
            quadCtr = SphereCentroid(quadCoords)
            lagQuadCtr = SphereCentroid(lagQuadCoords)
            
            edgePts(:,1) = pointAlongSphereVector(aParticles%physCoord(parentEdgeMidpoint(1)), quadCtr, -oosqrt5)
            edgePts(:,2) = pointAlongSphereVector(aParticles%physCoord(parentEdgeMidpoint(1)), quadCtr, oosqrt5)
            edgePts(:,3) = pointAlongSphereVector(quadCtr, aParticles%physCoord(parentEdgeMidpoint(3)), -oosqrt5)
            edgePts(:,4) = pointAlongSphereVector(quadCtr, aParticles%physCoord(parentEdgeMidpoint(3)), oosqrt5)
            edgePts(:,5) = pointAlongSphereVector(aParticles%physCoord(parentEdgeMidpoint(2)), quadCtr, -oosqrt5)
            edgePts(:,6) = pointAlongSphereVector(aParticles%physCoord(parentEdgeMidpoint(2)), quadCtr, -oosqrt5)
            edgePts(:,7) = pointAlongSphereVector(quadCtr, aParticles%physCoord(parentEdgeMidpoint(4)), -oosqrt5)
            edgePts(:,8) = pointAlongSphereVector(quadCtr, aParticles%physCoord(parentEdgeMidpoint(4)), oosqrt5)
            lagEdgePts(:,1) = pointAlongSphereVector(aParticles%lagCoord(parentEdgeMidpoint(1)), lagQuadCtr, -oosqrt5)
            lagEdgePts(:,2) = pointAlongSphereVector(aParticles%lagCoord(parentEdgeMidpoint(1)), lagQuadCtr, oosqrt5)
            lagEdgePts(:,3) = pointAlongSphereVector(lagQuadCtr, aParticles%lagCoord(parentEdgeMidpoint(3)), -oosqrt5)
            lagEdgePts(:,4) = pointAlongSphereVector(lagQuadCtr, aParticles%lagCoord(parentEdgeMidpoint(3)), oosqrt5)
            lagEdgePts(:,5) = pointAlongSphereVector(aParticles%lagCoord(parentEdgeMidpoint(2)), lagQuadCtr, -oosqrt5)
            lagEdgePts(:,6) = pointAlongSphereVector(aParticles%lagCoord(parentEdgeMidpoint(2)), lagQuadCtr, -oosqrt5)
            lagEdgePts(:,7) = pointAlongSphereVector(lagQuadCtr, aParticles%lagCoord(parentEdgeMidpoint(4)), -oosqrt5)
            lagEdgePts(:,8) = pointAlongSphereVector(lagQuadCtr, aParticles%lagCoord(parentEdgeMidpoint(4)), oosqrt5)            
        else
            quadCtr = EuclideanCentroid(quadCoords)
            lagQuadCtr = EuclideanCentroid(lagQuadCoords)
            
            edgePts(:,1) = pointAlongChordVector(aParticles%physCoord(parentEdgeMidpoint(1)), quadCtr, -oosqrt5)
            edgePts(:,2) = pointAlongChordVector(aParticles%physCoord(parentEdgeMidpoint(1)), quadCtr, oosqrt5)
            edgePts(:,3) = pointAlongChordVector(quadCtr, aParticles%physCoord(parentEdgeMidpoint(3)), -oosqrt5)
            edgePts(:,4) = pointAlongChordVector(quadCtr, aParticles%physCoord(parentEdgeMidpoint(3)), oosqrt5)
            edgePts(:,5) = pointAlongChordVector(aParticles%physCoord(parentEdgeMidpoint(2)), quadCtr, -oosqrt5)
            edgePts(:,6) = pointAlongChordVector(aParticles%physCoord(parentEdgeMidpoint(2)), quadCtr, -oosqrt5)
            edgePts(:,7) = pointAlongChordVector(quadCtr, aParticles%physCoord(parentEdgeMidpoint(4)), -oosqrt5)
            edgePts(:,8) = pointAlongChordVector(quadCtr, aParticles%physCoord(parentEdgeMidpoint(4)), oosqrt5)
            lagEdgePts(:,1) = pointAlongChordVector(aParticles%lagCoord(parentEdgeMidpoint(1)), lagQuadCtr, -oosqrt5)
            lagEdgePts(:,2) = pointAlongChordVector(aParticles%lagCoord(parentEdgeMidpoint(1)), lagQuadCtr, oosqrt5)
            lagEdgePts(:,3) = pointAlongChordVector(lagQuadCtr, aParticles%lagCoord(parentEdgeMidpoint(3)), -oosqrt5)
            lagEdgePts(:,4) = pointAlongChordVector(lagQuadCtr, aParticles%lagCoord(parentEdgeMidpoint(3)), oosqrt5)
            lagEdgePts(:,5) = pointAlongChordVector(aParticles%lagCoord(parentEdgeMidpoint(2)), lagQuadCtr, -oosqrt5)
            lagEdgePts(:,6) = pointAlongChordVector(aParticles%lagCoord(parentEdgeMidpoint(2)), lagQuadCtr, -oosqrt5)
            lagEdgePts(:,7) = pointAlongChordVector(lagQuadCtr, aParticles%lagCoord(parentEdgeMidpoint(4)), -oosqrt5)
            lagEdgePts(:,8) = pointAlongChordVector(lagQuadCtr, aParticles%lagCoord(parentEdgeMidpoint(4)), oosqrt5) 
        endif
        nParticles = aParticles%N
        call aParticles%insert(quadCtr, lagQuadCtr)
        do i=1,8
            edgePtInds(i) = nParticles+i+1
            call aParticles%insert(edgePts(:,i), lagEdgePts(:,i))
        enddo
        nEdges = anEdges%N
        call anEdges%insert(parentEdgeMidpoint(1), nParticles+1, self%N+1, self%N+2, edgePtInds(1:2))
        call anEdges%insert(nParticles+1, parentEdgeMidpoint(3), self%N+4, self%N+3, edgePtInds(3:4))
        call anEdges%insert(parentEdgeMidpoint(2), nParticles+1, self%N+2, self%N+3, edgePtInds(5:6))
        call anEdges%insert(nParticles+1, parentEdgeMidpoint(4), self%N+1, self%N+4, edgePtInds(7:8))
        
        newFaceEdges(2,1) = nEdges+1
        newFaceEdges(4,2) = nEdges+1

        newFaceEdges(2,4) = nEdges+2
        newFaceEdges(4,3) = nEdges+2

        newFaceEdges(3,2) = nEdges+3
        newFaceEdges(1,3) = nEdges+3

        newFaceEdges(3,1) = nEdges+4
        newFaceEdges(1,4) = nEdges+4
        
        !
        !   center particles
        !
        nParticles = aParticles%N
        do i = 1, 4
            do j=1,4
                quadCoords(:,j) = aParticles%physCoord(newFaceVerts(mod(3*j+9,12)+1,i))
                lagQuadCoords(:,j) = aParticles%lagCoord(newFaceVerts(mod(3*j+9,12)+1,i))
            enddo
            newPhysCenters = self%calcInteriorPts(quadCoords)
            newLagCenters = self%calcInteriorPts(lagQuadCoords)
            !
            !   reposition parent particle
            !
            call aParticles%replace(self%centerParticles(i,index), newPhysCenters(:,i), newLagCenters(:,i))
            !
            !   insert new interior particles
            !
            do j=1,4
                if (i/=j) then
                    call aParticles%insert(newPhysCenters(:,j), newLagCenters(:,j))
                    self%centerParticles(j,self%N+i) = aParticles%N
                endif
            enddo
        enddo
        
        
       !
       !    create child faces
       !
       do i=1,4
           self%vertices(:,self%N+i) = newFaceVerts(:,i)
           self%edges(:,self%N+i) = newFaceEdges(:,i)
           self%children(i, index) = self%N+i
           self%parent(self%N+i) = index
       enddo
       self%N = self%N + 4
       self%N_Active = self%N_Active + 3
       self%hasChildren(index) = .TRUE.
              
    class default
        call logMessage(log, ERROR_LOGGING_LEVEL, trim(logkey)//" divideCubicFace error : ", " cubic edges required.")
        return 
    end select
end subroutine

pure function calcInteriorPts(self, vertXyz)
    class(QuadCubicFaces), intent(in) :: self
    real(kreal), dimension(3,4) :: calcInteriorPts
    real(kreal), dimension(3,4), intent(in) :: vertXyz
    calcInteriorPts(:,1) = bilinearMap(vertXyz, -oosqrt5, oosqrt5)
    calcInteriorPts(:,2) = bilinearMap(vertXyz, -oosqrt5, -oosqrt5)
    calcInteriorPts(:,3) = bilinearMap(vertXyz, oosqrt5, -oosqrt5)
    calcInteriorPts(:,4) = bilinearMap(vertXyz, oosqrt5, oosqrt5)
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
    do i=1, size(self%vertices,1)
        physCentroid = physCentroid + aParticles%physCoord(self%vertices(i, index))
    enddo
    physCentroid = physCentroid / real(size(self%vertices,1), kreal)
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
    integer(kint) :: i
    real(kreal) :: norm
    
    lagCentroid = dzero
    do i=1, size(self%vertices,1)
        lagCentroid = lagCentroid + aParticles%lagCoord(self%vertices(i, index))
    enddo
    lagCentroid = lagCentroid / real(size(self%vertices,1), kreal)
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
    if (self%N + 4 >= self%N_Max) then
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
	    aParticles%weight(self%centerParticles(1,self%N+i)) = self%triFaceArea(self%N+i, aParticles)
	enddo
	! special case for child 4: re-use parent face centerParticles
	self%centerParticles(1,self%N+4) = self%centerParticles(1,index)
	self%vertices(:,self%N+4) = newFaceVerts(:,4)
	self%edges(:,self%N+4) = newFaceEdges(:,4)
	self%children(4,index) = self%N+4
	self%parent(self%N+4) = index
	aParticles%weight(self%centerParticles(1,self%N+4)) = self%triFaceArea(self%N+4, aParticles)
	
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