module CubicEdgesModule

use NumberKindsModule
use LoggerModule
use EdgesModule
use ParticlesModule
use SphereGeomModule
use PlaneGeomModule

implicit none
private

public CubicEdges
public New, Delete, Copy
public DivideEdge, InsertEdge

type, extends(Edges) :: CubicEdges
    integer(kint), allocatable :: interiorParticles(:,:)
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

interface InsertEdge
    module procedure insertPrivate
end interface

interface DivideEdge
    module procedure divideCubicEdge
end interface
!
!----------------
! Logging
!----------------
!
logical(klog), save :: logInit = .FALSE.
type(Logger) :: log
character(len=28), save :: logKey = 'CubicEdges'
integer(kint), parameter :: logLevel = DEBUG_LOGGING_LEVEL

contains

subroutine newPrivate(self, nMax)
    type(CubicEdges), intent(out) :: self
    integer(kint), intent(in) :: nMax
    
    if (.NOT. logInit) call InitLogger(log, procRank)
    
    if (nMax <= 0) then
        call LogMessage(log, ERROR_LOGGING_LEVEL, logkey, " New CubicEdges error: invalid nMax")
        return
    endif
    
    call New(self%edges, nMax)
    
    allocate(self%interiorParticles(2, nMax))
    self%interiorParticles = 0
end subroutine

subroutine deletePrivate(self)
    type(CubicEdges), intent(inout) :: self
    
    if (allocated(self%interiorParticles)) deallocate(self%interiorParticles)
    call Delete(self%edges)
end subroutine

subroutine copyPrivate(self, other)
    type(CubicEdges), intent (inout) :: self
    type(CubicEdges), intent(in) :: other
    !
    integer(kint) :: j
    
    call Copy(self%edges, other%edges)
    
    do j=1, other%N
        self%interiorParticles(:,j) = other%interiorParticles(:,j)
    enddo
    
end subroutine

subroutine divideCubicEdge(self, edgeIndex, aParticles)
    type(CubicEdges), intent(inout) :: self
    integer(kint), intent(in) :: edgeIndex
    type(Particles), intent(inout) :: aParticles
    !
    real(kreal), dimension(3) :: midPt, lagMidPt, origXyz, destXyz, lagOrig, lagDest
    real(kreal), dimension(3,4) :: newCoords, newLagCoords
    integer(kint) :: j, particleInsertIndex
    integer(kint), dimension(4) :: parentParticles
    
    if ( self%N + 2 > self%N_Max ) then
        call LogMessage(log, ERROR_LOGGING_LEVEL, logkey, " divideCubicEdge : out of memory.")
        return
    endif
    
    origXyz = PhysCoord(aParticles, self%orig(edgeIndex))
    lagOrig = LagCoord(aParticles, self%orig(edgeIndex))
    destXyz = PhysCoord(aParticles, self%dest(edgeIndex))
    lagDest = LagCoord(aParticles, self%dest(edgeIndex))
    
    parentParticles(1) = self%orig(edgeIndex)
    parentParticles(2:3) = self%interiorParticles(:, edgeIndex)
    parentParticles(4) = self%dest(edgeIndex)
    
    if (aParticles%geomKind == SPHERE_GEOM) then
        midPt = SphereMidpoint(origXyz, destXyz)
        lagMidPt = SphereMidpoint(lagOrig, lagDest)
        
        newCoords(:,1) = pointAlongSphereVector(origXyz, midPt, -oosqrt5)
        newCoords(:,2) = pointAlongSphereVector(origXyz, midPt, oosqrt5)
        newCoords(:,3) = pointAlongSphereVector(midPt, destXyz, -oosqrt5)
        newCoords(:,4) = pointAlongSphereVector(midPt, destXyz, oosqrt5)
        
        newLagCoords(:,1) = pointAlongSphereVector(lagOrig, lagMidpt, -oosqrt5)
        newLagCoords(:,2) = pointAlongSphereVector(lagOrig, lagMidpt, oosqrt5)
        newLagCoords(:,3) = pointAlongSphereVector(lagMidpt, lagDest, -oosqrt5)
        newLagCoords(:,4) = pointAlongSphereVector(lagMidpt, lagDest, oosqrt5)
    else
        midPt = 0.5_kreal * (origXyz + destXyz)
        lagMidPt = 0.5_kreal * (lagOrig + lagDest)
        
        newCoords(:,1) = pointAlongChordVector(origXyz, midPt, -oosqrt5)
        newCoords(:,2) = pointAlongChordVector(origXyz, midPt, oosqrt5)
        newCoords(:,3) = pointAlongChordVector(midPt, destXyz, -oosqrt5)
        newCoords(:,4) = pointAlongChordVector(midPt, destXyz, oosqrt5)
        
        newLagCoords(:,1) = pointAlongChordVector(lagOrig, lagMidpt, -oosqrt5)
        newLagCoords(:,2) = pointAlongChordVector(lagOrig, lagMidpt, oosqrt5)
        newLagCoords(:,3) = pointAlongChordVector(lagMidpt, lagDest, -oosqrt5)
        newLagCoords(:,4) = pointAlongChordVector(lagMidpt, lagDest, oosqrt5)
    endif
    
    particleInsertIndex = aParticles%N + 1
    call ReplaceParticle(aParticles, parentParticles(2), newCoords(:,1), newLagCoords(:,1))
    call InsertParticle(aParticles, newCoords(:,2), newLagCoords(:,2)) ! particeInsertIndex
    call InsertParticle(aParticles, midPt, lagMidPt) ! particeInsertIndex +1
    call InsertParticle(aParticles, newCoords(:,3), newLagCoords(:,3)) ! particeInsertIndex + 2
    call ReplaceParticle(aParticles, parentParticles(3), newCoords(:,4), newLagCoords(:,4))
    
    self%hasChildren(edgeIndex) = .TRUE.
    self%child1(edgeIndex) = self%N + 1
    self%child2(edgeIndex) = self%N + 2
    self%parent(self%N + 1) = edgeIndex
    self%parent(self%N + 2) = edgeIndex
    
    self%orig(self%N + 1) = parentParticles(1)
    self%dest(self%N + 1) = particleInsertIndex + 1
    self%leftFace(self%N + 1) = self%leftFace(edgeIndex)
    self%rightFace(self%N + 1) = self%rightFace(edgeIndex)
    self%interiorParticles(:, self%N+1) = (/ parentParticles(2), particleInsertIndex /)
    
    self%orig(self%N +2) = particleInsertIndex + 1
    self%dest(self%N +2) = parentParticles(4)
    self%leftFace(self%N+2) = self%leftFace(edgeIndex)
    self%rightFace(self%N+2) = self%rightFace(edgeIndex)
    self%interiorParticles(:, self%N+2) = (/ particleInsertIndex +2, parentParticles(3) /)
    
    self%N = self%N + 2
end subroutine

subroutine insertPrivate(self, origIndex, destIndex, leftFace, rightFace, interiorIndices)
    type(CubicEdges), intent(inout) :: self
    integer(kint), intent(in) :: origIndex, destIndex, leftFace, rightFace
    integer(kint), dimension(2), intent(in) :: interiorIndices
    !
    integer(kint) :: nn
    
    if ( self%N >= self%N_Max ) then
        call LogMessage(log, ERROR_LOGGING_LEVEL, logKey, " InsertEdge : out of memory. ")
        return
    endif
    
    nn = self%N
    
    self%orig(nn+1) = origIndex
    self%dest(nn+1) = destIndex
    self%leftFace(nn+1) = leftFace
    self%rightFace(nn+1) = rightFace
    self%interiorParticles(:,nn+1) = interiorIndices
    
    self%N = nn+1
end subroutine

end module