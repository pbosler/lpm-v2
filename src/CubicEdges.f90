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
character(len=28), save :: logKey = 'Edges'
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
end subroutine


end module