module ParticlesOOModule

use NumberKindsModule
use OutputWriterModule
use LoggerModule

implicit none
private
public Particles

type Particles 
    real(kreal), allocatable :: x(:)   !< physical coordinate
	real(kreal), allocatable :: y(:)   !< physical coordinate
	real(kreal), allocatable :: z(:)   !< physical coordinate
	real(kreal), allocatable :: x0(:)  !< Lagrangian coordinate
	real(kreal), allocatable :: y0(:)  !< Lagrangian coordinate
	real(kreal), allocatable :: z0(:)  !< Lagrangian coordinate
	real(kreal), allocatable :: weight(:)
	integer(kint) :: N = 0
	integer(kint) :: N_Max = 0
	integer(kint) :: geomkind = 0
	
	contains
	    procedure :: init
	    final :: deletePrivate
	    procedure :: copy
	    procedure :: logStats
!	    procedure :: printDebugInfo
	    procedure :: totalWeight
	    procedure :: resetLagrangianCoordinates
	    procedure :: physCoord
	    procedure :: lagCoord
	    procedure :: insert
	    procedure :: replace
end type
!
!----------------
! Logging
!----------------
!
logical(klog), save :: logInit = .FALSE.
type(Logger) :: log
character(len=28), save :: logKey = 'ParticlesLog'
integer(kint), parameter :: logLevel = DEBUG_LOGGING_LEVEL
character(len=MAX_STRING_LENGTH) :: logString

contains

subroutine init(self, nMax, geomKind)
    class(Particles), intent(inout) :: self
    integer(kint), intent(in) :: nMax
    integer(kint), intent(in) :: geomkind
    
    if (.NOT. logInit ) call InitLogger(log, procRank)
    
    ! error checking
	if ( nMax <= 0 ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL, logKey, " invalid nMax.")
		return
	endif
	if ( geomKind /= SPHERE_GEOM .AND. geomKind /= PLANAR_GEOM ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL, logKey, " invalid geometry.")
		return
	endif
	
	self%N_Max = nMax
	self%N = 0
	self%geomKind = geomKind
	allocate(self%x(nMax))
	self%x = 0.0_kreal
	allocate(self%y(nMax))
	self%y = 0.0_kreal
	allocate(self%x0(nMax))
	self%x0 = 0.0_kreal
	allocate(self%y0(nMax))
	self%y0 = 0.0_kreal
	allocate(self%weight(nMax))
	self%weight = 0.0_kreal
	if ( geomKind /= PLANAR_GEOM ) then
		allocate(self%z(nMax))
		self%z = 0.0
		allocate(self%z0(nMax))
		self%z0 = 0.0
	endif	
end subroutine

subroutine deletePrivate(self)
    type(Particles), intent(inout) :: self
    if ( allocated(self%x)) deallocate(self%x)
	if ( allocated(self%y)) deallocate(self%y)
	if ( allocated(self%x0)) deallocate(self%x0)
	if ( allocated(self%y0)) deallocate(self%y0)
	if ( allocated(self%weight)) deallocate(self%weight)
	if ( allocated(self%z)) deallocate(self%z)
	if ( allocated(self%z0)) deallocate(self%z0)
end subroutine

subroutine copy(self, other)
    class(Particles), intent(inout) :: self
    class(Particles), intent(in) :: other
    
    if ( self%N_Max < other%N ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL,logkey//" CopyParticles ERROR : ", " not enough memory.")
		return
	endif
	if ( allocated(other%z) .AND. .NOT. allocated(self%z) ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL,logkey//" CopyParticles ERROR : ", " dimension mismatch.")
		return
	endif
	
	self%x(1:other%N) = other%x(1:other%N)
	self%x0(1:other%N) = other%x0(1:other%N)
	self%y(1:other%N) = other%y(1:other%N)
	self%y0(1:other%N) = other%y0(1:other%N)
	self%weight(1:other%N) = other%weight(1:other%N)
	if (allocated(self%z)) then
	    self%z(1:other%N) = other%z(1:other%n)
	    self%z0(1:other%N) = other%z(1:other%n)
	endif
end subroutine

subroutine logStats(self, alog)
    class(Particles), intent(in) :: self
    type(Logger), intent(inout) :: alog
    
    call LogMessage(aLog, TRACE_LOGGING_LEVEL, logKey, " Particles Stats:")
    call StartSection(aLog)
    call LogMessage(aLog, TRACE_LOGGING_LEVEL, "particles.N = ", self%N )
	call LogMessage(aLog, TRACE_LOGGING_LEVEL, "particles.N_Max = ", self%N_Max )
    call LogMessage(aLog, TRACE_LOGGING_LEVEL, "n particles with nonzero weight = ", count(self%weight > dzero))
	call LogMessage(aLog, TRACE_LOGGING_LEVEL, "max x = ", maxval(self%x) )
	call LogMessage(aLog, TRACE_LOGGING_LEVEL, "min x = ", minval(self%x) )
	call LogMessage(aLog, TRACE_LOGGING_LEVEL, "max y = ", maxval(self%y) )
	call LogMessage(aLog, TRACE_LOGGING_LEVEL, "min y = ", minval(self%y) )
	call LogMessage(alog, TRACE_LOGGING_LEVEL, "min weight = ", minval(self%weight))
	call LogMessage(aLog, TRACE_LOGGING_LEVEL, "max weight = ", maxval(self%weight))
	call LogMessage(aLog, TRACE_LOGGING_LEVEL, "total weight = ", sum(self%weight))
    call EndSection(aLog)
end subroutine

pure function totalWeight(self)
    real(kreal) :: totalWeight
    class(Particles), intent(in) :: self
    totalWeight = sum(self%weight(1:self%N))
end function

subroutine resetLagrangianCoordinates(self)
    class(Particles), intent(inout) :: self
    self%x0(1:self%N) = self%x(1:self%N)
    self%y0(1:self%N) = self%y(1:self%N)
    if (allocated(self%z)) self%z0(1:self%N) = self%z(1:self%N)
end subroutine

pure function physCoord(self, index)
    real(kreal), dimension(3) :: physCoord
    class(Particles), intent(in) :: self
    integer(kint), intent(in) :: index
    physCoord(1) = self%x(index)
    physCoord(2) = self%y(index)
    if (allocated(self%z)) then
        physCoord(3) = self%z(index)
    else
        physCoord(3) = dzero
    endif
end function

pure function lagCoord(self, index)
    real(kreal), dimension(3) :: lagCoord
    class(Particles), intent(in) :: self
    integer(kint), intent(in) :: index
    lagCoord(1) = self%x(index)
    lagCoord(2) = self%y(index)
    if (allocated(self%z)) then
        lagCoord(3) = self%z(index)
    else
        lagCoord(3) = dzero
    endif
end function

subroutine insert(self, physX, lagX)
    class(Particles), intent(inout) :: self
    real(kreal), dimension(:), intent(in) :: physX, lagX
    
    if (self%N+1 >= self%N_Max) then
        call LogMessage(log, ERROR_LOGGING_LEVEL, logkey, " insert error : not enough memory.")
        return
    endif
    
    self%x(self%N+1) = physX(1)
    self%y(self%N+1) = physX(2)
    self%x0(self%N+1) = lagX(1)
    selF%y0(self%N+1) = lagX(2)
    if (self%geomkind /= PLANAR_GEOM) then
        self%z(self%N+1) = physX(3)
        self%z0(self%N+1) = lagX(3)
    endif
    self%N = self%N + 1
end subroutine

subroutine replace(self, index, physx, lagx)
    class(Particles), intent(inout) :: self
    integer(kint), intent(in) :: index
    real(kreal), dimension(3), intent(in) :: physx, lagx
    
    self%x(index) = physx(1)
    self%y(index) = physx(2)
    self%x0(index) = lagx(1)
    self%y0(index) = lagx(2)
    if (self%geomkind /= PLANAR_GEOM) then
        self%z(index) = physx(3)
        self%z0(index) = lagx(3)
    endif
end subroutine

!> @brief Initializes a logger for the Particles module
!> 
!> Output is controlled both by message priority and by MPI Rank
!> @param aLog Target Logger object
!> @param rank Rank of this processor
subroutine InitLogger(aLog,rank)
	type(Logger), intent(out) :: aLog
	integer(kint), intent(in) :: rank
	write(logKey,'(A,A,I0.3,A)') trim(logKey),'_',rank,' : '
	if ( rank == 0 ) then
		call New(aLog,logLevel)
	else
		call New(aLog,ERROR_LOGGING_LEVEL)
	endif
	logInit = .TRUE.
end subroutine

end module