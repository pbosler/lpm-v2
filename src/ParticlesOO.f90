module ParticlesOOModule
!------------------------------------------------------------------------------
! Lagrangian Particle Method (LPM) version 2.1
!------------------------------------------------------------------------------
!> @file
!> Provides the primitive Particles data structure that defines the spatial discretization of LPM.
!> @author
!> Peter Bosler, Sandia National Laboratories, Center for Computing Research
!
!> @defgroup Particles Particles
!> @brief Provides a vectorized data structure that defines the particles for LPM spatial discretization.
!>
!> Particles are combined with the @ref Field data structure to define scalar and vector fields over a spatial domain.
!> Particles may be combined with a mesh object (e.g., PolyMesh2d) or an unstructured data object (e.g., a quadtree)
!> to facilitate interpolation, differentiation, quadrature, etc.
!>
!> The particles data structures is a "structure of arrays," so that all information about a particle i is located at
!> index i in the relevant arrays.
!> For example, the x-coordinate in physical space of particle i is `aParticles%%x(i)`.
!> Its Lagrangian coordinates are given by `aParticles%%x0(i)`, `aParticles%%y0(i)`, etc.
!>
!> Particles that contribute to a midpoint quadrature rule are called "active," and this status is recorded in the `isActive` array.
!> Active particles have nonzero area (2d) or volume (3d).
!> Particles that do not contribute to the midpoint rule (i.e., vertex particles) are considered "passive."
!> Passive particles have zero area or volume.
!>
!> Each particle has physical coordinates @f$ (x(t),y(t),z(t)) @f$ and Lagrangian coordinates @f$ (x_0, y_0, z_0) @f$;
!> area (for 2d models) or volume (for 3d models).
!>
!> Particle sets used with a mesh are initialized in the @ref PolyMesh2d module.
!>
!>
!> @{
!
!------------------------------------------------------------------------------
use NumberKindsModule
use OutputWriterModule
use UtilitiesModule
use LoggerModule

implicit none
private
public Particles

!> @brief Class used to define a spatial discretization that may move in physical space.
type Particles
    real(kreal), allocatable :: x(:)   !< physical coordinate
	real(kreal), allocatable :: y(:)   !< physical coordinate
	real(kreal), allocatable :: z(:)   !< physical coordinate
	real(kreal), allocatable :: x0(:)  !< Lagrangian coordinate
	real(kreal), allocatable :: y0(:)  !< Lagrangian coordinate
	real(kreal), allocatable :: z0(:)  !< Lagrangian coordinate
	real(kreal), allocatable :: weight(:) !< weight (e.g., area, volume, mass) assigned to each particle
	integer(kint) :: N = 0 !< number of particles currently in use
	integer(kint) :: N_Max = 0 !< maximum number of particles allowed in memory
	integer(kint) :: geomkind = 0 !< ! geometry identifier e.g., @ref NumberKinds::planar_geom

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
	    procedure :: writeMatlab
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

!> @brief Initialize an empty particle set with space for up to nMax particles.
!> All values are set to zero until separately set.
!>
!> @param self Target particles object
!> @param nMax amount of space to reserve
!> @param geomKind Geometry kind (e.g., spherical or planar) as defined in @ref NumberKinds
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

!> @brief Deletes and frees memory associated with a particles object
!> @param self Target Particles object
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

!> @brief Performs a deep copy of one particles object into another
!> @param self Target Particles object
!> @param other Source Particles object
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

!> @brief Writes particles information to console using a loggermodule::logger object for formatting.
!> @param self
!> @param aLog
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

!> @brief Computes the sum of the particles%weight array
!> @Warning: only for the midpoint rule is this equivalent to a numerical integration.
pure function totalWeight(self)
    real(kreal) :: totalWeight
    class(Particles), intent(in) :: self
    totalWeight = sum(self%weight(1:self%N))
end function

!> @brief Redefines the Lagrangian coordinate to be equal to the physical coordinate
subroutine resetLagrangianCoordinates(self)
    class(Particles), intent(inout) :: self
    self%x0(1:self%N) = self%x(1:self%N)
    self%y0(1:self%N) = self%y(1:self%N)
    if (allocated(self%z)) self%z0(1:self%N) = self%z(1:self%N)
end subroutine

!> @brief Returns a particle's physical coordinate vector
!> @param self
!> @param index
!> @return PhysCoord coordinate vector
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

!> @brief Returns a particle's Lagrangian coordinate vector
!> @param self
!> @param index
!> @return LagCoord coordinate vector
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

!> @brief Inserts a single particle into a particles object, at the end of the arrays.
!> @param[inout] self a set of particles
!> @param[in] physX physical coordinate vector of particle to add
!> @param[in] lagX Lagrangian coordinate vector of particle to add
subroutine insert(self, physX, lagX)
    class(Particles), intent(inout) :: self
    real(kreal), dimension(:), intent(in) :: physX, lagX

    if (self%N+1 > self%N_Max) then
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

!> @brief Re-locates a single particle.
!> @param[inout] self a set of particles
!> @param[in] particle to be overwritten.
!> @param[in] physX physical coordinate vector of particle to add
!> @param[in] lagX Lagrangian coordinate vector of particle to add
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

!> @brief Writes particles to a script .m file readable by Matlab
!> @param self
!> @param fileunit
!> @todo WriteToCSV for python
subroutine writeMatlab(self, fileunit)
    class(Particles), intent(in) :: self
    integer(kint), intent(in) :: fileunit
    !
    integer(kint) :: i

    if (allocated(self%z)) then
        write(fileunit,'(A)',advance='no') 'xyz = ['
        do i=1,self%N-1
            write(fileunit,*) self%x(i), ',', self%y(i),',',self%z(i),';'
        enddo
        write(fileunit,*) self%x(self%N), ',', self%y(self%N),',',self%z(self%N),'];'
        write(fileunit,'(A)', advance='no') 'xyz0 = ['
        do i=1,self%N-1
            write(fileunit,*) self%x0(i), ',', self%y0(i),',',self%z0(i),';'
        enddo
        write(fileunit,*) self%x0(self%N), ',', self%y0(self%N),',',self%z0(self%N),'];'
    else
        write(fileunit,'(A)',advance='no') 'xy = ['
        do i=1,self%N-1
            write(fileunit,*) self%x(i), ',', self%y(i),';'
        enddo
        write(fileunit,*) self%x(self%N), ',', self%y(self%N),'];'
        write(fileunit,'(A)', advance='no') 'xyz0 = ['
        do i=1,self%N-1
            write(fileunit,*) self%x0(i), ',', self%y0(i),';'
        enddo
        write(fileunit,*) self%x0(self%N), ',', self%y0(self%N),'];'
    endif
    write(fileunit,'(A)',advance='no') 'weight = ['
    do i=1,self%N-1
        write(fileunit,*) self%weight(i), ','
    enddo
    write(fileunit,*) self%weight(self%N), '];'
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

!>@}
end module
