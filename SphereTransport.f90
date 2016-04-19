module SphereTransportModule
!> @file SphereTransport.f90
!> Data structure for representing solutions of the Barotropic Vorticity Equation (BVE) on the surface of a rotating sphere
!> @author Peter Bosler, Sandia National Laboratories Center for Computing Research
!> 
!>
!> @defgroup SphereTransport SphereTransport
!> @brief Data structure for representing solutions of the advection equation on the surface of a rotating sphere.
!>
!> @{
use NumberKindsModule
use OutputWriterModule
use LoggerModule
use ParticlesModule
use EdgesModule
use FacesModule
use PolyMesh2dModule
use FieldModule
use MPISetupModule
use SphereGeomModule

implicit none

include 'mpif.h'

private
!
!----------------
! Module types and public declarations
!----------------
!
public TransportMesh, New, Delete, Copy
public AddTracers
public SetVelocityOnMesh, SetTracerOnMesh, SetInitialDensityOnMesh, SetDivergenceOnMesh
public LogStats
public OutputToVTK
public TracerMass

type TransportMesh
	type(PolyMesh2d) :: mesh
	type(Field) :: density
	type(Field) :: velocity
	type(Field) :: divergence
	type(Field), dimension(:), allocatable :: tracers
	real(kreal) :: radius = 1.0_kreal
	logical(klog) :: nonzeroDivergence = .TRUE.
	type(MPISetup) :: mpiParticles
	
	contains
		final :: deletePrivate
end type

!
!----------------
! Module interfaces
!----------------
!
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
	module procedure logStatsPrivate
end interface

interface SetVelocityOnMesh
	module procedure setVelocityFromFunction
end interface

interface SetTracerOnMesh
	module procedure setScalarTracerOnMesh
	module procedure setVectorTracerOnMesh
end interface


!
!----------------
! Logging
!----------------
!
logical(klog), save :: logInit = .FALSE.
type(Logger) :: log
character(len=28), save :: logKey = 'TRANS'
integer(kint), parameter :: logLevel = DEBUG_LOGGING_LEVEL

contains
!
!----------------
! private methods
!----------------
!
subroutine SetInitialDensityOnMesh( self, densFn )
	type(TransportMesh), intent(inout) :: self
	procedure(scalarFnOf3DSpace), optional :: densFn
	integer(kint) :: i
	
	if ( present(densFn) ) then
		do i = 1, self%mesh%particles%N
			call InsertScalarToField( self%density, densFn( self%mesh%particles%x0(i), self%mesh%particles%y0(i), &
				 self%mesh%particles%z0(i)) ) 
		enddo
	else
		do i = 1, self%mesh%particles%N
			call InsertScalarToField(self%density, 1.0_kreal )
		enddo
	endif
end subroutine

subroutine SetDivergenceOnMesh(self, divFn, t)
	type(TransportMesh), intent(inout) :: self
	procedure(scalarFnOf3DSpaceAndTime), optional :: divFn
	real(kreal), intent(in), optional :: t
	integer(kint) :: i
	
	self%divergence%N = self%mesh%particles%N
	if ( .NOT. present(divFn) ) then
		call SetFieldToZero(self%divergence)
	else
		do i = 1, self%mesh%particles%N
			self%divergence%scalar(i) = divFn(self%mesh%particles%x(i), self%mesh%particles%y(i), &
				 self%mesh%particles%z(i), t ) 
		enddo
	endif
end subroutine

subroutine AddTracers(self, nTracers, tracerDims )
	type(TransportMesh), intent(inout) :: self
	integer(kint), intent(in) :: nTracers
	integer(kint), dimension(nTracers), intent(in) :: tracerDims
	integer(kint) :: i
	
	allocate(self%tracers(nTracers))
	do i = 1, nTracers
		call New(self%tracers(i), tracerDims(i), self%mesh%particles%N_Max )
	enddo
end subroutine

subroutine OutputToVTK( self, filename )
	type(TransportMesh), intent(in) :: self
	character(len=*), intent(in) :: filename
	integer(kint) :: i, writeStat
	
	open( unit=WRITE_UNIT_1, file=filename, status='REPLACE', action='WRITE', iostat=writeStat)
		if ( writeStat /= 0 ) then
			call LogMessage(log,ERROR_LOGGING_LEVEL, trim(logKey)//" OutputToVTK ERROR opening file = ", trim(filename))
			return
		endif
	
		call WriteVTKPoints(self%mesh%particles, WRITE_UNIT_1)
		call WriteFacesToVTKPolygons( self%mesh%faces, WRITE_UNIT_1)
		
		call WriteVTKPointDataSectionHeader(WRITE_UNIT_1, self%mesh%particles%N)
		
		call WriteVTKLagCoords( self%mesh%particles, WRITE_UNIT_1)
		
		call WriteFieldToVTKPointData( self%density, WRITE_UNIT_1 )
		call WriteFieldToVTKPointData( self%velocity, WRITE_UNIT_1)
		call WriteFieldToVTKPointData( self%divergence, WRITE_UNIT_1 )
		if ( allocated(self%tracers) ) then
			do i = 1, size(self%tracers)
				call WriteFieldToVTKPointData(self%tracers(i), WRITE_UNIT_1)
			enddo
		endif
		
		call WriteFaceAreaToVTKCellData( self%mesh%faces, self%mesh%particles, WRITE_UNIT_1)
	close(WRITE_UNIT_1)
end subroutine


!
!----------------
! private methods
!----------------
!
subroutine newPrivate( self, meshSeed, initNest, maxNest, amrLimit, sphereRadius, hasDivergence )
	type(TransportMesh), intent(out) :: self
	integer(kint), intent(in) :: meshSeed
	integer(kint), intent(in) :: initNest
	integer(kint), intent(in) :: maxNest
	integer(kint), intent(in) :: amrLimit
	real(kreal), intent(in) :: sphereRadius
	logical(klog), intent(in) :: hasDivergence
	
	if ( .NOT. logInit ) call InitLogger(log, procRank)
	
	if ( .NOT. ( meshSeed == ICOS_TRI_SPHERE_SEED .OR. meshSeed == CUBED_SPHERE_SEED ) ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL, trim(logKey), " new transport mesh ERROR: invalid meshSeed")
		return
	endif
	
	call New(self%mesh, meshSeed, initNest, maxNest, amrLimit, sphereRadius )
	
	call New( self%density, 1, self%mesh%particles%N_Max, "density", "kg/m^3")
	call New( self%divergence, 1, self%mesh%particles%N_Max, "divergence", "1/s")
	call New( self%velocity, 3, self%mesh%particles%N_Max, "velocity", "m/s")
	
	call New(self%mpiParticles, self%mesh%particles%N, numProcs)
	
	self%radius = sphereRadius
	self%nonzeroDivergence = hasDivergence
end subroutine

subroutine deletePrivate(self)
	type(TransportMesh), intent(inout) :: self
	integer(kint) :: i
	
	call Delete(self%mpiParticles)
	call Delete(self%velocity)
	call Delete(self%divergence)
	call Delete(self%density)
	call Delete(self%mesh)
	if ( allocated(self%tracers) ) then
		do i = 1, size(self%tracers)
			call Delete(self%tracers(i))
		enddo
		deallocate(self%tracers)
	endif
end subroutine 

subroutine copyPrivate(self, other)
	type(TransportMesh), intent(inout) :: self
	type(TransportMesh), intent(in) :: other
	integer(kint) :: i
	
	call Copy(self%mesh, other%mesh)
	call Copy(self%density, other%density)
	call Copy(self%divergence, other%divergence)
	call Copy(self%velocity, other%velocity)
	call Copy(self%mpiParticles, other%mpiParticles)
	self%radius = other%radius
	
	if ( allocated(other%tracers)) then
		if ( allocated(self%tracers)) then
			if ( size(self%tracers) /= size(other%tracers) ) then
				call LogMessage(log, ERROR_LOGGING_LEVEL, logkey, &
					" copy TransportMesh ERROR : tracers not allocated properly.")
			else
				do i = 1, size(other%tracers)
					call Copy(self%tracers(i), other%tracers(i))
				enddo
			endif
		else
			call LogMessage(log, WARNING_LOGGING_LEVEL, logKey, " copy TransportMesh WARNING : tracers not copied.")
		endif
	endif
end subroutine

subroutine logStatsPrivate(self, aLog )
	type(TransportMesh), intent(in) :: self
	type(Logger), intent(inout) :: aLog
	integer(kint) :: i
	
	call LogMessage(aLog, TRACE_LOGGING_LEVEL, "TransportMesh ", "stats : ")
	call LogStats(self%mesh, aLog)
	call LogStats(self%density, aLog)
	call LogStats(self%divergence, aLog)
	call LogStats(self%velocity, aLog)
	if ( allocated(self%tracers)) then
		do i = 1, size(self%tracers)
			call LogStats(self%tracers(i), aLog)
		enddo
	endif	
end subroutine

subroutine setVelocityFromFunction( self, velFn, t )
	type(TransportMesh), intent(inout) :: self
	procedure(vectorFnOf3DSpaceAndTime) :: velFn
	real(kreal), intent(in) :: t
	integer(kint) :: i
	real(kreal), dimension(3) :: vel
	
	self%velocity%N = self%mesh%particles%N
	do i = 1, self%mesh%particles%N
		vel = velFn( self%mesh%particles%x(i), self%mesh%particles%y(i), self%mesh%particles%z(i), t)
		self%velocity%xComp(i) = vel(1)
		self%velocity%yComp(i) = vel(2)
		self%velocity%zComp(i) = vel(3)
	enddo
end subroutine

subroutine setScalarTracerOnMesh(self, tracerID, tracerFn)
	type(TransportMesh), intent(inout) :: self
	integer(kint), intent(in) :: tracerID
	procedure(scalarFnOf3DSpace) :: tracerFn
	integer(kint) :: i
	
	self%tracers(tracerID)%N = self%mesh%particles%N
	do i = 1, self%tracers(tracerID)%N
		self%tracers(tracerID)%scalar(i) = tracerFn( self%mesh%particles%x0(i), self%mesh%particles%y0(i), &
			self%mesh%particles%z0(i) )
	enddo
end subroutine

subroutine setVectorTracerOnMesh(self, tracerID, tracerFn )
	type(TransportMesh), intent(inout) :: self
	integer(kint), intent(in) :: tracerID
	procedure(vectorFnOf3DSpace) :: tracerFn
	integer(kint) :: i
	real(kreal), dimension(3) :: vec
	
	self%tracers(tracerID)%N = self%mesh%particles%N
	do i = 1, self%tracers(tracerID)%N
		vec = tracerFn( self%mesh%particles%x0(i), self%mesh%particles%y0(i), self%mesh%particles%z0(i) )
		self%tracers(tracerID)%xComp(i) = vec(1)
		self%tracers(tracerID)%yComp(i) = vec(2)
		self%tracers(tracerID)%zComp(i) = vec(3)
	enddo
end subroutine

function FluidMass(self)
	real(kreal) :: FluidMass
	type(TransportMesh), intent(in) :: self
	integer(kint) :: i
	
	FluidMass = 0.0_kreal
	do i = 1, self%mesh%particles%N
		if ( self%mesh%particles%isActive(i) ) then
			FluidMass = FluidMass + self%density%scalar(i) * self%mesh%particles%area(i)
		endif
	enddo
end function

function TracerMass(self, tracerID)
	real(kreal) :: TracerMass
	type(TransportMesh), intent(in) :: self
	integer(kint), intent(in) :: tracerID
	integer(kint) :: i
	
	TracerMass = 0.0_kreal
	do i = 1, self%mesh%particles%N
		if ( self%mesh%particles%isActive(i) ) then
			TracerMass = TracerMass + self%tracers(tracerID)%scalar(i) * self%density%scalar(i) * &
					 self%mesh%particles%area(i)
		endif
	enddo
end function

!> @brief Initializes a logger for the BVE module
!> 
!> Output is controlled both by message priority and by MPI Rank
!> @param aLog Target Logger object
!> @param rank Rank of this processor
subroutine InitLogger(aLog,rank)
	type(Logger), intent(out) :: aLog
	integer(kint), intent(in) :: rank
	write(logKey,'(A,A,I0.2,A)') trim(logKey),'_',rank,' : '
	if ( rank == 0 ) then
		call New(aLog,logLevel)
	else
		call New(aLog,ERROR_LOGGING_LEVEL)
	endif
	logInit = .TRUE.
end subroutine


!> @}
end module