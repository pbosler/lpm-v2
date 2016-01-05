module SphereSWEModule

use NumberKindsModule
use OutputWriterModule
use LoggerModule
use ParticlesModule
use EdgesModule
use FacesModule
use PolyMesh2dModule
use FieldModule
use MPISetupModule
use PSEDirectSumModule
use SphereGeomModule

implicit none

include 'mpif.h'

private

!
!----------------
! Module types and public declarations
!----------------
!
public SWEMesh, New, Delete, Copy
public AddTracers
public LogStats
public OutputToVTK
public SetInitialHOnMesh, SetInitialDivergenceOnMesh, SetInitialPotVortOnMesh
public SetInitialVelocityOnMesh, SetInitialRelVortOnMesh
public SetBottomHeightOnMesh

type SWEMesh
	type(PolyMesh2d) :: mesh
	type(Field) :: potVort
	type(Field) :: relVort
	type(Field) :: divergence
	type(Field) :: velocity
	type(Field) :: h
	type(FIeld) :: hBottom
	type(Field), pointer :: tracers(:) => null()
	real(kreal) :: radius
	real(kreal) :: rotationRate = 0.0_kreal
	real(kreal) :: g = 1.0_kreal
	real(kreal) :: pseEps = 0.0_kreal
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



!
!----------------
! Logging
!----------------
!
logical(klog), save :: logInit = .FALSE.
type(Logger) :: log
character(len=28), save :: logKey = 'SWE'
integer(kint), parameter :: logLevel = DEBUG_LOGGING_LEVEL

contains

!
!----------------
! public methods
!----------------
!
subroutine newPrivate( self, meshSeed, initNest, maxNest, amrLimit, sphereRadius, rotationRate, g)
	type(SWEMesh), intent(out) :: self
	integer(kint), intent(in) :: meshSeed
	integer(kint), intent(in) :: initNest
	integer(kint), intent(in) :: maxNest
	integer(kint), intent(in) :: amrLimit
	real(kreal), intent(in) :: sphereRadius
	real(kreal), intent(in) :: rotationRate
	real(kreal), intent(in) :: g
	
	if ( .NOT. logInit ) call InitLogger(log, procRank)
	
	if ( meshSeed /= ICOS_TRI_SPHERE_SEED .AND. meshSeed /= CUBED_SPHERE_SEED ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL, trim(logKey), "new SWEMesh ERROR : invalid meshSeed.")
		return
	endif
	
	call New(self%mesh, meshSeed, initNest, maxNest, amrLimit, sphereRadius)
	
	call New(self%potVort, 1, self%mesh%particles%N_Max, "potVort", "1/s")
	call New(self%relVort, 1, self%mesh%particles%N_Max, "relVort", "1/s")
	call New(self%divergence, 1, self%mesh%particles%N_Max, "divergence", "1/s")
	call New(self%h, 1, self%mesh%particles%N_Max, "h", "m")
	call New(self%hBottom, 1, self%mesh%particles%N_Max, "hBottom", "m")
	call New(self%velocity, 3, self%mesh%particles%N_Max, "velocity", "m/s")
	
	call New(self%mpiParticles, self%mesh%particles%N, numProcs)
	
	self%radius = sphereRadius
	self%rotationRate = rotationRate
	self%g = g
end subroutine

subroutine deletePrivate( self )
	type(SWEMesh), intent(inout) :: self
	integer(kint) :: i
	
	call Delete(self%mpiParticles)
	call Delete(self%velocity)
	call Delete(self%h)
	call Delete(self%divergence)
	call Delete(self%relVort)
	call Delete(self%potVort)
	call Delete(self%mesh)
	call Delete(self%hBottom)
	if ( associated(self%tracers)) then
		do i = 1, size(self%tracers)
			call Delete(self%tracers(i))
		enddo
		deallocate(self%tracers)
	endif
end subroutine

subroutine copyPrivate(self, other)
	type(SWEMesh), intent(inout) :: self
	type(SWEMesh), intent(in) :: other
	integer(kint) :: i
	
	call Copy(self%mesh, other%mesh)
	call Copy(self%h, other%h)
	call Copy(self%relVort, other%relVort)
	call Copy(self%potVort, other%potVort)
	call Copy(self%divergence, other%divergence)
	call Copy(self%velocity, other%velocity)
	call Copy(self%hBottom, other%hBottom)
	self%radius = other%radius
	self%rotationRate = other%rotationRate
	self%g = other%g
	call Copy(self%mpiParticles, other%mpiParticles)
	if ( associated(other%tracers) ) then
		if ( associated(self%tracers)) then
			if ( size(self%tracers) < size(other%tracers) ) then
				call LogMessage(log, ERROR_LOGGING_LEVEL, logkey, " copy SWEMesh ERROR : tracers not allocated properly.")
			else
				do i = 1, size(other%tracers)
					call Copy(self%tracers(i), other%tracers(i))
				enddo
			endif
		else
			allocate(self%tracers(size(other%tracers)))
			do i = 1, size(other%tracers)
				call New(self%tracers(i), other%tracers(i)%nDim, other%mesh%particles%N)
				call Copy(self%tracers(i), other%tracers(i))
			enddo
		endif
	endif
end subroutine

subroutine AddTracers(self, nTracers, tracerDims)
	type(SWEMesh), intent(inout) :: self
	integer(kint), intent(in) :: nTracers
	integer(kint), intent(in), dimension(:) :: tracerDims
	integer(kint) :: i
	
	if ( size(tracerDims) /= nTracers ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL, trim(logkey)//" AddTracers ERROR : ", &
			" must specify dimension of each tracer field.")
		return
	endif
	
	allocate(self%tracers(nTracers))
	do i = 1, nTracers
		call New(self%tracers(i), tracerDims(i), self%mesh%particles%N_Max)
	enddo	
end subroutine

subroutine logStatsPrivate( self, aLog )
	type(SWEMesh), intent(in) :: self
	type(Logger), intent(inout) :: aLog
	!
	integer(kint) :: i
	
	call StartSection(aLog, "SWEMesh ", "stats:")
	call LogMessage(aLog, TRACE_LOGGING_LEVEL, "rotationRate = ", self%rotationRate)
	call LogMessage(aLog, TRACE_LOGGING_LEVEL, "gravity = ", self%g )
	call LogMessage(aLog, TRACE_LOGGING_LEVEL, "radius = ", self%radius )
	call LogMessage(aLog, TRACE_LOGGING_LEVEL, "pseEps = ", self%pseEps)
	call LogStats(self%mesh, aLog)
	call LogStats(self%velocity, aLog )
	call LogStats(self%relVort, aLog )
	call LogStats(self%potVort, aLog )
	call LogStats(self%h, aLog )
	call LogStats(self%divergence, aLog )
	call LogStats(self%hBottom, aLog )
	if ( associated(self%tracers) ) then
		do i = 1, size(self%tracers)
			call LogStats(self%tracers(i), aLog)
		enddo
	endif
	call EndSection(aLog)
end subroutine

subroutine OutputToVTK( self, filename )
	type(SWEMesh), intent(in) :: self
	character(len=*), intent(in) :: filename
	!
	integer(kint) :: i, writeStat
	
	open(unit=WRITE_UNIT_1, file=filename, status='REPLACE', action='WRITE', iostat = writeStat)
		if ( writeStat /= 0 ) then
			call LogMessage(log, ERROR_LOGGING_LEVEL, trim(logkey)//" OutputToVTK ERROR writing to file = ", trim(filename))
			return
		endif
		!
		!	write points and topology
		!
		call WriteVTKPoints(self%mesh%particles, WRITE_UNIT_1)
		call WriteFacesToVTKPolygons(self%mesh%faces, WRITE_UNIT_1)
		!
		!	write vtk point data
		!
		call WriteVTKPointDataSectionHeader(WRITE_UNIT_1, self%mesh%particles%N)
		call WriteVTKLagCoords( self%mesh%particles, WRITE_UNIT_1)
		call WriteFieldToVTKPointData( self%potVort, WRITE_UNIT_1)
		call WriteFieldToVTKPointData( self%relVort, WRITE_UNIT_1)
		call WriteFieldToVTKPointData( self%divergence, WRITE_UNIT_1)
		call WriteFieldToVTKPointData( self%h, WRITE_UNIT_1)
		call WriteFieldToVTKPointData( self%hBottom, WRITE_UNIT_1)
		call WriteFieldToVTKPointData( self%velocity, WRITE_UNIT_1)
		if ( associated( self%tracers ) ) then
			do i = 1, size(self%tracers)
				call WriteFieldToVTKPointData( self%tracers(i), WRITE_UNIT_1)
			enddo
		endif
		!
		!	write vtk cell data
		!
		call WriteFaceAreaToVTKCellData( self%mesh%faces, self%mesh%particles, WRITE_UNIT_1)
		
	close(WRITE_UNIT_1)
end subroutine

subroutine SetInitialRelVortOnMesh( self, vortFn )
	type(SWEMesh), intent(inout) :: self
	procedure(scalarFnOf3DSpace) :: vortFn
	!
	integer(kint) :: i
	
	do i = 1, self%mesh%particles%N
		call InsertScalarToField(self%relVort, &
				vortFn( self%mesh%particles%x(i), self%mesh%particles%y(i), self%mesh%particles%z(i)))
	enddo
end subroutine

subroutine SetScalarTracerOnMesh( self, tracerID, tracerFn)
	type(SWEMesh), intent(inout) :: self
	integer(kint), intent(in) :: tracerID
	procedure(scalarFnOf3DSpace) :: tracerFn
	integer(kint) :: i
	call SetFieldToZero(self%tracers(tracerID))
	self%tracers(tracerId)%N = self%mesh%particles%N
	do i = 1, self%mesh%particles%N
		self%tracers(tracerID)%scalar(i) = &
			tracerFn( self%mesh%particles%x0(i), self%mesh%particles%y0(i), self%mesh%particles%z0(i))
	enddo
end subroutine

subroutine SetVectorTracerOnMesh( self, tracerID, tracerFn )
	type(SWEMesh), intent(inout) :: self
	integer(kint), intent(in) :: tracerID
	procedure(vectorFnOf3DSpace) :: tracerFn
	integer(kint) :: i
	real(kreal), dimension(3) :: vec
	call SetFieldToZero(self%tracers(tracerID))
	self%tracers(tracerID)%N = self%mesh%particles%N
	do i = 1, self%mesh%particles%N
		vec = tracerFn( self%mesh%particles%x0(i), self%mesh%particles%y0(i), self%mesh%particles%z0(i) )
		self%tracers(tracerID)%xComp(i) = vec(1)
		self%tracers(tracerID)%yComp(i) = vec(2)
		self%tracers(tracerID)%zComp(i) = vec(3)
	enddo
end subroutine

subroutine SetInitialHOnMesh( self, hFn )
	type(SWEMesh), intent(inout) :: self
	procedure(scalarFnOf3DSpace) :: hFn
	integer(kint) :: i
	do i = 1, self%mesh%particles%N
		call InsertScalarToField( self%h, &
			hFn( self%mesh%particles%x(i), self%mesh%particles%y(i), self%mesh%particles%z(i)))
	enddo
end subroutine

subroutine SetBottomHeightOnMesh( self, topographyFn)
	type(SWEMesh), intent(inout) :: self
	procedure(scalarFnOf3DSpace) :: topographyFn
	!
	integer(kint) :: i
	
	self%hBottom%N = self%mesh%particles%N
	do i = 1, self%mesh%particles%N
		self%hBottom%scalar(i) = topographyFn( self%mesh%particles%x(i), self%mesh%particles%y(i), self%mesh%particles%z(i))
	enddo
end subroutine

subroutine SetInitialDivergenceOnMesh( self, divFn )
	type(SWEMesh), intent(inout) :: self
	procedure(scalarFnOf3DSpace) :: divFn
	integer(kint) :: i
	do i = 1, self%mesh%particles%N
		call InsertScalarToField( self%divergence, &
			divFn( self%mesh%particles%x(i), self%mesh%particles%y(i), self%mesh%particles%z(i)))
	enddo
end subroutine

subroutine SetInitialVelocityOnMesh( self, velFn )
	type(SWEMesh), intent(inout) :: self
	procedure(vectorFnOf3DSpace) :: velFn
	integer(kint) :: i
	do i = 1, self%mesh%particles%N
		call InsertVectorToField( self%velocity, &
			velFn( self%mesh%particles%x(i), self%mesh%particles%y(i), self%mesh%particles%z(i)))
	enddo
end subroutine

subroutine SetInitialPotVortOnMesh( self )
	type(SWEMesh), intent(inout) :: self
	integer(kint) :: i
	real(kreal) :: pv, xi(3)
	do i = 1, self%mesh%particles%N
		xi = PhysCoord(self%mesh%particles, i)
		pv = (self%relVort%scalar(i) + 2.0_kreal * self%rotationRate * sin(Latitude(xi)) ) / self%h%scalar(i)
		call InsertScalarToField(self%potVort, pv)
	enddo
end subroutine

!
!----------------
! private methods
!----------------
!
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

end module