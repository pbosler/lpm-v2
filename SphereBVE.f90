module SphereBVEModule

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
public BVEMesh, New, Delete, Copy
public AddTracers
public LogStats
public OutputToVTK
public SetInitialVorticityOnMesh, SetVelocityOnMesh, SetStreamFunctionsOnMesh
public TotalKE, TotalEnstrophy
public MaxCirculationMagnitudePerFace

type BVEMesh
	type(PolyMesh2d) :: mesh
	type(Field) :: absVort
	type(Field) :: relVort
	type(Field) :: relStream
	type(Field) :: absStream
	type(Field) :: velocity
	type(Field), dimension(:), allocatable :: tracers
	real(kreal) :: radius = 1.0_kreal
	real(kreal) :: rotationRate = 0.0_kreal
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
	module procedure setVelocityFromVorticity
end interface

!
!----------------
! Logging
!----------------
!
logical(klog), save :: logInit = .FALSE.
type(Logger) :: log
character(len=28), save :: logKey = 'BVE'
integer(kint), parameter :: logLevel = DEBUG_LOGGING_LEVEL

contains
!
!----------------
! public methods
!----------------
!
subroutine newPrivate( self, meshSeed, initNest, maxNest, amrLimit, sphereRadius, rotationRate )
	type(BVEMesh), intent(out) :: self
	integer(kint), intent(in) :: meshSeed
	integer(kint), intent(in) :: initNest
	integer(kint), intent(in) :: maxNest
	integer(kint), intent(in) :: amrLimit
	real(kreal), intent(in) :: sphereRadius
	real(kreal), intent(in) :: rotationRate
	
	if ( .NOT. logInit) call InitLogger(log, procRank)
	
	if ( meshSeed /= ICOS_TRI_SPHERE_SEED .AND. meshSeed /= CUBED_SPHERE_SEED ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL, trim(logKey), "new BVEMesh ERROR : invalid meshSeed.")
		return
	endif
	
	call New(self%mesh, meshSeed, initNest, maxNest, amrLimit, sphereRadius)
	
	call New(self%absVort, 1, self%mesh%particles%N_Max, "absVort", "1/s")
	call New(self%relVort, 1, self%mesh%particles%N_Max, "relVort", "1/s")
	call New(self%absStream, 1, self%mesh%particles%N_Max, "absStreamFn", "m^2/s")
	call New(self%relStream, 1, self%mesh%particles%N_Max, "relStreamFn", "m^2/s")
	call New(self%velocity, 3, self%mesh%particles%N_Max, "velocity", "m/s")
	call New(self%mpiParticles, self%mesh%particles%N, numProcs)
	
	self%radius = sphereRadius
	self%rotationRate = rotationRate	
end subroutine

subroutine deletePrivate( self )
	type(BVEMesh), intent(inout) :: self
	integer(kint) :: i
	
	call Delete(self%mpiParticles)
	call Delete(self%velocity)
	call Delete(self%relStream)
	call Delete(self%absStream)
	call Delete(self%absVort)
	call Delete(self%relVort)
	call Delete(self%mesh)
	if ( allocated(self%tracers) ) then
		do i = 1, size(self%tracers)
			call Delete(self%tracers(i))
		enddo
		deallocate(self%tracers)
	endif
end subroutine

subroutine copyPrivate(self, other)
	type(BVEMesh), intent(inout) :: self
	type(BVEMesh), intent(in) :: other
	integer(kint) :: i
	
	call Copy(self%mesh, other%mesh)
	call Copy(self%relVort, other%relVort)
	call Copy(self%absVort, other%absVort)
	call Copy(self%velocity, other%velocity)
	call Copy(self%absStream, other%absStream)
	call Copy(self%relStream, other%relStream)
	call Copy(self%mpiParticles, other%mpiParticles)
	self%rotationRate = other%rotationRate
	self%radius = other%radius
	
	if ( allocated(other%tracers)) then
		if ( allocated(self%tracers)) then
			if ( size(self%tracers) /= size(other%tracers) ) then
				call LogMessage(log, ERROR_LOGGING_LEVEL, logkey, " copy BVEMesh ERROR : tracers not allocated properly.")
			else
				do i = 1, size(other%tracers)
					call Copy(self%tracers(i), other%tracers(i))
				enddo
			endif
		else
			call LogMessage(log, WARNING_LOGGING_LEVEL, logKey, " copy BVEMesh ERROR : tracers not copied.")
		endif
	endif
end subroutine

subroutine AddTracers(self, nTracers, tracerDims)
	type(BVEMesh), intent(inout) :: self
	integer(kint), intent(in) :: nTracers
	integer(kint), dimension(:), intent(in) :: tracerDims
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

subroutine logStatsPrivate( self, aLog)
	type(BVEMEsh), intent(in) :: self
	type(Logger), intent(inout) :: aLog
	!
	integer(kint) :: i
	
	call LogMessage(aLog, TRACE_LOGGING_LEVEL, "BVEMesh ", "stats : ")
	call LogStats(self%mesh, aLog)
	call LogStats(self%velocity, aLog)
	call LogStats(self%relVort, aLog)
	call LogStats(self%absVort, aLog)
	call LogStats(self%relStream, aLog)
	call LogStats(self%absStream, aLog)
	if ( allocated( self%tracers ) ) then
		do i = 1, size(self%tracers)
			call LogStats(self%tracers(i), aLog)
		enddo
	endif
end subroutine

subroutine OutputToVTK(self, filename)
	type(BVEMesh), intent(in) :: self
	character(len=*), intent(in) :: filename
	integer(kint) :: i, writeStat
	
	open(unit=WRITE_UNIT_1, file=filename, status='REPLACE', action='WRITE', iostat=writeStat)
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
		call WriteFieldToVTKPointData( self%relVort, WRITE_UNIT_1)
		call WriteFieldToVTKPointData( self%absVort, WRITE_UNIT_1)
		call WriteFieldToVTKPointData( self%velocity, WRITE_UNIT_1)
		call WriteFieldToVTKPointData( self%relStream, WRITE_UNIT_1)
		call WriteFieldToVTKPointData( self%absStream, WRITE_UNIT_1)
		if ( allocated( self%tracers ) ) then
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

subroutine SetInitialVorticityOnMesh(self, relVortFn )
	type(BVEMesh), intent(inout) :: self
	procedure(scalarFnOf3DSpace) :: relVortFn
	integer(kint) :: i
	real(kreal) :: zeta
	
	self%relVort%N = self%mesh%particles%N
	self%absVort%N = self%mesh%particles%N
	do i = 1, self%mesh%particles%N
		zeta = relVortFn( self%mesh%particles%x0(i), self%mesh%particles%y0(i), self%mesh%particles%z0(i) )
		!call InsertScalarToField(self%relVort, zeta )
		!call InsertScalarToField(self%absVort, zeta + 2.0_kreal*self%rotationRate*self%mesh%particles%z0(i) / self%radius)
		self%relVort%scalar(i) = zeta
		self%absVort%scalar(i) = zeta + 2.0_kreal * self%rotationRate * self%mesh%particles%z0(i) / self%radius
	enddo
end subroutine

subroutine setVelocityFromFunction(self, velFn )
	type(BVEMesh), intent(inout) :: self
	procedure(vectorFnOf3DSpace) :: velFn 
	integer(kint) :: i
	
	do i = 1, self%mesh%particles%N
		call InsertVectorToField(self%velocity, &
			velFn( self%mesh%particles%x0(i), self%mesh%particles%y0(i), self%mesh%particles%z0(i)) )
	enddo
end subroutine

subroutine SetScalarTracerOnMesh(self, tracerId, tracerFn )
	type(BVEMesh), intent(inout) :: self
	integer(kint), intent(in) :: tracerId
	procedure(scalarFnOf3DSpace) :: tracerFn
	integer(kint) :: i
	
	call SetFieldToZero(self%tracers(tracerID))
	self%tracers(tracerId)%N = self%mesh%particles%N
	do i = 1, self%mesh%particles%N
		self%tracers(tracerID)%scalar(i) = tracerFn( self%mesh%particles%x0(i), self%mesh%particles%y0(i), &
												     self%mesh%particles%z0(i) )
	enddo
end subroutine

subroutine SetVectorTracerOnMesh( self, tracerID, tracerFn )
	type(BVEMesh), intent(inout) :: self
	integer(kint), intent(in) :: tracerId
	procedure(vectorFnOf3DSpace) :: tracerFn
	integer(kint) :: i
	real(kreal), dimension(3) :: vec
	
	call SetFieldToZero(self%tracers(tracerID))
	self%tracers(tracerID)%N = self%mesh%particles%N
	do i = 1, self%mesh%particles%N
		vec = tracerFn(self%mesh%particles%x0(i), self%mesh%particles%y0(i), self%mesh%particles%z0(i))
		self%tracers(tracerID)%xComp(i) = vec(1)
		self%tracers(tracerID)%yComp(i) = vec(2)
		self%tracers(tracerID)%zComp(i) = vec(3)
	enddo
end subroutine

function TotalKE( self )
	real(kreal) :: TotaLKE
	type(BVEMesh), intent(in) :: self
	integer(kint) :: i
	real(kreal) :: magSq

	TotalKE = 0.0_kreal
	do i = 1, self%mesh%particles%N
		if ( self%mesh%particles%isActive(i) ) then
			magSq = self%velocity%xComp(i)**2 + self%velocity%yComp(i)**2 + self%velocity%zComp(i)**2
			TotalKE = TotalKE + magSq * self%mesh%particles%area(i)
		endif
	enddo
	TotalKE = 0.5_kreal * TotalKE
end function

function TotalEnstrophy( self )
	real(kreal) :: TotalEnstrophy
	type(BVEMesh), intent(in) :: self
	integer(kint) :: i
	
	TotalEnstrophy = 0.0_kreal
	do i = 1, self%mesh%particles%N
		if ( self%mesh%particles%isActive(i) ) then
			TotalEnstrophy = TotalEnstrophy + self%relVort%scalar(i)**2 * self%mesh%particles%area(i)
		endif
	enddo
	TotalEnstrophy = 0.5_kreal * TotalEnstrophy 
end function

subroutine SetStreamFunctionsOnMesh( self )
	type(BVEMesh), intent(inout) :: self
	integer(kint) :: i, j, mpiErrCode
	real(kreal) :: greensKernel 
	real(kreal), dimension(3) :: xi, xj
	
	self%relStream%n = self%mesh%particles%n
	self%absStream%n = self%mesh%particles%n
	do i = self%mpiParticles%indexStart(procRank), self%mpiParticles%indexEnd(procRank)
		self%relStream%scalar(i) = 0.0_kreal
		self%absStream%scalar(i) = 0.0_kreal
		xi = PhysCoord(self%mesh%particles, i)
		do j = 1, i - 1
			if ( self%mesh%particles%isActive(j) ) then
				xj = PhysCoord(self%mesh%particles, j)
				greensKernel = - dlog( self%radius * self%radius - sum(xi * xj)) / (4.0_kreal * PI * self%radius)
				self%relStream%scalar(i) = self%relStream%scalar(i) + greensKernel * &
					self%relVort%scalar(j) * self%mesh%particles%area(j)
				self%absStream%scalar(i) = self%absStream%scalar(i) + greensKernel * &
					self%absVort%scalar(j) * self%mesh%particles%area(j)
			endif
		enddo
		do j = i + 1, self%mesh%particles%N
			if ( self%mesh%particles%isActive(j) ) then
				xj = PhysCoord(self%mesh%particles, j)
				greensKernel = - dlog( self%radius * self%radius - sum(xi * xj)) / (4.0_kreal * PI * self%radius)
				self%relStream%scalar(i) = self%relStream%scalar(i) + greensKernel * &
					self%relVort%scalar(j) * self%mesh%particles%area(j)
				self%absStream%scalar(i) = self%absStream%scalar(i) + greensKernel * & 
					self%absVort%scalar(j) * self%mesh%particles%area(j)
			endif
		enddo
	enddo
	
	do i = 0, numProcs - 1
		call MPI_BCAST(self%relStream%scalar(self%mpiParticles%indexStart(i):self%mpiParticles%indexEnd(i)), &
			self%mpiParticles%messageLength(i), MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)
		call MPI_BCAST(self%absStream%scalar(self%mpiParticles%indexStart(i):self%mpiParticles%indexEnd(i)), &
			self%mpiParticles%messageLength(i), MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)			
	enddo
end subroutine

subroutine setVelocityFromVorticity( self )
	type(BVEMesh), intent(inout) :: self
	integer(kint) :: i, j, mpiErrCode
	real(kreal) :: strength
	real(kreal), dimension(3) :: xi, xj
	
	self%velocity%n = self%mesh%particles%N
	do i = self%mpiParticles%indexStart(procRank), self%mpiParticles%indexEnd(procRank)
		self%velocity%xComp(i) = 0.0_kreal
		self%velocity%yComp(i) = 0.0_kreal
		self%velocity%zComp(i) = 0.0_kreal
		xi = PhysCoord(self%mesh%particles, i)
		do j = 1, i - 1
			if ( self%mesh%particles%isActive(j) ) then
				xj = PhysCoord(self%mesh%particles, j)
				strength = -self%relVort%scalar(j)*self%mesh%particles%area(j) / &
					(4.0_kreal * PI * self%radius * (self%radius * self%radius - sum( xi * xj )))
				self%velocity%xComp(i) = self%velocity%xComp(i) + (xi(2)*xj(3) - xi(3)*xj(2)) * strength 
				self%velocity%yComp(i) = self%velocity%yComp(i) + (xi(3)*xj(1) - xi(1)*xj(3)) * strength 
				self%velocity%zComp(i) = self%velocity%zComp(i) + (xi(1)*xj(2) - xi(2)*xj(1)) * strength
			endif
		enddo
		do j = i + 1, self%mesh%particles%N
			if ( self%mesh%particles%isActive(j) ) then
				xj = PhysCoord(self%mesh%particles, j)
				strength = -self%relVort%scalar(j)*self%mesh%particles%area(j) / &
					(4.0_kreal * PI * self%radius * (self%radius * self%radius - sum( xi * xj )))
				self%velocity%xComp(i) = self%velocity%xComp(i) + (xi(2)*xj(3) - xi(3)*xj(2)) * strength
				self%velocity%yComp(i) = self%velocity%yComp(i) + (xi(3)*xj(1) - xi(1)*xj(3)) * strength 
				self%velocity%zComp(i) = self%velocity%zComp(i) + (xi(1)*xj(2) - xi(2)*xj(1)) * strength
			endif
		enddo
	enddo
	
	do i = 0, numProcs - 1
		call MPI_BCAST(self%velocity%xComp(self%mpiParticles%indexStart(i):self%mpiParticles%indexEnd(i)), &
					   self%mpiParticles%messageLength(i), MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)
		call MPI_BCAST(self%velocity%yComp(self%mpiParticles%indexStart(i):self%mpiParticles%indexEnd(i)), &
					   self%mpiParticles%messageLength(i), MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)
		call MPI_BCAST(self%velocity%zComp(self%mpiParticles%indexStart(i):self%mpiParticles%indexEnd(i)), &
					   self%mpiParticles%messageLength(i), MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)					   					   
	enddo
end subroutine

function MaxCirculationMagnitudePerFace( self )
	real(kreal) :: MaxCirculationMagnitudePerFace
	type(BVEMesh), intent(in) :: self
	!
	integer(kint) :: i
	real(kreal) :: testCirc
	integer(kint) :: particleIndex
	
	MaxCirculationMagnitudePerFace = 0.0_kreal
	do i = 1, self%mesh%faces%N
		if ( .NOT. self%mesh%faces%hasChildren(i) ) then
			particleIndex = self%mesh%faces%centerParticle(i)
			testCirc = abs( self%relVort%scalar(particleIndex) ) * self%mesh%particles%area(particleIndex)
			if ( testCirc > MaxCirculationMagnitudePerFace ) MaxCirculationMagnitudePerFace = testCirc
		endif
	enddo
end function

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

