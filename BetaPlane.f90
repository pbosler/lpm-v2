module BetaPlaneMeshModule
!> @file BetaPlane.f90
!> Data structure for  a planar incompressible flow with a background rotation
!> @author Peter Bosler, Sandia National Laboratories Center for Computing Research
!> 
!>
!> @defgroup BetaPlane BetaPlane
!> Data structure for  a planar incompressible flow with a background rotation
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
use PlaneGeomModule

implicit none

include 'mpif.h'

private
!
!----------------
! Module types and public declarations
!----------------
!
public BetaPlaneMesh, New, Delete, Copy
public AddTracers
public LogStats
public OutputToVTK
public SetInitialVorticityOnMesh, SetVelocityOnMesh, SetStreamFunctionsOnMesh
public TotalKE, TotalEnstrophy
public MaxCirculationMagnitudePerFace
public ResetPeriodicBounds

type BetaPlaneMesh
	type(PolyMesh2D) :: mesh
	type(Field) :: absVort
	type(Field) :: relVort
	type(Field) :: relStream
	type(Field) :: absStream
	type(Field) :: velocity
	type(Field), dimension(:), allocatable :: tracers
	real(kreal) :: f0 = 0.0_kreal
	real(kreal) :: beta = 0.0_kreal
	type(MPISetup) :: mpiParticles
	logical(klog) :: useAMR = .FALSE.
	
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

interface ResetPeriodicBounds
	module procedure resetPeriodicBoundsMesh
	module procedure resetPeriodicBoundsArray
end interface


!
!----------------
! Logging
!----------------
!
logical(klog), save :: logInit = .FALSE.
type(Logger) :: log
character(len=28), save :: logKey = 'BetaPlane'
integer(kint), parameter :: logLevel = DEBUG_LOGGING_LEVEL
character(len=MAX_STRING_LENGTH) :: logString

contains
!
!----------------
! public methods
!----------------
!
subroutine newPrivate( self, initNest, maxNest, amrLimit, baseLat, rotRate )
	type(BetaPlaneMesh), intent(out) :: self
	integer(kint), intent(in) :: initNest
	integer(kint), intent(in) :: maxNest
	integer(kint), intent(in) :: amrLimit
	real(kreal), intent(in) :: baseLat
	real(kreal), intent(in), optional :: rotRate
	!
	integer(kint), parameter :: meshSeed = BETA_PLANE_SEED
	real(kreal), parameter :: ampFactor = 1.0_kreal
	
	if ( .NOT. logInit) call InitLogger(log, procRank)
	
	call New(self%mesh, meshSeed, initNest, maxNest, amrLimit, ampFactor)
	
	if ( abs(baseLat) > PI/2.0_kreal) then
		call LogMessage(log, ERROR_LOGGING_LEVEL, trim(logKey)//" NewBetaPlaneMesh ERROR : ", "invalid baseLat")
	else
		if ( present(rotRate) ) then
			self%f0 = 2.0_kreal * rotRate * sin( baseLat )
			self%beta = 2.0_kreal * rotRate * cos( baseLat )
		else
			self%f0 = 2.0_kreal * (2.0_kreal * PI) * sin( baseLat )
			self%beta = 2.0_kreal * (2.0_kreal * PI) * cos( baseLat )
		endif
	endif
	
	call New(self%absVort, 1, self%mesh%particles%N_Max, "absVort", "1/time")
	call New(self%relVort, 1, self%mesh%particles%N_Max, "relVort", "1/time")
	call New(self%absStream, 1, self%mesh%particles%N_Max, "absStream", "area/time")
	call New(self%relStream, 1, self%mesh%particles%N_Max, "relStream", "area/time")
	call New(self%velocity, 2, self%mesh%particles%N_Max, "velocity", "distance/time")
	
	call New(self%mpiParticles, self%mesh%particles%N,  numProcs)	
end subroutine

subroutine deletePrivate(self)
	type(BetaPlaneMesh), intent(inout) :: self
	!
	integer(kint) :: i
	
	call Delete(self%mpiParticles)
	call Delete(self%velocity)
	call Delete(self%relStream)
	call Delete(self%absStream)
	call Delete(self%relVort)
	call Delete(self%absVort)
	call Delete(self%mesh)
	
	if ( allocated(self%tracers) ) then
		do i = 1, size(self%tracers)
			call Delete(self%tracers(i))
		enddo
		deallocate(self%tracers)
	endif
end subroutine

subroutine copyPrivate( self, other )
	type(BetaPlaneMesh), intent(inout) :: self
	type(BetaPlaneMesh), intent(in) :: other
	
	self%f0 = other%f0
	self%beta = other%beta
	call Copy(self%mesh, other%mesh)
	call Copy(self%absVort, other%absVort)
	call Copy(self%relVort, other%relVort)
	call Copy(self%absStream, other%absStream)
	call Copy(self%relStream, other%relStream)
	call Copy(self%velocity, other%velocity)
	call Copy(self%mpiParticles, other%mpiParticles)
end subroutine

subroutine AddTracers(self, nTracers, tracerDims)
	type(BetaPlaneMesh), intent(inout) :: self
	integer(kint), intent(in) :: nTracers
	integer(kint), dimension(:), intent(in) :: tracerDims
	!
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
	type(BetaPlaneMesh), intent(in) :: self
	type(Logger), intent(inout) :: aLog
	call LogMessage(aLog, TRACE_LOGGING_LEVEL, "BVEMesh ", "stats : ")
	call LogStats(self%mesh, aLog)
	call LogStats(self%velocity, aLog)
	call LogStats(self%relVort, aLog)
	call LogStats(self%absVort, aLog)
	call LogStats(self%relStream, aLog)
	call LogStats(self%absStream, aLog)
end subroutine

subroutine OutputToVTK(self, filename)
	type(BetaPlaneMesh), intent(in) :: self
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

subroutine SetInitialVorticityOnMesh( self, relVortFn )
	type(BetaPlaneMesh), intent(inout) :: self
	procedure(scalarFnof2DSpace) :: relVortFn
	!
	integer(kint) :: i
	real(kreal) :: zeta
	
	self%relVort%N = self%mesh%particles%N
	self%absVort%N = self%mesh%particles%N
	do i = 1, self%mesh%particles%N
		zeta = relVortFn( self%mesh%particles%x0(i), self%mesh%particles%y0(i) )
		self%relVort%scalar(i) = zeta
		self%absVort%scalar(i) = zeta + self%f0 + self%beta * self%mesh%particles%y0(i)
	enddo	
end subroutine

subroutine setVelocityFromFunction(self, velFn )
	type(BetaPlaneMesh), intent(inout) :: self
	procedure(vectorFnOf2DSpace) :: velFn
	!
	integer(kint) :: i
	real(kreal), dimension(2) :: vel
	
	self%velocity%N = self%mesh%particles%N
	do i = 1, self%mesh%particles%N
		vel = velFn( self%mesh%particles%x(i), self%mesh%particles%y(i) )
		self%velocity%xComp(i) = vel(1)
		self%velocity%yComp(i) = vel(2)
	enddo
end subroutine

subroutine SetScalarTracerOnMesh( self, tracerID, tracerFn )
	type(BetaPlaneMesh), intent(inout) :: self
	integer(kint), intent(in) :: tracerID
	procedure(scalarFnof2DSpace) :: tracerFn
	!
	integer(kint) :: I
	
	call SetFieldToZero(self%tracers(tracerID))
	self%tracers(tracerID)%N = self%mesh%particles%N
	do i = 1, self%mesh%particles%N
		self%tracers(tracerID)%scalar(i) = tracerFn( self%mesh%particles%x0(i), self%mesh%particles%y0(i))
	enddo
end subroutine

subroutine SetVectorTracerOnMesh( self, tracerID, tracerFn)
	type(BetaPlaneMesh), intent(inout) :: self
	integer(kint), intent(in) :: tracerID
	procedure(vectorFnOf2DSpace) :: tracerFn
	!
	integer(kint) :: i
	real(kreal), dimension(2) :: vec
	
	call SetFieldToZero( self%tracers(tracerID))
	self%tracers(tracerID)%N = self%mesh%particles%N
	do i = 1, self%mesh%particles%N
		vec = tracerFn(self%mesh%particles%x0(i), self%mesh%particles%y0(i))
		self%tracers(tracerID)%xComp(i) = vec(1)
		self%tracers(tracerID)%yComp(i) = vec(2)
	enddo
end subroutine

function TotalKE( self )
	real(kreal) :: TotalKE
	type(BetaPlaneMesh), intent(in) :: self
	!
	integer(kint) :: i
	real(kreal) :: magSq
	
	TotalKE = 0.0_kreal
	do i = 1, self%mesh%particles%N
		if ( self%mesh%particles%isActive(i) ) then
			magSq = self%velocity%xCOmp(i)**2 + self%velocity%yComp(i)**2
			TotalKE = TotalKE + magSq * self%mesh%particles%area(i)
		endif
	enddo
	TotalKE = 0.5_kreal * TotalKE
end function

function TotalEnstrophy(self)
	real(kreal) :: TotalEnstrophy
	type(BetaPlaneMesh), intent(in) :: self
	!
	integer(kint) :: i
	
	TotalEnstrophy = 0.0_kreal
	do i = 1, self%mesh%particles%N
		if ( self%mesh%particles%isActive(i) ) then
			TotalEnstrophy = TotalEnstrophy + self%relVort%scalar(i)**2 * self%mesh%particles%area(i)
		endif
	enddo
	TotalEnstrophy = 0.5_kreal * TotalEnstrophy
end function

function MaxCirculationMagnitudePerFace( self )
	real(kreal) :: MaxCirculationMagnitudePerFace
	type(BetaPlaneMesh), intent(in) :: self
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

subroutine setVelocityFromVorticity( self )
	type(BetaPlaneMesh), intent(inout) :: self
	!
	integer(kint) :: i, j, mpiErrCode
	real(kreal) :: strength
	real(kreal), dimension(3) :: xi, xj
	
	self%velocity%N = self%mesh%particles%N
	do i = self%mpiParticles%indexStart(procRank), self%mpiParticles%indexEnd(procRank)
		self%velocity%xComp(i) = 0.0_kreal
		self%velocity%yComp(i) = 0.0_kreal
		xi = PhysCoord(self%mesh%particles, i)
		do j = 1, i - 1
			if ( self%mesh%particles%isActive(j) ) then
				xj = PhysCoord(self%mesh%particles, j)
				strength = 0.5_kreal * self%relVort%scalar(j) * self%mesh%particles%area(j) / &
					( cosh( 2.0_kreal * PI * (xi(2) - xj(2) )) - cos( 2.0_kreal * PI * ( xi(1) - xj(1) )) )
				self%velocity%xComp(i) = self%velocity%xComp(i) - sinh( 2.0_kreal * PI * ( xi(2) - xj(2) ) ) * strength 
				self%velocity%yComp(i) = self%velocity%yComp(i) + sin( 2.0_kreal * PI * ( xi(1) - xj(1)) ) * strength
			endif
		enddo
		do j = i+1, self%mesh%particles%N
			if ( self%mesh%particles%isActive(j) ) then
				xj = PhysCoord(self%mesh%particles, j)
				strength = 0.5_kreal * self%relVort%scalar(j) * self%mesh%particles%area(j) / &
					( cosh( 2.0_kreal * PI * (xi(2) - xj(2) )) - cos( 2.0_kreal * PI * ( xi(1) - xj(1) )) )
				self%velocity%xComp(i) = self%velocity%xComp(i) - sinh( 2.0_kreal * PI * ( xi(2) - xj(2) ) ) * strength 
				self%velocity%yComp(i) = self%velocity%yComp(i) + sin( 2.0_kreal * PI * ( xi(1) - xj(1)) ) * strength
			endif
		enddo
	enddo
	
	do i = 0, numProcs - 1
		call MPI_BCAST(self%velocity%xComp(self%mpiParticles%indexStart(i):self%mpiParticles%indexEnd(i)), &
					   self%mpiParticles%messageLength(i), MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)
		call MPI_BCAST(self%velocity%yComp(self%mpiParticles%indexStart(i):self%mpiParticles%indexEnd(i)), &
					   self%mpiParticles%messageLength(i), MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)
	enddo
end subroutine

subroutine SetStreamFunctionsOnMesh( self ) 
	type(BetaPlaneMesh), intent(inout) :: self
	!
	integer(kint) :: i, j, mpiErrCode
	real(kreal) :: greensKernel
	real(kreal), dimension(3) :: xi, xj
	
	self%absStream%N = self%mesh%particles%N
	self%relStream%N = self%mesh%particles%N
	do i = self%mpiParticles%indexStart(procRank), self%mpiParticles%indexEnd(procRank)
		self%absStream%scalar(i) = 0.0_kreal
		self%relStream%scalar(i) = 0.0_kreal
		xi = PhysCoord(self%mesh%particles, i)
		do j = 1, i - 1
			if ( self%mesh%particles%isActive(j) ) then
				xj = PhysCoord(self%mesh%particles, j)
				greensKernel = dlog( cosh(2.0_kreal * PI * ( xi(2) - xj(2) )) - cos( 2.0_kreal * PI * (xi(1) - xj(1)))) / &
					(4.0_kreal * PI )
				self%relStream%scalar(i) = self%relStream%scalar(i) + greensKernel * &
					self%relVort%scalar(j) * self%mesh%particles%area(j)
				self%absStream%scalar(i) = self%absStream%scalar(i) + greensKernel * &
					self%absVort%scalar(j) * self%mesh%particles%area(j)
			endif
		enddo
		do j = i + 1, self%mesh%particles%N
			if ( self%mesh%particles%isActive(j) ) then
				xj = PhysCoord(self%mesh%particles, j)
				greensKernel = dlog( cosh(2.0_kreal * PI * ( xi(2) - xj(2) )) - cos( 2.0_kreal * PI * (xi(1) - xj(1)))) / &
					(4.0_kreal * PI )
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

subroutine resetPeriodicBoundsArray( x )
	real(kreal), dimension(:), intent(inout) :: x
	!
	integer(kint) :: i
	
	do i = 1, size(x)
		if ( x(i) < 0.0_kreal ) then
			x(i) = x(i) + 1.0_kreal
		endif
		if ( x(i) > 1.0_kreal ) then
			x(i) = x(i) - 1.0_kreal
		endif
	enddo
end subroutine

subroutine resetPeriodicBoundsMesh( self )
	type(BetaPlaneMesh), intent(inout) :: self
	!
	integer(kint) :: i
	
	do i = 1, self%mesh%particles%N
		if ( self%mesh%particles%x(i) < 0.0_kreal ) then
			self%mesh%particles%x(i) = self%mesh%particles%x(i) + 1.0_kreal
		endif
		if ( self%mesh%particles%x(i) > 1.0_kreal ) then
			self%mesh%particles%x(i) = self%mesh%particles%x(i) - 1.0_kreal
		endif
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

!> @}
end module