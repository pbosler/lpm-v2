module PlanarIncompressibleModule

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
public PlaneMeshIncompressible, New, Delete, Copy
public AddTracers
public LogStats
public OutputToVTK
public SetInitialVorticityOnMesh, SetVelocityOnMesh, SetStreamFunctionOnMesh
public TotalKE, TotalEnstrophy
public MaxCirculationPerFace
public SetAbsoluteTolerances

type PlaneMeshIncompressible
	type(PolyMesh2d) :: mesh
	type(Field) :: vorticity
	type(Field) :: streamFn
	type(Field) :: velocity
	type(Field), dimension(:), allocatable :: tracers
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

!
!----------------
! Logging
!----------------
!
logical(klog), save :: logInit = .FALSE.
type(Logger) :: log
character(len=28), save :: logKey = 'PlaneIncomp'
integer(kint), parameter :: logLevel = DEBUG_LOGGING_LEVEL
character(len=MAX_STRING_LENGTH) :: logString

contains
!
!----------------
! public methods
!----------------
!
subroutine newPrivate( self, initNest, maxNest, meshSeed, amrLimit, meshRadius )
	type(PlaneMeshIncompressible), intent(out) :: self
	integer(kint), intent(in) :: initNest
	integer(kint), intent(in) :: maxNest
	integer(kint), intent(in) :: meshSeed
	integer(kint), intent(in) :: amrLimit
	real(kreal), intent(in) :: meshRadius
	
	if ( .NOT. logInit ) call InitLogger( log, procRank )
	
	call New(self%mesh, meshSeed, initNest, maxNest, amrLimit, meshRadius)
	
	self%useAMR = (( amrLimit > 0 ) .AND. (maxNest > initNest) )
	
	call New(self%vorticity, 1, self%mesh%particles%N_Max, "vorticity", "1/time")
	call New(self%streamFn, 1, self%mesh%particles%N_Max, "streamFn", "area/time")
	call New(self%velocity, 2, self%mesh%particles%N_Max, "velocity", "dist/time")
	call New(self%mpiParticles, self%mesh%particles%N, numProcs )
end subroutine

subroutine deletePrivate(self)
	type(PlaneMeshIncompressible), intent(inout) :: self
	!
	integer(kint) :: i
	
	call Delete(self%mpiParticles)
	call Delete(self%velocity)
	call Delete(self%streamFn)
	call Delete(self%vorticity)
	call Delete(self%mesh)
	if ( allocated(self%tracers)) then
		do i = 1, size(self%tracers)
			call Delete(self%tracers(i))
		enddo
		deallocate(self%tracers)
	endif
end subroutine

subroutine copyPrivate( self, other )
	type(PlaneMeshIncompressible), intent(inout) :: self
	type(PlaneMeshIncompressible), intent(in) :: other
	!
	integer(kint) :: i
	
	call Copy(self%mesh, other%mesh)
	call Copy(self%vorticity, other%vorticity)
	call Copy(self%streamFn, other%streamFn)
	call Copy(self%velocity, other%velocity)
	call Copy(self%mpiParticles, other%mpiParticles)
	if ( allocated(self%tracers) .AND. allocated(other%tracers)) then
		do i = 1, size(other%tracers)
			call Copy(self%tracers(i), other%tracers(i))
		enddo
	endif
end subroutine

subroutine AddTracers(self, nTracers, tracerDims)
	type(PlaneMeshIncompressible), intent(inout) :: self
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

subroutine logStatsPrivate(self, alog)
	type(PlaneMeshIncompressible), intent(in) :: self
	type(Logger), intent(inout) :: aLog
	
	call StartSection(aLog, "PlaneMeshIncompressible Stats:")
	call LogStats(self%mesh, alog)
	call LogStats(self%vorticity, alog)
	call LogStats(self%streamFn, alog)
	call LogStats(self%velocity, alog)
	call EndSection(aLog)
end subroutine

subroutine SetAbsoluteTolerances( self, circTol, lagVarTol )
	type(PlaneMeshIncompressible), intent(in) :: self
	real(kreal), intent(inout) :: circTol
	real(kreal), intent(inout) :: lagVarTol
	!
	circTol = circTol * MaxCirculationPerFace(self)
	lagVarTol = lagVarTol * MaxLagrangianVariationPerFace(self%mesh)
end subroutine

subroutine OutputToVTK(self, filename)
	type(PlaneMeshIncompressible), intent(in) :: self
	character(len=*), intent(in) :: filename
	!
	integer(kint) :: i, writeStat
	
	open( unit=WRITE_UNIT_1, file=filename, status='REPLACE', action='WRITE', iostat=writeStat)
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
		call WriteFieldToVTKPointData(self%vorticity, WRITE_UNIT_1)
		call WriteFieldToVTKPointData(self%streamFn, WRITE_UNIT_1)
		call WriteFieldToVTKPointData(self%velocity, WRITE_UNIT_1)
		if ( allocated(self%tracers)) then
			do i = 1, size(self%tracers)
				call WriteFieldToVTKPointData(self%tracers(i), WRITE_UNIT_1)
			enddo
		endif
		!
		!	write vtk cell data
		!
		call WriteFaceAreaToVTKCellData( self%mesh%faces, self%mesh%particles, WRITE_UNIT_1)
	close(WRITE_UNIT_1)
end subroutine

subroutine SetInitialVorticityOnMesh( self, vortFn )
	type(PlaneMeshIncompressible), intent(inout) :: self
	procedure(scalarFnOf2DSpace) :: vortFn
	!
	integer(kint) :: i
	
	self%vorticity%N = self%mesh%particles%N
	do i = 1, self%mesh%particles%N
		self%vorticity%scalar(i) = vortFn( self%mesh%particles%x0(i), self%mesh%particles%y0(i) )
	enddo
end subroutine

subroutine setVelocityFromFunction( self, velFn )
	type(PlaneMeshIncompressible), intent(inout) :: self
	procedure(vectorFnOf2DSpace) :: velFn
	!
	integer(kint) :: i
	real(kreal), dimension(2) :: velVec
	
	self%velocity%N = self%mesh%particles%N
	do i = 1, self%mesh%particles%N
		velVec = velFn( self%mesh%particles%x(i), self%mesh%particles%y(i))
		self%velocity%xComp(i) = velVec(1)
		self%velocity%yComp(i) = velVec(2)
	enddo
end subroutine

subroutine SetScalarTracerOnMesh( self, tracerID, tracerFn )
	type(PlaneMeshIncompressible), intent(inout) :: self
	integer(kint), intent(in) :: tracerID
	procedure(scalarFnOf2DSpace) :: tracerFn
	!
	integer(kint) :: i
	
	self%tracers(tracerID)%N = self%mesh%particles%N
	do i = 1, self%Mesh%particles%n
		self%tracers(tracerID)%scalar(i) = tracerFn( self%mesh%particles%x0(i), self%mesh%particles%y0(i))
	enddo
end subroutine

subroutine SetVectorTracerOnMesh( self, tracerID, tracerFn)
	type(PlaneMeshIncompressible), intent(inout) :: self
	integer(kint), intent(in) :: tracerID
	procedure(vectorFnOf2DSpace) :: tracerFn
	!
	integer(kint) :: i
	real(kreal), dimension(2) :: vecT
	
	self%tracers(tracerID)%N = self%mesh%particles%N
	do i = 1, self%mesh%particles%N
		vecT = tracerFn( self%mesh%particles%x0(i), self%mesh%particles%y0(i))
		self%tracers(tracerID)%xComp(i) = vecT(1)
		self%tracers(tracerID)%yComp(i) = vecT(2)
	enddo
end subroutine

function TotalKE( self )
	real(kreal) :: TotalKE
	type(PlaneMeshIncompressible), intent(in) :: self
	!
	integer(kint) :: i
	real(kreal) :: spdSq
	
	TotalKE = 0.0_kreal
	do i = 1, self%mesh%particles%N
		if ( self%mesh%particles%isActive(i) ) then
			spdSq = self%velocity%xComp(i)**2 + self%velocity%yComp(i)**2
			TotalKE = TotalKE + spdSq * self%mesh%particles%area(i)
		endif
	enddo
	TotalKE = 0.5_kreal * TotalKE
end function

function TotalEnstrophy(self)
	real(kreal) :: TotalEnstrophy
	type(PlaneMeshIncompressible), intent(in) :: self
	!
	integer(kint) :: i
	
	TotalEnstrophy = 0.0_kreal
	do i = 1, self%mesh%particles%N
		if ( self%mesh%particles%isActive(i) ) then
			TotalEnstrophy = TotalEnstrophy + self%vorticity%scalar(i)**2 * self%mesh%particles%area(i)
		endif
	enddo
	TotalEnstrophy = 0.5_kreal * TotalEnstrophy
end function

function MaxCirculationPerFace(self)
	real(kreal) :: MaxCirculationPerFace
	type(PlaneMeshIncompressible), intent(in) :: self
	!
	integer(kint) :: i
	real(kreal) :: testCirc
	integer(kint) :: particleIndex
	
	MaxCirculationPerFace = 0.0_kreal
	do i = 1, self%mesh%faces%N
		if ( .NOT. self%mesh%faces%hasChildren(i) ) then
			particleIndex = self%mesh%faces%centerParticle(i)
			testCirc = abs( self%vorticity%scalar(particleIndex) ) * self%mesh%particles%area(particleIndex)
			if ( testCirc > MaxCirculationPerFace) MaxCirculationPerFace = testCirc
		endif
	enddo
end function

subroutine setVelocityFromVorticity( self )
	type(PlaneMeshIncompressible), intent(inout) :: self
	!
	integer(kint) :: i, j, mpiErrCode
	real(kreal) :: xi, yi, xj, yj, strength
	
	self%velocity%N = self%mesh%particles%N
	do i = self%mpiParticles%indexStart(procRank), self%mpiParticles%indexEnd(procRank)
		self%velocity%xComp(i) = 0.0_kreal
		self%velocity%yComp(i) = 0.0_kreal
		xi = self%mesh%particles%x(i)
		yi = self%mesh%particles%y(i)
		do j = 1, i - 1
			if ( self%mesh%particles%isActive(j) ) then
				xj = self%mesh%particles%x(j)
				yj = self%mesh%particles%y(j)
				strength = self%vorticity%scalar(j) * self%mesh%particles%area(j) / &
					( 2.0_kreal * PI * ( (xi - xj)**2 + (yi - yj)**2 ))
				self%velocity%xComp(i) = self%velocity%xComp(i) - ( yi - yj ) * strength
				self%velocity%yComp(i) = self%velocity%yComp(i) + ( xi - xj ) * strength
			endif
		enddo
		do j = i + 1, self%mesh%particles%N
			if ( self%mesh%particles%isActive(j) ) then
				xj = self%mesh%particles%x(j)
				yj = self%mesh%particles%y(j)
				strength = self%vorticity%scalar(j) * self%mesh%particles%area(j) / &
					( 2.0_kreal * PI * ( (xi - xj)**2 + (yi - yj)**2 ))
				self%velocity%xComp(i) = self%velocity%xComp(i) - ( yi - yj ) * strength
				self%velocity%yComp(i) = self%velocity%yComp(i) + ( xi - xj ) * strength
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

subroutine SetStreamFunctionOnMesh(self)
	type(PlaneMeshIncompressible), intent(inout) :: self
	!
	integer(kint) :: i, j, mpiErrCode
	real(kreal) :: greensKernel, xi, xj, yi, yj
	
	self%streamFn%N = self%mesh%particles%N
	do i = self%mpiParticles%indexStart(procRank), self%mpiParticles%indexEnd(procRank)
		self%streamFn%scalar(i) = 0.0_kreal
		xi = self%mesh%particles%x(i)
		yi = self%mesh%particles%y(i)
		do j = 1, i - 1
			if ( self%mesh%particles%isActive(j)) then
				xj = self%mesh%particles%x(j)
				yj = self%mesh%particles%y(j)
				greensKernel = dlog( sqrt((xi - xj)**2 + (yi - yj)**2 ) )/(2.0_kreal * PI)
				self%streamFn%scalar(i) = self%streamFn%scalar(i) + greensKernel * self%vorticity%scalar(j) * &
					self%mesh%particles%area(j)
			endif
		enddo
		do j = i + 1, self%mesh%particles%N
			if ( self%mesh%particles%isActive(j)) then
				xj = self%mesh%particles%x(j)
				yj = self%mesh%particles%y(j)
				greensKernel = dlog( sqrt((xi - xj)**2 + (yi - yj)**2 ) )/(2.0_kreal * PI)
				self%streamFn%scalar(i) = self%streamFn%scalar(i) + greensKernel * self%vorticity%scalar(j) * &
					self%mesh%particles%area(j)
			endif
		enddo
	enddo
	
	do i = 0, numProcs - 1
		call MPI_BCAST(self%streamFn%scalar( self%mpiParticles%indexStart(i):self%mpiParticles%indexEnd(i)), &
				self%mpiParticles%messageLength(i), MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)
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
