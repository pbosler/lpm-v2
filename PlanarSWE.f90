module PlanarSWEModule

use NumberKindsModule
use OutputWriterModule
use LoggerModule
use ParticlesModule
use EdgesModule, only : MaxEdgeLength
use FacesModule
use PolyMesh2dModule
use FieldModule
use MPISetupModule
use PSEDirectSumModule

implicit none

include 'mpif.h'

private
public SWEMesh, New, Delete
public AddTracers
public LogStats
public OutputToVTKFile
public SetBottomHeightOnMesh, SetInitialDivergenceOnMesh, SetInitialHOnMesh, SetInitialPotVortOnMesh
public SetInitialVelocityOnMesh, SetInitialVorticityOnMesh

type SWEMesh
	type(PolyMesh2d) :: mesh
	type(Field) :: relVort
	type(Field) :: potVort
	type(Field) :: divergence
	type(Field) :: velocity
	type(Field) :: h
	type(Field) :: hBottom
	type(Field), pointer :: tracers(:) => null()
	real(kreal) :: f0 = 0.0_kreal
	real(kreal) :: beta = 0.0_kreal
	real(kreal) :: g = 0.0_kreal
	real(kreal) :: pseEps = 0.0_kreal
	type(MPISetup) :: mpiParticles
	contains
		final :: deletePrivate	
end type

interface New
	module procedure newPrivate
end interface

interface Copy
	module procedure copyPrivate
end interface

interface Delete
	module procedure deletePrivate
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

!interface
!	real(8) function scalarFnOf2DSpace( x, y)
!		real(8), intent(in) :: x, y
!	end function
!end interface
!
!interface
!	function vectorFnOf2DSpace( x, y)
!		real(8), dimension(2) :: vectorFnOf2DSpace
!		real(8), intent(in) :: x, y
!	end function
!end interface



contains
!
!----------------
! public methods
!----------------
!

subroutine newPrivate(self, meshSeed, initNest, maxNest, amrLimit, meshRadius, f0, beta, g )
	type(SWEMesh), intent(out) :: self
	integer(kint), intent(in) :: meshSeed
	integer(kint), intent(in) :: initNest
	integer(kint), intent(in) :: maxNest
	integer(kint), intent(in) :: amrLimit
	real(kreal), intent(in) :: meshRadius
	real(kreal), intent(in) :: f0, beta, g
	
	if ( .NOT. logInit ) call InitLogger(log, procRank)
	
	if ( meshSeed /= TRI_HEX_SEED .AND. meshSeed /= QUAD_RECT_SEED ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL, trim(logKey), "new SWEMesh ERROR : invalid meshSeed.")
		return
	endif
	
	call New(self%mesh, meshSeed, initNest, maxNest, amrLimit, meshRadius )
	
	call New(self%relVort, 1, self%mesh%particles%N, "relVort", "1/s")
	call New(self%potVort, 1, self%mesh%particles%N, "potVort", "1/s")
	call New(self%divergence, 1, self%mesh%particles%N, "divergence", "1/s")	
	call New(self%velocity, 2, self%mesh%particles%N, "velocity", "m/s")
	call New(self%h, 1, self%mesh%particles%N, "h", "m")
	call New(self%hBottom, 1, self%mesh%particles%N, "hB", "m")
	call New(self%mpiParticles, self%mesh%particles%N, numProcs)
	
	self%f0 = f0
	self%beta = beta
	self%g = g
end subroutine

subroutine SWEComputeVelocity(self, velocity, x, y, relVort, divergence, area )
	type(SWEMesh), intent(inout) :: self
	type(Field), intent(inout) :: velocity
	real(kreal), intent(in) :: x(:), y(:), relVort(:), divergence(:), area(:)
	!
	integer(kint) :: i, j, mpiErrCode
	real(kint) :: denom
	
	call SetFieldToZero(velocity)
	velocity%N = self%mesh%particles%N
	
	do i = self%mpiParticles%indexStart(procRank), self%mpiParticles%indexEnd(procRank)
		do j = 1, self%mesh%particles%N
			if ( self%mesh%particles%isActive(j) ) then
				denom = 2.0_kreal * PI * ( (x(i) - x(j))**2 + (y(i) - y(j))**2)
				velocity%xComp(i) = velocity%xComp(i) + (( y(i) - y(j) ) * relVort(j) + (x(i) - x(j)) * divergence(j)) * &
					area(j) / denom
				velocity%yComp(i) = velocity%yComp(i) + (-(x(i) - x(j) ) * relVort(j) + (y(i) - y(j)) * divergence(j)) * &
					area(j) / denom
			endif
		enddo
	enddo
	
	do i = 0, numProcs - 1
		call MPI_BCAST(velocity%xComp(self%mpiParticles%indexStart(i):self%mpiParticles%indexEnd(i)), &
				self%mpiParticles%messageLength(i), MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)
		call MPI_BCAST(velocity%yComp(self%mpiParticles%indexStart(i):self%mpiParticles%indexEnd(i)), &
				self%mpiParticles%messageLength(i), MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)				
	enddo
end subroutine

subroutine copyPrivate(self, other)
	type(SWEMesh), intent(inout) :: self
	type(SWEMesh), intent(in) :: other
end subroutine

subroutine deletePrivate(self)
	type(SWEMesh), intent(inout) :: self
	!
	integer(kint) :: i
	
	call Delete(self%mpiParticles)
	call Delete(self%hBottom)
	call Delete(self%relVort)
	call Delete(self%potVort)
	call Delete(self%divergence)
	call Delete(self%velocity)
	call Delete(self%h)
	call Delete(self%mesh)
	if ( associated(self%tracers)) then
		do i = 1, size(self%tracers)
			call Delete(self%tracers(i))
		enddo
		deallocate(self%tracers)
	endif
end subroutine

subroutine logStatsPrivate(self, aLog)
	type(SWEMesh), intent(in) :: self
	type(Logger), intent(inout) :: aLog
	call LogMessage(aLog, TRACE_LOGGING_LEVEL, "SWE Mesh "," stats : ")
	call LogStats(self%mesh, aLog)
	call LogStats(self%h, aLog)
	call LogMessage(aLog, TRACE_LOGGING_LEVEL, "Coriolis f0 = ", self%f0)
	call LogMessage(aLog, TRACE_LOGGING_LEVEL, "Coriolis beta = ", self%beta)
	call LogStats(self%relVort, aLog)
	call LogStats(self%potVort, aLog)
	call LogStats(self%divergence, aLog)
	call LogStats(self%velocity, aLog) 
	call LogStats(self%hBottom, aLog)
end subroutine

subroutine SetCoriolis(self, f0, beta)
	type(SWEMesh), intent(inout) :: self
	real(kreal), intent(in) :: f0, beta
	self%f0 = f0
	self%beta = beta
end subroutine

subroutine AddTracers( self, nTracers, tracerDims )
	type(SWEMesh), intent(inout) :: self
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
		call New(self%tracers(i), tracerDims(i), self%mesh%particles%N)
	enddo	
end subroutine

subroutine OutputToVTKFile( self, filename )
	type(SWEMesh), intent(in) :: self
	character(len=*), intent(in) :: filename
	!
	integer(kint) :: i, writeStat
	
	open(unit=WRITE_UNIT_1, file=filename, status="REPLACE", action="WRITE", iostat=writeStat)
		if ( writeStat /= 0 ) then
			call LogMessage(log, ERROR_LOGGING_LEVEL, trim(logkey)//" OutputToVTK ERROR writing to file = ", trim(filename))
			return
		endif
		!
		!	write points and topology
		!
		call WriteVTKFileHeader(WRITE_UNIT_1)
		write(WRITE_UNIT_1, '(A,I8,A)') "POINTS ", self%mesh%particles%N, " double"
		do i = 1, self%mesh%particles%N
			write(WRITE_UNIT_1,*) self%mesh%particles%x(i), self%mesh%particles%y(i), &
								  self%h%scalar(i) + self%hBottom%scalar(i)
		enddo
		call WriteFacesToVTKPolygons( self%mesh%faces, WRITE_UNIT_1)
		!
		!	write particle data
		!
		call WriteVTKPointDataSectionHeader(WRITE_UNIT_1, self%mesh%particles%N)
		call WriteVTKLagCoords(self%mesh%particles, WRITE_UNIT_1) ! always has to be first point data call
		call WriteFieldToVTKPointData(self%h, WRITE_UNIT_1)
		call WriteFieldToVTKPointData(self%relVort, WRITE_UNIT_1)
		call WriteFieldToVTKPointData(self%potVort, WRITE_UNIT_1)
		call WriteFieldToVTKPointData(self%divergence, WRITE_UNIT_1)
		call WriteFieldToVTKPointData(self%velocity, WRITE_UNIT_1)
		call WriteFieldToVTKPointData(self%hBottom, WRITE_UNIT_1)
		if ( associated(self%tracers)) then
			do i = 1, size(self%tracers)
				call WriteFieldToVTKPointData(self%tracers(i), WRITE_UNIT_1)
			enddo
		endif
		!
		!	write face data
		!
		call WriteFaceAreaToVTKCellData(self%mesh%faces, self%mesh%particles, WRITE_UNIT_1) ! always has to be first face data call
	close(WRITE_UNIT_1)
end subroutine

subroutine OutputToToMatlabMFile( self, filename )
	type(SWEMesh), intent(in) :: self
	character(len=*), intent(in) :: filename
	!
	integer(kint) :: writeStat
	open(unit=WRITE_UNIT_1, file=filename, status="REPLACE", action="WRITE", iostat=writeStat)
		if ( writeStat /= 0 ) then
			call LogMessage(log, ERROR_LOGGING_LEVEL, trim(logkey)//" OutputToVTK ERROR writing to file = ", trim(filename))
			return
		endif
	close(WRITE_UNIT_1)
end subroutine

subroutine SetInitialVorticityOnMesh(self, vorticityFn )
	type(SWEMesh), intent(inout) :: self
	procedure(scalarFnOf2DSpace) :: vorticityFn
	!
	integer(kint) :: i
	
	do i = 1, self%mesh%particles%N
		call InsertScalarToField(self%relVort, vorticityFn( self%mesh%particles%x(i), self%mesh%particles%y(i) ) )
	enddo	
end subroutine

subroutine SetInitialDivergenceOnMesh(self, divergenceFn )
	type(SWEMesh), intent(inout) :: self
	procedure(scalarFnOf2DSpace) :: divergenceFn
	!
	integer(kint) :: i
	
	do i = 1, self%mesh%particles%N
		call InsertScalarToField(self%divergence, divergenceFn( self%mesh%particles%x(i), self%mesh%particles%y(i) ) )
	enddo	
end subroutine

subroutine SetInitialHOnMesh(self, hFn )
	type(SWEMesh), intent(inout) :: self
	procedure(scalarFnOf2DSpace) :: hFn
	!
	integer(kint) :: i
	
	do i = 1, self%mesh%particles%N
		call InsertScalarToField(self%h, hFn( self%mesh%particles%x(i), self%mesh%particles%y(i) ))
	enddo
end subroutine

subroutine SetBottomHeightOnMesh( self, depthFn)
	type(SWEMesh), intent(inout) :: self
	procedure(scalarFnOf2DSpace) :: depthFn
	!
	integer(kint) :: i
	
	self%hBottom%N = self%mesh%particles%N
	do i = 1, self%mesh%particles%N
		self%hBottom%scalar(i) = depthFn( self%mesh%particles%x(i), self%mesh%particles%y(i) )
	enddo		
end subroutine

subroutine SetInitialVelocityOnMesh(self, velFn)
	type(SWEMesh), intent(inout) :: self
	procedure(vectorFnOf2DSpace) :: velFn
	!
	integer(kint) :: i
	
	do i = 1, self%mesh%particles%N
		call InsertVectorToField(self%velocity, velFn( self%mesh%particles%x(i), self%mesh%particles%y(i)))
	enddo
end subroutine

subroutine SetInitialPotVortOnMesh( self )
	type(SWEMesh), intent(inout) :: self
	!
	integer(kint) :: i
	real(kreal) :: pv
	
	do i = 1, self%mesh%particles%N
		if ( self%h%scalar(i) > 0.0_kreal ) then
			pv = (self%relVort%scalar(i) + self%f0 + self%beta*self%mesh%particles%y(i)) / self%h%scalar(i)
		else
			pv = 0.0_kreal
		endif
		call InsertScalarToField( self%potVort, pv )
	enddo
end subroutine

!
!----------------
! private methods
!----------------
!

subroutine InitLogger(aLog,rank)
! Initialize a logger for this module and processor
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
