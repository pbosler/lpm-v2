module PlanarIncompressibleModule
!> @file PlanarIncompressible.f90
!> Data structure for solving the fluid equations for two-dimensional (planar) inviscid, incompressible flow.
!> @author Peter Bosler, Sandia National Laboratories Center for Computing Research
!>
!>
!> @defgroup PlanarIncompressible PlanarIncompressible
!> @brief Data structure for solving the fluid equations for two-dimensional (in the plane) inviscid, incompressible flow.
!> 
!> Inviscid planar incompressible flow may be characterized completely by the flow's vorticity and stream function.
!> This module provides a data structure for discretizing the flow using a collection of point vortices.
!> Each point vortex carries vorticity and is advected by the velocity induced by all other vortices.   
!>
!> The spatial domain is discretized in space and time by a @ref PolyMesh2d object.
!> Associated @ref Field objects for vorticity, velocity, and the stream function are provided, as well as a set 
!> of optional @ref Field objects for passive tracers carried by the flow.  
!>
!> In addition to the variables carried by the base @ref Particles object (e.g., physical and Lagrangian coordinates), 
!> the PlaneMeshIncompressible data type adds fields for vorticity @f$ \zeta(\vec{x},t) @f$, which is materially conserved, 
!> velocity @f$ \vec{u}(\vec{x},t) @f$, and the stream function @f$ \psi(\vec{x},t) @f$.
!> These variables are related by the equations
!> @f{align*}{
!>	\frac{\partial \zeta}{\partial t} + \vec{u}\cdot\nabla \zeta & = \frac{D\zeta}{Dt} = 0, \\
!>  \nabla^2 \psi & = - \zeta, \\
!>  \zeta &= \nabla \times \vec{u}.
!> @f}
!> These equations are solved using vortex methods (particles = point vortices).  For more information on these 
!> numerical techniques, see
!> * A. J. Chorin and J. E. Marsden, _A Mathematical Introduction to Fluid Mechanics_, 3rd edition, Springer, 2000.
!> * G.-F. Cottet and P. D. Koumoutsakos, _Vortex Methods_, Cambridge University Press, 2000.
!> * A. J. Majda and A. L. Bertozzi, _Vorticity and Incompressible Flow_, Cambridge University Press, 2002.
!>
!> Currently, only free boundary conditions are supported.  
!> The flow is integrated in time using a @ref PlanarIncompressibleSolver object.
!> 
!> The data structure is designed to use a replicated data algorithm (via the @ref MPISetup module) for efficiency in 
!> parallel computing environments.
!> 
!> Output routines suitable for reading/dispaly with [Paraview](www.paraview.org) are provided. 
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
public PlaneMeshIncompressible, New, Delete, Copy
public AddTracers
public LogStats
public OutputToVTK
public SetInitialVorticityOnMesh, SetVelocityOnMesh, SetStreamFunctionOnMesh
public TotalKE, TotalEnstrophy
public MaxCirculationPerFace
public SetAbsoluteTolerances

type PlaneMeshIncompressible
	type(PolyMesh2d) :: mesh !< @ref PolyMesh2d for spatial discretization
	type(Field) :: vorticity !< scalar @ref Field
	type(Field) :: streamFn !< scalar @ref Field
	type(Field) :: velocity !< vector @ref Field
	type(Field), dimension(:), allocatable :: tracers  !< allocatable array of scalar and vector @ref Field objects
	type(MPISetup) :: mpiParticles !< @ref MPISetup to distribute all particles over the available MPI ranks
	logical(klog) :: useAMR = .FALSE. !< .TRUE. if mesh will be adaptively refined
	
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

!> @brief Allocates memory and initializes a PlanarIncompressible mesh with accompanying @ref Field objects for physical variables.  
!> Tracers must be added separately, if desired.
!> 
!> @param[out] self PlanarIncompressible mesh
!> @param[in] initNest initial recursion level for @ref PolyMesh2d 
!> @param[in] maxNest maximum recursion level for @ref PolyMesh2d
!> @param[in] meshSeed mesh seed integer as defined by @ref NumberKinds
!> @param[in] amrLimit maximum number of refinements beyond initNest allowed
!> @param[in] meshRadius maximum spatial extent of mesh from origin
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

!> @brief Deletes and frees memory associated with a PlanarIncompressible mesh
!> @param self Planar Incompressible mesh
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

!> @brief Performs a deep copy of one PlanarIncompressible mesh into another.  Both objects must have already been allocated.
!> @param[out] self Target PlanarIncompressible mesh
!> @param[in] other Source PlanarIncompressible mesh
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

!> @brief Adds memory for passive tracers to a PlanarIncompressible mesh.
!> @param[inout] self PlanarIncompressible mesh
!> @param[in] nTracers number of tracers to add
!> @param[in] tracerDims array of values (either 1 or 2) corresponding to the dimensions of each tracer.  Scalar tracers will have dimension 1, vector tracers will have dimension 2.
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

!> @brief Output basic information about a PlanarIncompressible mesh to a @ref Logger
!> @param[in] self PlanarIncompressible mesh
!> @param[inout] aLog @ref Logger
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

!> @brief Sets absolute tolerances for AMR criteria based on an initially uniform mesh with corresponding vorticity distribution.
!> @param[in] self PlanarIncompressible mesh
!> @param[inout] circTol On input, relative tolerance (between 0 and 1).  On output, absolute tolerance for circulation magnitude about a face.
!> @param[inout] lagVarTol On input, relative tolerance (between 0 and 1).  On output, absolute tolerance for Lagrangian variation magnitude on each face.
subroutine SetAbsoluteTolerances( self, circTol, lagVarTol )
	type(PlaneMeshIncompressible), intent(in) :: self
	real(kreal), intent(inout) :: circTol
	real(kreal), intent(inout) :: lagVarTol
	!
	circTol = circTol * MaxCirculationPerFace(self)
	lagVarTol = lagVarTol * MaxLagrangianVariationPerFace(self%mesh)
end subroutine

!> @brief Writes a PlanarIncompressible mesh to a legacy formatted .vtk file, including all @ref Field data for variables and tracers.
!> @param[in] self PlanarIncompressible mesh
!> @param[in] filename Name of output file
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

!> @brief Defines an initial vorticity distribution on a PlanarIncompressible mesh.
!> 
!> @param[inout] self PlanarIncompressible mesh
!> @param[in] vortFn Vorticity distribution function.  Must have same interface as numberkindsmodule::scalarFnOf2DSpace
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

!> @brief Defines an initial velocity distribution on a PlanarIncompressible mesh.
!> 
!> @param[inout] self PlanarIncompressible mesh
!> @param[in] velFn Velocity distribution function.  Must have same interface as numberkindsmodule::vectorFnOf2DSpace
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

!> @brief Defines an initial scalar tracer distribution on a PlanarIncompressible mesh.
!> 
!> @param[inout] self PlanarIncompressible mesh
!> @param[in] tracerID index of tracer in `PlanarIncompressiblMesh%%tracers(:)`
!> @param[in] tracerFn Tracer distribution function.  Must have same interface as numberkindsmodule::scalarFnOf2DSpace
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

!> @brief Defines an initial scalar tracer distribution on a PlanarIncompressible mesh.
!> 
!> @param[inout] self PlanarIncompressible mesh
!> @param[in] tracerID index of tracer in `PlanarIncompressiblMesh%%tracers(:)`
!> @param[in] tracerFn Tracer distribution function.  Must have same interface as numberkindsmodule::vectorFnOf2DSpace
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

!> @brief Computes the total kinetic energy (a conserved integral) on a PlanarIncompressible mesh.
!> @param[in] self PlanarIncompressible mesh
!> @return totalKE total kinetic energy
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

!> @brief Computes the total enstrophy (a conserved integral) on a PlanarIncompressible mesh.
!> @param[in] self PlanarIncompressible mesh
!> @return totalEnstrophy
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

!> @brief Computes the maximum circulation magnitude about each face and returns the maximum value for the whole mesh.
!> @param[in] self PlanarIncompressible mesh
!> @return MaximumCirculation magnitude
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

!> @brief Defines the velocity distribution on a PlanarIncompressible mesh using the Biot-Savart integral.
!> @param[inout] self PlanarIncompressible mesh
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

!> @brief Defines the stream function on a PlanarIncompressible mesh using the Green's function integral.
!> @param[inout] self PlanarIncompressible mesh
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
!> @brief Initializes a logger for the PlanarIncompressible module
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
