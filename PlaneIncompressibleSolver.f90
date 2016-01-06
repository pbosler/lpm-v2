module PlanarIncompressibleSolverModule
!> @file PlaneIncompressibleSolver.f90
!> Data structure for integrating a planar incompressible flow forward in time using a fourth-order Runge-Kutta method.
!> @author Peter Bosler, Sandia National Laboratories Center for Computing Research
!> 
!>
!> @defgroup PlanarIncompressibleSolver PlanarIncompressibleSolver
!> Data structure for integrating a planar incompressible flow forward in time using a fourth-order Runge-Kutta method.
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
use PlanarIncompressibleModule

implicit none

include 'mpif.h'

private
!
!----------------
! Module types and public declarations
!----------------
!
public PlaneSolver, New, Delete
public Timestep

type PlaneSolver
	real(kreal), allocatable, dimension(:) :: xStart !< starting x-coordinate of each particle
	real(kreal), allocatable, dimension(:) :: yStart !< starting y-coordinate of each particle
	real(kreal), allocatable, dimension(:) :: u !< x-component of velocity for each particle
	real(kreal), allocatable, dimension(:) :: v !< y-component of velocity for each particle
	real(kreal), allocatable, dimension(:) :: vort !< vorticity carried by each particle
	real(kreal), allocatable, dimension(:) :: area !< area represented by each particle (passive particles represent zero area)
	logical(klog), allocatable, dimension(:) :: mask !< mask(i) is .TRUE. if particle i is active
	
	real(kreal), allocatable, dimension(:) :: xIn !< x-coordinate input to RK4
	real(kreal), allocatable, dimension(:) :: xStage1 !< x-coordinates of each particle at RK4 stage 1
	real(kreal), allocatable, dimension(:) :: xStage2 !< x-coordinates of each particle at RK4 stage 2
	real(kreal), allocatable, dimension(:) :: xStage3 !< x-coordinates of each particle at RK4 stage 3
	real(kreal), allocatable, dimension(:) :: xStage4 !< x-coordinates of each particle at RK4 stage 4
	real(kreal), allocatable, dimension(:) :: yIn !< y-coordinate input to RK4
	real(kreal), allocatable, dimension(:) :: yStage1 !< y-coordinates of each particle at RK4 stage 1
	real(kreal), allocatable, dimension(:) :: yStage2 !< y-coordinates of each particle at RK4 stage 2
	real(kreal), allocatable, dimension(:) :: yStage3 !< y-coordinates of each particle at RK4 stage 3
	real(kreal), allocatable, dimension(:) :: yStage4 !< y-coordinates of each particle at RK4 stage 4
	
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

interface Timestep
	module procedure timestepPrivate
end interface

!
!----------------
! Logging
!----------------
!
logical(klog), save :: logInit = .FALSE.
type(Logger) :: log
character(len=28), save :: logKey = 'PlaneIncSolver'
integer(kint), parameter :: logLevel = DEBUG_LOGGING_LEVEL
character(len=MAX_STRING_LENGTH) :: logString

contains

!
!----------------
! public methods
!----------------
!

!> @brief Allocates memory for a new solver.  Unlike a mesh, which maintains extra space in memory for possible adaptive
!> refinement, solvers only allocate memory based on the current number of particles. Hence, a new solver needs to be 
!> created every time a mesh is refined.
!> 
!> @param[out] self Target solver
!> @param[in] plane @ref PlanarIncompressible mesh
subroutine newPrivate( self, plane )
	type(PlaneSolver), intent(out) :: self
	type(PlaneMeshIncompressible), intent(in) :: plane
	!
	integer(kint) :: nParticles
	
	if ( .NOT. logInit ) call InitLogger( log, procRank)
	
	nParticles = plane%mesh%particles%N
	
	allocate(self%xStart(nParticles))
	allocate(self%yStart(nParticles))
	allocate(self%area(nParticles))
	allocate(self%u(nParticles))
	allocate(Self%v(nParticles))
	allocate(self%mask(nParticles))
	allocate(self%vort(nParticles))
	
	allocate(self%xIn(nParticles))
	allocate(self%xStage1(nParticles))
	allocate(self%xStage2(nParticles))
	allocate(self%xStage3(nParticles))
	allocate(self%xStage4(nParticles))
	allocate(self%yIn(nParticles))
	allocate(self%yStage1(nParticles))
	allocate(self%yStage2(nParticles))
	allocate(self%yStage3(nParticles))
	allocate(self%yStage4(nParticles))
	
	self%xStart = plane%mesh%particles%x(1:nParticles)
	self%yStart = plane%mesh%particles%y(1:nParticles)
	self%u = plane%velocity%xComp(1:nParticles)
	self%v = plane%velocity%yComp(1:nParticles)
	self%area = plane%mesh%particles%area(nParticles)
	self%vort = plane%vorticity%scalar(nParticles)
	self%mask = plane%mesh%particles%isActive(nParticles)
end subroutine

!> @brief Deletes and frees memory associated with a PlanarIncompressible solver.
!> @param[inout] self solver
subroutine deletePrivate( self ) 
	type(PlaneSolver), intent(inout) :: self
	
	if ( allocated(self%xIn) ) then
		deallocate(Self%xStart)
		deallocate(self%yStart)
		deallocate(self%u)
		deallocate(self%v)
		deallocate(self%vort)
		deallocate(self%area)
		deallocate(self%mask)
		deallocate(self%xIn)
		deallocate(self%xStage1)
		deallocate(self%xStage2)
		deallocate(self%xStage3)
		deallocate(self%xStage4)
		deallocate(self%yIn)
		deallocate(self%yStage1)
		deallocate(self%yStage2)
		deallocate(self%yStage3)
		deallocate(self%yStage4)
	endif
end subroutine

!> @brief Advances a PlanarIncompressible mesh forward in time by one timestep using 4th order Runge-Kutta.
!> @param[inout] self PlanarIncompressible solver
!> @param[inout] plane PlanarIncompressible mesh
!> @param[in] dt timestep increment
subroutine timestepPrivate( self, plane, dt ) 
	type(PlaneSolver), intent(inout) :: self
	type(PlaneMeshIncompressible), intent(inout) :: plane
	real(kreal), intent(in) :: dt
	!
	integer(kint) :: i, nParticles
	
	nParticles = plane%mesh%particles%N
	
	!
	! RK Stage 1
	!
	do i = 1, nParticles
		self%xStage1(i) = dt * self%u(i)
		self%yStage1(i) = dt * self%v(i)
	enddo

	!
	! RK Stage 2
	!
	do i = 1, nParticles
		self%xIn(i) = self%xStart(i) + 0.5_kreal * self%xStage1(i)
		self%yIn(i) = self%yStart(i) + 0.5_kreal * self%yStage1(i)
	enddo
	
	call planarIncompressibleVelocity( self%xStage2, self%yStage2, self%xIn, self%yIn, &
			 self%vort, self%area, self%mask, plane%mpiParticles)
	
	do i = 1, nParticles
		self%xStage2(i) = dt * self%xStage2(i)
		self%yStage2(i) = dt * self%yStage2(i)
	enddo
			 
	!
	! RK Stage 3
	!
	do i = 1, nParticles
		self%xIn(i) = self%xStart(i) + 0.5_kreal * self%xStage2(i)
		self%yIn(i) = self%yStart(i) + 0.5_kreal * self%yStage2(i)
	enddo
	
	call planarIncompressibleVelocity( self%xStage3, self%yStage3, self%xIn, self%yIn, &
			 			 self%vort, self%area, self%mask, plane%mpiParticles)

	
	do i = 1, nParticles
		self%xStage3(i) = dt * self%xStage3(i)
		self%yStage3(i) = dt * self%yStage3(i)
	enddo
	
	!
	! RK Stage 4
	!
	do i = 1, nParticles
		self%xIn(i) = self%xStart(i) + self%xStage3(i)
		self%yIn(i) = self%yStart(i) + self%yStage3(i)
	enddo	
	
	call planarIncompressibleVelocity( self%xStage4, self%yStage4, self%xIn, self%yIn, &
			 			 self%vort, self%area, self%mask, plane%mpiParticles)

	
	do i = 1, nParticles
		self%xStage4(i) = dt * self%xStage4(i)
		self%yStage4(i) = dt * self%yStage4(i)
	enddo
	
	!
	! RK update
	!
	do i = 1, nParticles
		self%xStart(i) = self%xStart(i) + self%xStage1(i) / 6.0_kreal + self%xStage2(i)/3.0_kreal + &
						 self%xStage3(i) / 3.0_kreal + self%xStage4(i) / 6.0_kreal
		self%yStart(i) = self%yStart(i) + self%yStage1(i) / 6.0_kreal + self%yStage2(i)/3.0_kreal + &
						 self%yStage3(i) / 3.0_kreal + self%yStage4(i) / 6.0_kreal
	enddo
	
	plane%mesh%particles%x(1:nParticles) = self%xStart
	plane%mesh%particles%y(1:nParticles) = self%yStart
	
!	call SetVelocityOnMesh( plane )
!	self%u = plane%velocity%xComp(1:nParticles)
!	self%v = plane%velocity%yComp(1:nParticles)
	
	call planarIncompressibleVelocity( self%u, self%v, self%xStart, self%yStart, self%vort, &
		self%area, self%mask, plane%mpiParticles)
	
	plane%velocity%xComp(1:nParticles) = self%u
	plane%velocity%yComp(1:nParticles) = self%v
	
	call SetStreamFunctionOnMesh( plane )
end subroutine

!
!----------------
! private methods
!----------------
!

!> @brief Computes velocity based on the given vorticity of all point vortices (particles) using the Biot-Savart integral.
!> Integral is computed in parallel as a direct summation.
!>
!> @param[out] u x-component of velocity at each particle
!> @param[out] v y-component of velocity at each particle
!> @param[in] x x-coordinate of each particle
!> @param[in] y y-coordinate of each particle
!> @param[in] vortIn vorticity carried by each particle
!> @param[in] areaIn area carried by each particle
!> @param[in] activeMask .TRUE. for active particles, .FALSE. for passive particles
!> @param[in] mpiParticles @ref MPISetup object to distribute integral over MPI processes
subroutine planarIncompressibleVelocity( u, v, x, y, vortIn, areaIn, activeMask, mpiParticles )
	real(kreal), dimension(:), intent(out) :: u
	real(kreal), dimension(:), intent(out) :: v
	real(kreal), dimension(:), intent(in) :: x
	real(kreal), dimension(:), intent(in) :: y
	real(kreal), dimension(:), intent(in) :: vortIn
	real(kreal), dimension(:), intent(in) :: areaIn
	logical(klog), dimension(:), intent(in) :: activeMask
	type(MPISetup), intent(in) :: mpiParticles
	!
	integer(kint) :: i, j, mpiErrCode
	real(kreal) :: strength
	
	do i = mpiParticles%indexStart(procRank), mpiParticles%indexEnd(procRank)
		u(i) = 0.0_kreal
		v(i) = 0.0_kreal
		do j = 1, i - 1
			if ( activeMask(j) ) then
				strength = vortIn(j) * areaIn(j) / ( 2.0_kreal * PI * ( (x(i) - x(j))**2 + (y(i) - y(j))**2))
				u(i) = u(i) - ( y(i) - y(j) ) * strength
				v(i) = v(i) + ( x(i) - x(j) ) * strength
			endif
		enddo
		do j = i + 1, size(x)
			if ( activeMask(j) ) then
				strength = vortIn(j) * areaIn(j) / ( 2.0_kreal * PI * ( (x(i) - x(j))**2 + (y(i) - y(j))**2))
				u(i) = u(i) - ( y(i) - y(j) ) * strength
				v(i) = v(i) + ( x(i) - x(j) ) * strength
			endif
		enddo
	enddo
	
	do i = 0, numProcs - 1
		call MPI_BCAST( u(mpiParticles%indexStart(i):mpiParticles%indexEnd(i)), mpiParticles%messageLength(i), &
				MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)
		call MPI_BCAST( v(mpiParticles%indexStart(i):mpiParticles%indexEnd(i)), mpiParticles%messageLength(i), &
				MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)
	enddo
end subroutine

!> @brief Initializes a logger for the PlanarIncompressible solver module
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