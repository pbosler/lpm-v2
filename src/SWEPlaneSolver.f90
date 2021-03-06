module SWEPlaneSolverModule
!> @file SWEPlaneSolver.f90
!> Data structure for solving the Shallow Water Equations (SWE) in the beta plane.
!> @author Peter Bosler, Sandia National Laboratories Center for Computing Research
!>
!>
!> @defgroup SWEPlaneSolver SWEPlaneSolver
!> Data structure for solving the Shallow Water Equations (SWE) in the beta plane.
!> Used with the @ref PlanarSWE module.
!>
!> The shallow water equations are presented in the detailed description of the @ref PlanarSWE module.  @n
!>
!> Integrals are used to compute the velocity @f$\vec{u}@f$, the double dot product @f$ \nabla\vec{u}:\nabla\vec{u} @f$, and a PSE integral (see @ref PSEDirectSum)
!> approximates the Laplacian of the fluid surface.@n
!> As with the other integrals in LPM, these integrals are computed as a parallel direct summation across all MPI ranks.@n
!> See sweplanesolvermodule::sweplanerhsintegrals.
!>
!> All other terms are computed as ODEs along each particle trajectory; see sweplanesolvermodule::timestepprivate.
!>
!> @{
use NumberKindsModule
use LoggerModule
use UtilitiesModule
use ParticlesModule
use PolyMesh2dModule
use FieldModule
use MPISetupModule
use PSEDirectSumModule
use PlanarSWEModule

implicit none

include 'mpif.h'

private
public SWESolver, New, Delete
public Timestep
public LogStats

!
!----------------
! types and module variables
!----------------
!
type SWESolver
	real(kreal), allocatable :: xStart(:) !< starting x-coordinate of each particle
	real(kreal), allocatable :: yStart(:) !< starting y-coordinate of each particle
	real(kreal), allocatable :: relVortStart(:) !< starting relative vorticity of each particle
	real(kreal), allocatable :: areaStart(:) !< starting area of each particle
	real(kreal), allocatable :: divStart(:) !< starting divergence of each particle
	real(kreal), allocatable :: hStart(:) !< starting depth at each particle
	real(kreal), allocatable :: u(:) !< x-component of fluid velocity
	real(kreal), allocatable :: v(:) !< y-component of fluid velocity
	logical(klog), allocatable :: mask(:) !< mask(i) is .TRUE. if particle i is active
	real(kreal), allocatable :: potVort(:) !< potential vorticity of each particle
	real(kreal), allocatable :: doubleDot(:) !< double dot product of velocity at each particle
	real(kreal), allocatable :: lapSurf(:) !< Laplacian of fluid surface at each particle

	real(kreal), allocatable :: xIn(:)  !< x-coordinate input to RK4
	real(kreal), allocatable :: xStage1(:) !< x-coordinates of each particle at RK4 stage 1
	real(kreal), allocatable :: xStage2(:) !< x-coordinates of each particle at RK4 stage 2
	real(kreal), allocatable :: xStage3(:) !< x-coordinates of each particle at RK4 stage 3
	real(kreal), allocatable :: xStage4(:) !< x-coordinates of each particle at RK4 stage 4
	real(kreal), allocatable :: yIn(:) !< y-coordinate input to RK4
	real(kreal), allocatable :: yStage1(:) !< y-coordinates of each particle at RK4 stage 1
	real(kreal), allocatable :: yStage2(:) !< y-coordinates of each particle at RK4 stage 2
	real(kreal), allocatable :: yStage3(:) !< y-coordinates of each particle at RK4 stage 3
	real(kreal), allocatable :: yStage4(:) !< y-coordinates of each particle at RK4 stage 4
	real(kreal), allocatable :: areaIn(:) !< area input to RK4
	real(kreal), allocatable :: areaStage1(:) !< area of each particle at RK4 stage 1
	real(kreal), allocatable :: areaStage2(:) !< area of each particle at RK4 stage 2
	real(kreal), allocatable :: areaStage3(:) !< area of each particle at RK4 stage 3
	real(kreal), allocatable :: areaStage4(:) !< area of each particle at RK4 stage 4
	real(kreal), allocatable :: relVortIn(:) !< relative vorticity input to RK4
	real(kreal), allocatable :: relVortStage1(:) !< relative vorticity of each particle at RK4 stage 1
	real(kreal), allocatable :: relVortStage2(:) !< relative vorticity of each particle at RK4 stage 2
	real(kreal), allocatable :: relVortStage3(:) !< relative vorticity of each particle at RK4 stage 3
	real(kreal), allocatable :: relVortStage4(:) !< relative vorticity of each particle at RK4 stage 4
	real(kreal), allocatable :: divIn(:) !< divergence input to RK4
	real(kreal), allocatable :: divStage1(:) !< divergence of each particle at RK4 stage 1
	real(kreal), allocatable :: divStage2(:) !< divergence of each particle at RK4 stage 2
	real(kreal), allocatable :: divStage3(:) !< divergence of each particle at RK4 stage 3
	real(kreal), allocatable :: divStage4(:) !< divergence of each particle at RK4 stage 4
	real(kreal), allocatable :: hIn(:) !< depth input to RK4
	real(kreal), allocatable :: hStage1(:) !< depth of each particle at RK4 stage 1
	real(kreal), allocatable :: hStage2(:) !< depth of each particle at RK4 stage 2
	real(kreal), allocatable :: hStage3(:) !< depth of each particle at RK4 stage 3
	real(kreal), allocatable :: hStage4(:) !< depth of each particle at RK4 stage 4

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
character(len=28), save :: logKey = 'SWESolver'
integer(kint), parameter :: logLevel = DEBUG_LOGGING_LEVEL
character(len=MAX_STRING_LENGTH) :: logString

contains
!
!----------------
! public methods
!----------------
!

!> @brief Allocates memory for a new SWE Solver.  Initializes the solver to match the current mesh.
!> @param[out] self Target SWE Solver
!> @param[in] plane @ref PlanarSWE mesh
!> @param[in] topoFn bottom topography function, must conform to the numberkindsmodule::scalarFnOf2DSpace interface
subroutine newPrivate( self, plane, topoFn )
	type(SWESolver), intent(out) :: self
	type(SWEMesh), intent(in) :: plane
	procedure(scalarFnOf2dSpace) :: topoFn
	!
	integer(kint) :: i, nP

	if (.NOT. logInit ) call InitLogger(log, procRank)

	nP = plane%mesh%particles%N

	allocate(self%xStart(1:nP))
	allocate(self%yStart(1:nP))
	allocate(self%areaStart(1:nP))
	allocate(self%divStart(1:nP))
	allocate(self%relVortStart(1:nP))
	allocate(self%mask(1:nP))
	allocate(self%u(1:nP))
	allocate(self%v(1:nP))
	allocate(self%potVort(1:nP))
	allocate(self%hStart(1:nP))
	allocate(self%doubleDot(1:nP))
	allocate(self%lapSurf(1:nP))

	allocate(self%xIn(1:nP))
	allocate(self%xStage1(1:nP))
	allocate(self%xStage2(1:nP))
	allocate(self%xStage3(1:nP))
	allocate(self%xStage4(1:nP))
	allocate(self%yIn(1:nP))
	allocate(self%yStage1(1:nP))
	allocate(self%yStage2(1:nP))
	allocate(self%yStage3(1:nP))
	allocate(self%yStage4(1:nP))
	allocate(self%relVortIn(1:nP))
	allocate(self%relVortStage1(1:nP))
	allocate(self%relVortStage2(1:nP))
	allocate(self%relVortStage3(1:nP))
	allocate(self%relVortStage4(1:nP))
	allocate(self%divIn(1:nP))
	allocate(self%divStage1(1:nP))
	allocate(self%divStage2(1:nP))
	allocate(self%divStage3(1:nP))
	allocate(self%divStage4(1:nP))
	allocate(self%areaIn(1:nP))
	allocate(self%areaStage1(1:nP))
	allocate(self%areaStage2(1:nP))
	allocate(self%areaStage3(1:nP))
	allocate(self%areaStage4(1:nP))
	allocate(self%hIn(1:nP))
	allocate(self%hStage1(1:nP))
	allocate(self%hStage2(1:nP))
	allocate(self%hStage3(1:nP))
	allocate(self%hStage4(1:nP))

	self%xStart = plane%mesh%particles%x(1:nP)
	self%yStart = plane%mesh%particles%y(1:nP)
	self%areaStart = plane%mesh%particles%area(1:nP)
	self%relVortStart = plane%relVort%scalar(1:nP)
	self%potVort = plane%potVort%scalar(1:nP)
	self%divStart = plane%divergence%scalar(1:nP)
	self%hStart = plane%h%scalar(1:nP)
	self%u = plane%velocity%xComp(1:nP)
	self%v = plane%velocity%yComp(1:nP)
	self%mask = plane%mesh%particles%isActive(1:nP)

	call SWEPlaneRHSIntegrals( self%u, self%v, self%doubleDot, self%lapSurf, self%xStart, self%yStart, self%relVortStart, &
				self%divStart, self%hStart, topoFn, self%areaStart, self%mask, plane%pseEps, plane%mpiParticles)
end subroutine

subroutine logStatsPrivate(self, aLog)
	type(SWESolver), intent(in) :: self
	type(Logger), intent(inout) :: aLog

	call StartSection(aLog, trim(logKey)//" SWESolver Stats:")

	call LogMessage(aLog, TRACE_LOGGING_LEVEL, "size xStart = ", size(self%xStart))
	call LogMessage(aLog, TRACE_LOGGING_LEVEL, "max xStart = ", maxval(self%xStart))

	call LogMessage(aLog, TRACE_LOGGING_LEVEL, "size yStart = ", size(self%yStart))
	call LogMessage(aLog, TRACE_LOGGING_LEVEL, "max yStart = ", maxval(self%yStart))

	call LogMessage(aLog, TRACE_LOGGING_LEVEL, "size u = ", size(self%u))
	call LogMessage(aLog, TRACE_LOGGING_LEVEL, "max u = ", maxval(self%u))
	call LogMessage(aLog, TRACE_LOGGING_LEVEL, "size v = ", size(self%v))
	call LogMessage(aLog, TRACE_LOGGING_LEVEL, "max v = ", maxval(self%v))

	call LogMessage(aLog, TRACE_LOGGING_LEVEL, "size doubleDot = ", size(self%doubleDot))
	call LogMessage(aLog, TRACE_LOGGING_LEVEL, "max doubleDot = ", maxval(self%doubleDot))

	call LogMessage(aLog, TRACE_LOGGING_LEVEL, "size lapSurf = ", size(self%lapSurf))
	call LogMessage(aLog, TRACE_LOGGING_LEVEL, "max lapSurf = ", maxval(self%lapSurf))

	call LogMessage(aLog, TRACE_LOGGING_LEVEL, "size relVortStart = ", size(self%relVortStart))
	call LogMessage(aLog, TRACE_LOGGING_LEVEL, "max relVortStart = ", maxval(self%relVortStart))

	call LogMessage(aLog, TRACE_LOGGING_LEVEL, "size hStart = ", size(self%hStart))
	call LogMessage(aLog, TRACE_LOGGING_LEVEL, "max hStart = ", maxval(self%hStart))

	call LogMessage(aLog, TRACE_LOGGING_LEVEL, "size areaStart = ", size(self%areaStart))
	call LogMessage(aLog, TRACE_LOGGING_LEVEL, "max areaStart = ", maxval(self%areaStart))

	call EndSection(aLog)
end subroutine

!> @brief Deletes and frees memory associated with an SWESolver.
!> @param[inout] self Solver
subroutine deletePrivate(self)
	type(SWESolver), intent(inout) :: self
	if ( allocated(self%xStart)) then
		deallocate(self%xStart)
		deallocate(self%yStart)
		deallocate(self%areaStart)
		deallocate(self%divStart)
		deallocate(self%relVortStart)
		deallocate(self%mask)
		deallocate(self%u)
		deallocate(self%v)
		deallocate(self%potVort)
		deallocate(self%hStart)
		deallocate(self%doubleDot)
		deallocate(self%lapSurf)

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
		deallocate(self%relVortIn)
		deallocate(self%relVortStage1)
		deallocate(self%relVortStage2)
		deallocate(self%relVortStage3)
		deallocate(self%relVortStage4)
		deallocate(self%divIn)
		deallocate(self%divStage1)
		deallocate(self%divStage2)
		deallocate(self%divStage3)
		deallocate(self%divStage4)
		deallocate(self%areaIn)
		deallocate(self%areaStage1)
		deallocate(self%areaStage2)
		deallocate(self%areaStage3)
		deallocate(self%areaStage4)
		deallocate(self%hIn)
		deallocate(self%hStage1)
		deallocate(self%hStage2)
		deallocate(self%hStage3)
		deallocate(self%hStage4)
	endif
end subroutine

!> @brief Advances a SWE mesh forward in time by one timestep using 4th order Runge-Kutta.
!> @param[inout] self SWE solver
!> @param[inout] plane @ref PlanarSWE mesh
!> @param[in] dt timestep increment
!> @param[in] topoFn bottom topography function, must conform to the numberkindsmodule::scalarFnOf2DSpace interface
subroutine timestepPrivate( self, plane, dt, topoFn )
	type(SWESolver), intent(inout) :: self
	type(SWEMesh), intent(inout) :: plane
	real(kreal), intent(in) :: dt
	procedure(scalarFnOf2dSpace) :: topoFn
	!
	integer(kint) :: i, nP

	nP = plane%mesh%particles%N
	!
	!	RK Stage 1
	!
	do i = 1, nP
		self%xStage1(i) = dt * self%u(i)
		self%yStage1(i) = dt * self%v(i)
		self%relVortStage1 = dt * ( - (self%relVortStart(i) + plane%f0 + plane%beta * self%yStart(i)) * self%divStart(i) - &
			plane%beta * self%v(i) )
		self%divStage1 = dt * ( - self%doubleDot(i) + ( plane%f0 + plane%beta * self%yStart(i)) * self%relVortStart(i) - &
			plane%g * self%lapSurf(i) )
		self%hStage1(i) = dt * ( - self%hStart(i) * self%divStart(i))
		self%areaStage1(i) = dt * ( self%areaStart(i) * self%divStart(i))
	enddo

	!
	!	RK Stage 2
	!
	do i = 1, nP
		self%xIn(i) = self%xStart(i) + 0.5_kreal * self%xStage1(i)
		self%yIn(i) = self%yStart(i) + 0.5_kreal * self%yStage1(i)
		self%relVortIn(i) = self%relVortStart(i) + 0.5_kreal * self%relVortStage1(i)
		self%divIn(i) = self%divStart(i) + 0.5_kreal * self%divStage1(i)
		self%areaIn(i) = self%areaStart(i) + 0.5_kreal * self%areaStage1(i)
		self%hIn(i) = self%hStart(i) + 0.5_kreal * self%hStage1(i)
	enddo

	call SWEPlaneRHSIntegrals( self%u, self%v, self%doubleDot, self%lapSurf, self%xIn, self%yIn, self%relVortIn, &
				self%divIn, self%hIn, topoFn, self%areaIn, self%mask, plane%pseEps, plane%mpiParticles)

	do i = 1, nP
		self%xStage2(i) = dt * self%u(i)
		self%yStage2(i) = dt * self%v(i)
		self%relVortStage2(i) = dt * ( - (self%relVortIn(i) + plane%f0 + plane%beta * self%yIn(i)) * self%divIn(i) - &
			plane%beta * self%v(i) )
		self%divStage2(i) = dt * ( - self%doubleDot(i) + (plane%f0 + plane%beta * self%yIn(i)) * self%relVortIn(i)  - &
			plane%g * self%lapSurf(i) )
		self%hStage2(i) = dt * ( - self%hIn(i) * self%divIn(i) )
		self%areaStage2(i) = dt * ( self%areaIn(i) * self%divIn(i) )
	enddo

	!
	! RK Stage 3
	!
	do i = 1, nP
		self%xIn(i) = self%xStart(i) + 0.5_kreal * self%xStage2(i)
		self%yIn(i) = self%yStart(i) + 0.5_kreal * self%yStage2(i)
		self%relVortIn(i) = self%relVortStart(i) + 0.5_kreal * self%relVortStage2(i)
		self%divIn(i) = self%divStart(i) + 0.5_kreal * self%divStage2(i)
		self%areaIn(i) = self%areaStart(i) + 0.5_kreal * self%areaStage2(i)
		self%hIn(i) = self%hStart(i) + 0.5_kreal * self%hStage2(i)
	enddo

	call SWEPlaneRHSIntegrals( self%u, self%v, self%doubleDot, self%lapSurf, self%xIn, self%yIn, self%relVortIn, &
				self%divIn, self%hIn, topoFn, self%areaIn, self%mask, plane%pseEps, plane%mpiParticles)

	do i = 1, nP
		self%xStage3(i) = dt * self%u(i)
		self%yStage3(i) = dt * self%v(i)
		self%relVortStage3(i) = dt * ( - (self%relVortIn(i) + plane%f0 + plane%beta * self%yIn(i)) * self%divIn(i) - &
			plane%beta * self%v(i) )
		self%divStage3(i) = dt * ( - self%doubleDot(i) + (plane%f0 + plane%beta * self%yIn(i)) * self%relVortIn(i)  - &
			plane%g * self%lapSurf(i) )
		self%hStage3(i) = dt * ( - self%hIn(i) * self%divIn(i) )
		self%areaStage3(i) = dt * ( self%areaIn(i) * self%divIn(i) )
	enddo

	!
	! RK Stage 4
	!
	do i = 1, nP
		self%xIn(i) = self%xStart(i) + self%xStage3(i)
		self%yIn(i) = self%yStart(i) + self%yStage3(i)
		self%relVortIn(i) = self%relVortStart(i) + self%relVortStage3(i)
		self%divIn(i) = self%divStart(i) + self%divStage3(i)
		self%areaIn(i) = self%areaStart(i) + self%areaStage3(i)
		self%hIn(i) = self%hStart(i) + self%hStage3(i)
	enddo

	call SWEPlaneRHSIntegrals( self%u, self%v, self%doubleDot, self%lapSurf, self%xIn, self%yIn, self%relVortIn, &
				self%divIn, self%hIn, topoFn, self%areaIn, self%mask, plane%pseEps, plane%mpiParticles)

	do i = 1, nP
		self%xStage4(i) = dt * self%u(i)
		self%yStage4(i) = dt * self%v(i)
		self%relVortStage4(i) = dt * ( - (self%relVortIn(i) + plane%f0 + plane%beta * self%yIn(i)) * self%divIn(i) - &
			plane%beta * self%v(i) )
		self%divStage4(i) = dt * ( - self%doubleDot(i) + (plane%f0 + plane%beta * self%yIn(i)) * self%relVortIn(i)  - &
			plane%g * self%lapSurf(i) )
		self%hStage4(i) = dt * ( - self%hIn(i) * self%divIn(i) )
		self%areaStage4(i) = dt * ( self%areaIn(i) * self%divIn(i) )
	enddo

	!
	! RK Update
	!
	do i = 1, nP
		self%xStart(i) = self%xStart(i) + self%xStage1(i) / 6.0_kreal + self%xStage2(i) / 3.0_kreal + &
										  self%xStage3(i) / 3.0_kreal + self%xStage4(i) / 6.0_kreal
		self%yStart(i) = self%yStart(i) + self%yStage1(i) / 6.0_kreal + self%yStage2(i) / 3.0_kreal + &
										  self%yStage3(i) / 3.0_kreal + self%yStage4(i) / 6.0_kreal
		self%relVortStart(i) = self%relVortStart(i) + self%relVortStage1(i) / 6.0_kreal + self%relVortStage2(i) / 3.0_kreal + &
										  self%relVortStage3(i) / 3.0_kreal + self%relVortStage4(i) / 6.0_kreal
		self%hStart(i) = self%hStart(i) + self%hStage1(i) / 6.0_kreal + self%hStage2(i) / 3.0_kreal + &
										  self%hStage3(i) / 3.0_kreal + self%hStage4(i) / 6.0_kreal
		self%areaStart(i) = self%areaStart(i) + self%areaStage1(i) / 6.0_kreal + self%areaStage2(i) / 3.0_kreal + &
										  self%areaStage3(i) / 3.0_kreal + self%areaStage4(i) / 6.0_kreal
		self%divStart(i) = self%divStart(i) + self%divStage1(i) / 6.0_kreal + self%divStage2(i) / 3.0_kreal + &
										  self%divStage3(i) / 3.0_kreal + self%divStage4(i) / 6.0_kreal
	enddo

	call SWEPlaneRHSIntegrals( self%u, self%v, self%doubleDot, self%lapSurf, self%xStart, self%yStart, self%relVortStart, &
			self%divStart, self%hStart, topoFn, self%areaStart, self%mask, plane%pseEps, plane%mpiParticles)

	plane%mesh%particles%x(1:nP) = self%xStart
	plane%mesh%particles%y(1:nP) = self%yStart
	plane%mesh%particles%area(1:nP) = self%areaStart
	plane%relVort%scalar(1:nP) = self%relVortStart
	plane%divergence%scalar(1:nP) = self%divStart
	plane%h%scalar(1:nP) = self%hStart
	plane%velocity%xComp(1:nP) = self%u
	plane%velocity%yComp(1:nP) = self%v
	call SetBottomHeightOnMesh(plane, topoFn)
end subroutine

!
!----------------
! private methods
!----------------
!

!> @brief Computes the integral expressions on the right-hand side of the set of ODEs from a Lagrangian method of lines discretization; see @ref PlanarSWE.
!>
!> Integrals are used to compute the double dot product @f$ \nabla\vec{u}:\nabla\vec{u} @f$ and a PSE integral (see @ref PSEDirectSum)
!> approximates the Laplacian of the fluid surface.@n
!> As with the other integrals in LPM, these integrals are computed as a parallel direct summation across all MPI ranks.
!>
!> @param[out] u x-component of velocity
!> @param[out] v y-component of velocity
!> @param[out] doubleDot double dot product
!> @param[out] lapSurf Laplacian of fluid surface
!> @param[in] x x-coordinates of particles
!> @param[in] y y-coordinates of particles
!> @param[in] vort vorticity carried by particles
!> @param[in] div divergence carried by particles
!> @param[in] h fluid depth at particles
!> @param[in] topoFn bottom topography height function, must conform to numberkindsmodule::scalarFnOf2DSpace interface
!> @param[in] area area represented by each particle
!> @param[in] mask mask(i) is true if particle i is active
!> @param[in] pseEps PSE particle influence radius
!> @param[in] mpiParticles @ref MPISetup to distribute work of particle integration across all MPI ranks
subroutine SWEPlaneRHSIntegrals(u, v, doubleDot, lapSurf, x, y, vort, div, h, topoFn, area, mask, pseEps, mpiParticles)
	real(kreal), dimension(:), intent(out) :: u
	real(kreal), dimension(:), intent(out) :: v
	real(kreal), dimension(:), intent(out) :: doubleDot
	real(kreal), dimension(:), intent(out) :: lapSurf
	real(kreal), dimension(:), intent(in) :: x
	real(kreal), dimension(:), intent(in) :: y
	real(kreal), dimension(:), intent(in) :: vort
	real(kreal), dimension(:), intent(in) :: div
	real(kreal), dimension(:), intent(in) :: h
	procedure(scalarFnOf2dSpace) :: topoFn
	real(kreal), dimension(:), intent(in) :: area
	logical(klog), dimension(:), intent(in) :: mask
	real(kreal), intent(in) :: pseEps
	type(MPISetup), intent(in) :: mpiParticles
	!
	integer(kint) :: i, j, mpiErrCode
	real(kreal) :: rotStrength, potStrength, denom, ddotkernel, pseKin, lapKernel
	real(kreal) :: surfHeightI, surfHeightJ
	real(kreal) :: ux, uy, vx, vy, denom2, sqDist

	do i = mpiParticles%indexStart(procRank), mpiParticles%indexEnd(procRank)
		u(i) = 0.0_kreal
		v(i) = 0.0_kreal
		doubleDot(i) = 0.0_kreal
		lapSurf(i) = 0.0_kreal
		surfHeightI = h(i) + topoFn( x(i), y(i) )
		ux = 0.0_kreal
		uy = 0.0_kreal
		vx = 0.0_kreal
		vy = 0.0_kreal
		do j = 1, size(x)
			if ( mask(j) ) then
				surfHeightJ = h(j) + topoFn(x(j), y(j))
				sqdist = (x(i) - x(j))**2 + (y(i) - y(j))**2
				denom = 2.0_kreal * PI * sqdist
				denom2 = PI * sqdist * sqdist
				rotStrength = vort(j) * area(j) / denom
				potStrength = div(j) * area(j) / denom

				pseKin = sqrt( sqDist ) / pseEps
				lapKernel = bivariateLaplacianKernel8( pseKin ) / pseEps**2
				lapSurf(i) = lapSurf(i) + lapKernel * (surfHeightJ - surfHeightI) * area(j)

				if ( i == j ) cycle

				u(i) = u(i) - (y(i) - y(j)) * rotStrength + (x(i) - x(j)) * potStrength
				v(i) = v(i) + (x(i) - x(j)) * rotStrength + (y(i) - y(j)) * potStrength

				ux = ux + potStrength - ( (x(i) - x(j)) * ( (x(i) - x(j)) * div(j) - (y(i) - y(j)) * vort(j) ) ) * &
										 area(j) / denom2
				uy = uy - rotStrength - ( (y(i) - y(j)) * ( (x(i) - x(j)) * div(j) - (y(i) - y(j)) * vort(j) ) ) * &
										 area(j) / denom2
				vx = vx + rotStrength - ( (x(i) - x(j)) * ( (y(i) - y(j)) * div(j) + (x(i) - x(j)) * vort(j) ) ) * &
										 area(j) / denom2
				vy = vy + potStrength - ( (y(i) - y(j)) * ( (y(i) - y(j)) * div(j) + (x(i) - x(j)) * vort(j) ) ) * &
										 area(j) / denom2
			endif
		enddo
!		do j = i + 1, size(x)
!			if ( mask(j) ) then
!				surfHeightJ = h(j) + topoFn(x(j), y(j))
!				sqdist = (x(i) - x(j))**2 + (y(i) - y(j))**2
!				denom = 2.0_kreal * PI * sqdist
!				denom2 = PI * sqdist * sqdist
!				rotStrength = vort(j) * area(j) / denom
!				potStrength = div(j) * area(j) / denom
!
!				u(i) = u(i) - (y(i) - y(j)) * rotStrength + (x(i) - x(j)) * potStrength
!				v(i) = v(i) + (x(i) - x(j)) * rotStrength + (y(i) - y(j)) * potStrength
!
!				ux = ux + potStrength - ( (x(i) - x(j)) * ( (x(i) - x(j)) * div(j) - (y(i) - y(j)) * vort(j) ) ) * &
!					 			         area(j) / denom2
!				uy = uy - rotStrength - ( (y(i) - y(j)) * ( (x(i) - x(j)) * div(j) - (y(i) - y(j)) * vort(j) ) ) * &
!										 area(j) / denom2
!				vx = vx + rotStrength - ( (x(i) - x(j)) * ( (y(i) - y(j)) * div(j) + (x(i) - x(j)) * vort(j) ) ) * &
!										 area(j) / denom2
!				vy = vy + potStrength - ( (y(i) - y(j)) * ( (y(i) - y(j)) * div(j) + (x(i) - x(j)) * vort(j) ) ) * &
!										 area(j) / denom2
!
!				pseKin = sqrt( sqDist ) / pseEps
!				lapKernel = bivariateLaplacianKernel8( pseKin ) / pseEps**2
!				lapSurf(i) = lapSurf(i) + lapKernel * (surfHeightJ - surfHeightI) * area(j)
!			endif
!		enddo
		lapSurf(i) = lapSurf(i) / pseEps**2
		doubleDot(i) = ux * ux + 2.0_kreal * uy * vx + vy * vy
	enddo

	do i = 0, numProcs - 1
		call MPI_BCAST(u(mpiParticles%indexStart(i):mpiParticles%indexEnd(i)), mpiParticles%messageLength(i), &
				MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)
		call MPI_BCAST(v(mpiParticles%indexStart(i):mpiParticles%indexEnd(i)), mpiParticles%messageLength(i), &
				MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)
		call MPI_BCAST(doubleDot(mpiParticles%indexStart(i):mpiParticles%indexEnd(i)), mpiParticles%messageLength(i), &
				MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)
		call MPI_BCAST(lapSurf(mpiParticles%indexStart(i):mpiParticles%indexEnd(i)), mpiParticles%messageLength(i), &
				MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)
	enddo

!	call LogMessage(log, DEBUG_LOGGING_LEVEL, trim(logKey)//" max u = ", maxval(u) )
!	call LogMessage(log, DEBUG_LOGGING_LEVEL, trim(logKey)//" min u = ", minval(u) )
!	call LogMessage(log, DEBUG_LOGGING_LEVEL, trim(logKey)//" max v = ", maxval(v) )
!	call LogMessage(log, DEBUG_LOGGING_LEVEL, trim(logKey)//" min v = ", minval(v) )
!	call LogMessage(log, DEBUG_LOGGING_LEVEL, trim(logKey)//" max(ddot) = ", maxval(doubleDot) )
!	call LogMessage(log, DEBUG_LOGGING_LEVEL, trim(logKey)//" min(ddot) = ", minval(doubleDot) )
!	call LogMessage(log, DEBUG_LOGGING_LEVEL, trim(logKey)//" max(lapS) = ", maxval(lapSurf) )
!	call LogMessage(log, DEBUG_LOGGING_LEVEL, trim(logKey)//" min(lapS) = ", minval(lapSurf) )
end subroutine

!> @brief Initializes a logger for the planar SWE solver module
!>
!> Output is controlled both by message priority and by MPI Rank
!> @param aLog Target Logger object
!> @param rank Rank of this processor
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

!> @}
end module
