module BetaPlaneSolverModule
!> @file BetaPlaneSolver.f90
!> Data structure for integrating a planar incompressible flow with a background rotation forward in time using a fourth-order Runge-Kutta method.
!> @author Peter Bosler, Sandia National Laboratories Center for Computing Research
!> 
!>
!> @defgroup BetaPlaneSolver BetaPlaneSolver
!> Data structure for integrating a planar incompressible flow with rotation forward in time using a fourth-order Runge-Kutta method.
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
use BetaPlaneMeshModule

implicit none

include 'mpif.h'

private
!
!----------------
! Module types and public declarations
!----------------
!
public BetaPlaneSolver, New, Delete
public Timestep

type BetaPlaneSolver
	real(kreal), pointer, dimension(:) :: xIn => null()
	real(kreal), pointer, dimension(:) :: xStage1 => null()
	real(kreal), pointer, dimension(:) :: xStage2 => null()
	real(kreal), pointer, dimension(:) :: xStage3 => null()
	real(kreal), pointer, dimension(:) :: xStage4 => null()
	real(kreal), pointer, dimension(:) :: yIn => null()
	real(kreal), pointer, dimension(:) :: yStage1 => null()
	real(kreal), pointer, dimension(:) :: yStage2 => null()
	real(kreal), pointer, dimension(:) :: yStage3 => null()
	real(kreal), pointer, dimension(:) :: yStage4 => null()
	real(kreal), pointer, dimension(:) :: relVortIn => null()
	real(kreal), pointer, dimension(:) :: relVortStage1 => null()
	real(kreal), pointer, dimension(:) :: relVortStage2 => null()
	real(kreal), pointer, dimension(:) :: relVortStage3 => null()
	real(kreal), pointer, dimension(:) :: relVortStage4 => null()
	
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
character(len=28), save :: logKey = 'BVESolver'
integer(kint), parameter :: logLevel = DEBUG_LOGGING_LEVEL
character(len=MAX_STRING_LENGTH) :: logString

contains

!
!----------------
! public methods
!----------------
!
subroutine newPrivate( self, betaPlane )
	type(BetaPlaneSolver), intent(out) :: self
	type(BetaPlaneMesh), intent(in) :: betaPlane
	!
	integer(kint) :: nP
	
	if ( .NOT. logInit ) call InitLogger(log, procRank)
	
	nP = betaPlane%mesh%particles%N
	
	allocate(self%xIn(nP))
	allocate(self%xStage1(nP))
	allocate(self%xStage2(nP))
	allocate(self%xStage3(nP))
	allocate(self%xStage4(nP))
	allocate(self%yIn(nP))
	allocate(self%yStage1(nP))
	allocate(self%yStage2(nP))
	allocate(self%yStage3(nP))
	allocate(self%yStage4(nP))
	allocate(self%relVortIn(nP))
	allocate(self%relVortStage1(nP))
	allocate(self%relVortStage2(nP))
	allocate(self%relVortStage3(nP))
	allocate(self%relVortStage4(nP))
end subroutine

subroutine deletePrivate( self ) 
	type(BetaPlaneSolver), intent(inout) :: self

	if ( associated(self%xIn) ) then
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
	endif
end subroutine

subroutine timestepPrivate( self, betaPlane, dt )
	type(BetaPlaneSolver), intent(inout) :: self
	type(BetaPlaneMesh), intent(inout) :: betaPlane
	real(kreal), intent(in) :: dt
	!
	integer(kint) :: i, nP
	
	nP = betaPlane%mesh%particles%N
	!
	! RK Stage 1
	!
	self%xStage1 = dt * betaPlane%velocity%xComp(1:nP)
	self%yStage1 = dt * betaPlane%velocity%yComp(1:nP)
	self%relVortStage1 = - dt * betaPlane%beta * betaPlane%velocity%yComp(1:nP)
	
	!
	! RK Stage 2
	!
	self%xIn = betaPlane%mesh%particles%x(1:NP) + 0.5_kreal * self%xStage1
	self%yIn = betaPlane%mesh%particles%y(1:nP) + 0.5_kreal * self%yStage1
	self%relVortIn = betaPlane%relVort%scalar(1:nP) + 0.5_kreal * self%relVortStage1
	
	!call ResetPeriodicBounds(self%xIn)
	
	call BetaPlaneVelocity( self%xStage2, self%yStage2, self%xIn, self%yIn, self%relVortIn, &
		betaPlane%mesh%particles%area, betaPlane%mesh%particles%isActive, betaPlane%mpiParticles)
		
	self%relVortStage2 = - dt * betaPlane%beta * self%yStage2
	self%xStage2 = dt * self%xStage2
	self%yStage2 = dt * self%yStage2
	
	!
	! RK Stage 3
	!
	self%xIn = betaPlane%mesh%particles%x(1:NP) + 0.5_kreal * self%xStage2
	self%yIn = betaPlane%mesh%particles%y(1:nP) + 0.5_kreal * self%yStage2
	self%relVortIn = betaPlane%relVort%scalar(1:nP) + 0.5_kreal * self%relVortStage2
	
	!call ResetPeriodicBounds(self%xIn)
	
	call BetaPlaneVelocity( self%xStage3, self%yStage3, self%xIn, self%yIn, self%relVortIn, &
		betaPlane%mesh%particles%area, betaPlane%mesh%particles%isActive, betaPlane%mpiParticles)
		
	self%relVortStage3 = - dt * betaPlane%beta * self%yStage3
	self%xStage3 = dt * self%xStage3
	self%yStage3 = dt * self%yStage3
	
	!
	! RK Stage 4
	!
	self%xIn = betaPlane%mesh%particles%x(1:NP) + self%xStage3
	self%yIn = betaPlane%mesh%particles%y(1:nP) + self%yStage3
	self%relVortIn = betaPlane%relVort%scalar(1:nP) + self%relVortStage3
	
	!call ResetPeriodicBounds(self%xIn)
	
	call BetaPlaneVelocity( self%xStage4, self%yStage4, self%xIn, self%yIn, self%relVortIn, &
		betaPlane%mesh%particles%area, betaPlane%mesh%particles%isActive, betaPlane%mpiParticles)
		
	self%relVortStage4 = - dt * betaPlane%beta * self%yStage4
	self%xStage4 = dt * self%xStage4
	self%yStage4 = dt * self%yStage4
	
	!
	! RK Update
	!
	betaPlane%mesh%particles%x(1:nP) = betaPlane%mesh%particles%x(1:nP) + self%xStage1 / 6.0_kreal + &
		self%xStage2 / 3.0_kreal + self%xStage3 / 3.0_kreal + self%xStage4 / 6.0_kreal
	betaPlane%mesh%particles%y(1:nP) = betaPlane%mesh%particles%y(1:nP) + self%yStage1 / 6.0_kreal + &
		self%yStage2 / 3.0_kreal + self%yStage3 / 3.0_kreal + self%yStage4 / 6.0_kreal
	betaPlane%relVort%scalar(1:nP) = betaPlane%relVort%scalar(1:nP) + self%relVortStage1 / 6.0_kreal + &
		self%relVortStage2 / 3.0_kreal + self%relVortStage3 / 3.0_kreal + self%relVortStage4 / 6.0_kreal
	
	!call ResetPeriodicBounds( betaPlane )
	
	call SetVelocityOnMesh(betaPlane)
	call SetStreamFunctionsOnMesh(betaPlane)
end subroutine


!
!----------------
! private methods
!----------------
!
subroutine BetaPlaneVelocity( u, v, x, y, relVortIn, areaIn, activemask, mpiParticles )
	real(kreal), dimension(:), intent(out) :: u
	real(kreal), dimension(:), intent(out) :: v
	real(kreal), dimension(:), intent(in) :: x
	real(kreal), dimension(:), intent(in) :: y
	real(kreal), dimension(:), intent(in) :: relVortIn
	real(kreal), dimension(:), intent(in) :: areaIn
	logical(klog), dimension(:), intent(in) :: activemask
	type(MPISetup), intent(in) :: mpiParticles
	!
	integer(kint) :: i, j, mpiErrCode
	real(kreal) :: strength
	
	do i = mpiParticles%indexStart(procRank), mpiParticles%indexEnd(procRank)
		u(i) = 0.0_kreal
		v(i) = 0.0_kreal
		do j = 1, i - 1
			if ( activeMask(j) ) then
				strength = 0.5_kreal * relVortIn(j) * areaIn(j) / ( cosh(2.0_kreal*PI*(y(i) - y(j))) - &
					cos(2.0_kreal * PI * (x(i) - x(j))))
				u(i) = u(i) - sinh( 2.0_kreal * PI * (y(i) - y(j))) * strength
				v(i) = v(i) + sin(2.0_kreal * PI * (x(i) - x(j))) * strength
			endif
		enddo
		do j = i + 1, size(x)
			if ( activeMask(j) ) then
				strength = 0.5_kreal * relVortIn(j) * areaIn(j) / ( cosh(2.0_kreal*PI*(y(i) - y(j))) - &
					cos(2.0_kreal * PI * (x(i) - x(j))))
				u(i) = u(i) - sinh( 2.0_kreal * PI * (y(i) - y(j))) * strength
				v(i) = v(i) + sin(2.0_kreal * PI * (x(i) - x(j))) * strength
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