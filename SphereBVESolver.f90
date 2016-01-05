module SphereBVESolverModule

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
use SphereBVEModule

implicit none

include 'mpif.h'

private
!
!----------------
! Module types and public declarations
!----------------
!
public BVESolver, New, Delete
public Timestep

type BVESolver
	real(kreal), allocatable :: xStart(:)
	real(kreal), allocatable :: yStart(:)
	real(kreal), allocatable :: zStart(:)
	real(kreal), allocatable :: relVortStart(:)
	real(kreal), allocatable :: area(:)
	real(kreal), allocatable :: u(:)
	real(kreal), allocatable :: v(:)
	real(kreal), allocatable :: w(:)
	logical(klog), allocatable :: mask(:)
	
	real(kreal), allocatable :: xIn(:) 
	real(kreal), allocatable :: xStage1(:) 
	real(kreal), allocatable :: xStage2(:) 
	real(kreal), allocatable :: xStage3(:) 
	real(kreal), allocatable :: xStage4(:) 
	real(kreal), allocatable :: yIn(:) 
	real(kreal), allocatable :: yStage1(:) 
	real(kreal), allocatable :: yStage2(:) 
	real(kreal), allocatable :: yStage3(:) 
	real(kreal), allocatable :: yStage4(:) 
	real(kreal), allocatable :: zIn(:) 
	real(kreal), allocatable :: zStage1(:) 
	real(kreal), allocatable :: zStage2(:) 
	real(kreal), allocatable :: zStage3(:) 
	real(kreal), allocatable :: zStage4(:) 
	real(kreal), allocatable :: relVortIn(:) 
	real(kreal), allocatable :: relVortStage1(:) 
	real(kreal), allocatable :: relVortStage2(:) 
	real(kreal), allocatable :: relVortStage3(:) 
	real(kreal), allocatable :: relVortStage4(:) 
	
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

contains
!
!----------------
! public methods
!----------------
!
subroutine newPrivate(self, sphereBVE )
	type(BVESolver), intent(out) :: self
	type(BVEMesh), intent(in) :: sphereBVE
	integer(kint) :: nP
	
	if (.NOT. logInit ) call InitLogger(log, procRank)
	
	nP = sphereBVE%mesh%particles%N
	
	allocate(self%xStart(nP))
	allocate(self%yStart(nP))
	allocate(self%zStart(nP))
	allocate(self%u(nP))
	allocate(self%v(nP))
	allocate(self%w(nP))
	allocate(self%mask(nP))
	allocate(self%area(nP))
	allocate(self%relVortStart(nP))
	
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
	allocate(self%zIn(nP))
	allocate(self%zStage1(nP))
	allocate(self%zStage2(nP))
	allocate(self%zStage3(nP))
	allocate(self%zStage4(nP))
	allocate(self%relVortIn(nP))
	allocate(self%relVortStage1(nP))
	allocate(self%relVortStage2(nP))
	allocate(self%relVortStage3(nP))
	allocate(self%relVortStage4(nP))
	
	self%xStart(1:nP) = sphereBVE%mesh%particles%x(1:nP)
	self%yStart(1:nP) = sphereBVE%mesh%particles%y(1:nP)
	self%zStart(1:nP) = sphereBVE%mesh%particles%z(1:nP)
	self%u(1:nP) = sphereBVE%velocity%xComp(1:nP)
	self%v(1:nP) = sphereBVE%velocity%yComp(1:nP)
	self%w(1:nP) = sphereBVE%velocity%zComp(1:nP)
	self%relVortStart(1:nP) = sphereBVE%relVort%scalar(1:nP)
	self%area(1:nP) = sphereBVE%mesh%particles%area(1:nP)
	self%mask(1:nP) = sphereBVE%mesh%particles%isActive(1:nP)
	!$acc enter data declare create(self, sphereBVE%radius, sphereBVE%rotationRate) &
	!$acc copyin( self%xStart, self%yStart, self%zStart, self%relVortStart, self%mask, self%area) &
	!$acc copyin( self%u, self%v, self%w) &
	!$acc create( self%xIn, self%xStage1, self%xStage2, self%xStage3, self%xStage4) &
	!$acc create( self%yIn, self%yStage1, self%yStage2, self%yStage3, self%yStage4) &
	!$acc create( self%zIn, self%zStage1, self%zStage2, self%zStage3, self%zStage4) &
	!$acc create( self%relVortIn, self%relVortStage1, self%relVortStage2, self%relVortStage3, self%relVortStage4)
end subroutine

subroutine deletePrivate(self)
	type(BVESolver), intent(inout) :: self
	if ( allocated(self%xIn) ) then
		deallocate(self%xStart)
		deallocate(self%yStart)
		deallocate(self%zStart)
		deallocate(self%relVortStart)
		deallocate(self%u)
		deallocate(self%v)
		deallocate(self%w)
	
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
		deallocate(self%zIn)
		deallocate(self%zStage1)
		deallocate(self%zStage2)
		deallocate(self%zStage3)
		deallocate(self%zStage4)
		deallocate(self%relVortIn)
		deallocate(self%relVortStage1)
		deallocate(self%relVortStage2)
		deallocate(self%relVortStage3)
		deallocate(self%relVortStage4)
		
		!$acc exit data &
		!$acc delete(self%xStart, self%yStart, self%zStart, self%u, self%v, self%w) &
		!$acc delete(self%mask, self%area, sphereBVE%radius, sphereBVE%rotationRate) &
		!$acc delete(self%xIn, self%xStage1, self%xStage2, self%xStage3, self%xStage4) &
		!$acc delete(self%yIn, self%yStage1, self%yStage2, self%yStage3, self%yStage4) &
		!$acc delete(self%zIn, self%zStage1, self%zStage2, self%zStage3, self%zStage4) &
		!$acc delete(self%relVortIn, self%relVortStage1, self%relVortStage2, self%relVortStage3, self%relVortStage4) &
		!$acc delete(self)
	endif
end subroutine

subroutine timestepPrivate( self, sphereBVE, dt )
	type(BVESolver), intent(inout) :: self
	type(BVEMesh), intent(inout) :: sphereBVE
	real(kreal), intent(in) :: dt
	!
	integer(kint) :: i, nP
	!$acc routine(BVESphereVelocity) gang

	nP = sphereBVE%mesh%particles%N
	!
	!	RK Stage 1
	!
	!$acc parallel loop
	do i = 1, nP
		self%xStage1(i) = dt * self%u(i)
		self%yStage1(i) = dt * self%v(i)
		self%zStage1(i) = dt * self%w(i)
		self%relVortStage1(i) =  - dt * 2.0_kreal * sphereBVE%rotationRate * self%w(i) / sphereBVE%radius
	enddo
	
	!
	! RK Stage 2
	!
	!$acc parallel loop
	do i = 1, nP
		self%xIn(i) = self%xStart(i) + 0.5_kreal * self%xStage1(i)
		self%yIn(i) = self%yStart(i) + 0.5_kreal * self%yStage1(i)
		self%zIn(i) = self%zStart(i) + 0.5_kreal * self%zStage1(i)
		self%relVortIn(i) = self%relVortStart(i) + 0.5_kreal * self%relVortStage1(i)
	enddo

	call BVESphereVelocity(self%xStage2, self%yStage2, self%zStage2, self%xIn, self%yIn, self%zIn, &
			self%relVortIn, self%area, sphereBVE%radius, sphereBVE%rotationRate, &
			self%mask, sphereBVE%mpiParticles)
	
	!$acc parallel loop
	do i = 1, nP
		self%relVortStage2(i) = - dt * 2.0_kreal * sphereBVE%rotationRate * self%zStage2(i) / sphereBVE%radius
		self%xStage2(i) = dt * self%xStage2(i)
		self%yStage2(i) = dt * self%yStage2(i)
		self%zStage2(i) = dt * self%zStage2(i)
	enddo
	
	!
	! RK Stage 3
	!
	!$acc parallel loop
	do i = 1, nP
		self%xIn(i) = self%xStart(i) + 0.5_kreal * self%xStage2(i)
		self%yIn(i) = self%yStart(i) + 0.5_kreal * self%yStage2(i)
		self%zIn(i) = self%zStart(i) + 0.5_kreal * self%zStage2(i)
		self%relVortIn(i) = self%relVortStart(i) + 0.5_kreal * self%relVortStage2(i)
	enddo
	
	call BVESphereVelocity( self%xStage3, self%yStage3, self%zStage3, self%xIn, self%yIn, self%zIn, &
			self%relVortIn, self%area, sphereBVE%radius, sphereBVE%rotationRate, &
			self%mask, sphereBVE%mpiParticles)
	
	!$acc parallel loop
	do i = 1, nP
		self%relVortStage3(i) = - dt * 2.0_kreal * sphereBVE%rotationRate * self%zStage3(i) / sphereBVE%radius
		self%xStage3(i) = dt * self%xStage3(i)
		self%yStage3(i) = dt * self%yStage3(i)
		self%zStage3(i) = dt * self%zStage3(i)
	enddo
	
	!
	! RK Stage 4
	!
	!$acc parallel loop
	do i = 1, nP
		self%xIn(i) = self%xStart(i) + self%xStage3(i)
		self%yIn(i) = self%yStart(i) + self%yStage3(i)
		self%zIn(i) = self%zStart(i) + self%zStage3(i)
		self%relVortIn(i) = self%relVortStart(i) + self%relVortStage3(i)
	enddo
	
	call BVESphereVelocity( self%xStage4, self%yStage4, self%zStage4, self%xIn, self%yIn, self%zIn, &
			self%relVortIn, self%area, sphereBVE%radius, sphereBVE%rotationRate, &
			self%mask, sphereBVE%mpiParticles)
	
	!$acc parallel loop
	do i = 1, nP
		self%relVortStage4(i) = - dt * 2.0_kreal * sphereBVE%rotationRate * self%zStage4(i) / sphereBVE%radius
		self%xStage4(i) = dt * self%xStage4(i)
		self%yStage4(i) = dt * self%yStage4(i)
		self%zStage4(i) = dt * self%zStage4(i)
	enddo
	
	!
	!	rk4 update
	!	
	
	!$acc parallel loop
	do i = 1, nP
		self%xStart(i) = self%xStart(i) + self%xStage1(i) / 6.0_kreal + self%xStage2(i) / 3.0_kreal + &
			self%xStage3(i) / 3.0_kreal + self%xStage4(i) / 6.0_kreal
		self%yStart(i) = self%yStart(i) + self%yStage1(i) / 6.0_kreal + self%yStage2(i) / 3.0_kreal + &
			self%yStage3(i) / 3.0_kreal + self%yStage4(i) / 6.0_kreal
		self%zStart(i) = self%zStart(i) + self%zStage1(i) / 6.0_kreal + self%zStage2(i) / 3.0_kreal + &
			self%zStage3(i) / 3.0_kreal + self%zStage4(i) / 6.0_kreal
		self%relVortStart(i) = self%relVortStart(i) + self%relVortStage1(i) / 6.0_kreal + &
			self%relVortStage2(i) / 3.0_kreal + self%relVortStage3(i) / 3.0_kreal + self%relVortStage4(i) / 6.0_kreal
	enddo
	
	!$acc update host(self%xStart, self%yStart, self%zStart, self%relVortStart)
	
	sphereBVE%mesh%particles%x(1:nP) = self%xStart
	sphereBVE%mesh%particles%y(1:nP) = self%yStart
	sphereBVE%mesh%particles%z(1:nP) = self%zStart
	sphereBVE%relVort%scalar(1:nP) = self%relVortStart
	
	!call SetVelocityOnMesh(sphereBVE)	
	
	!self%u = sphereBVE%velocity%xComp(1:nP)
	!self%v = sphereBVE%velocity%yComp(1:nP)
	!self%w = sphereBVE%velocity%zComp(1:nP)
	!! !$acc update device (self%u, self%v, self%w)
	
	call BVESphereVelocity(self%u, self%v, self%w, self%xStart, self%yStart, self%zStart, self%relVortStart, &
							self%area, sphereBVE%radius, sphereBVE%rotationRate, self%mask, sphereBVE%mpiParticles)
	
	sphereBVE%velocity%xComp(1:nP) = self%u
	sphereBVE%velocity%yComp(1:nP) = self%v
	sphereBVE%velocity%zComp(1:nP) = self%w
	
	call SetStreamFunctionsOnMesh(sphereBVE)
end subroutine


!
!----------------
! private methods
!----------------
!
subroutine BVESphereVelocity( u, v, w, x, y, z, relVortIn, areaIn, sphereRadius, rotationRate, activeMask, mpiParticles )
	!$acc routine gang
	real(kreal), dimension(:), intent(out) :: u
	real(kreal), dimension(:), intent(out) :: v
	real(kreal), dimension(:), intent(out) :: w
	real(kreal), dimension(:), intent(in) :: x
	real(kreal), dimension(:), intent(in) :: y
	real(kreal), dimension(:), intent(in) :: z
	real(kreal), dimension(:), intent(in) :: relVortIn
	real(kreal), dimension(:), intent(in) :: areaIn
	real(kreal), intent(in) :: sphereRadius
	real(kreal), intent(in) :: rotationRate
	logical(klog), dimension(:), intent(in) :: activeMask
	type(MPISetup), intent(in) :: mpiParticles
	!
	integer(kint) :: i, j, mpiErrCode
	real(kreal) :: strength
	
	!$acc parallel loop gang
	do i = mpiParticles%indexStart(procRank), mpiParticles%indexEnd(procRank)
		u(i) = 0.0_kreal
		v(i) = 0.0_kreal
		w(i) = 0.0_kreal
		!$acc loop vector
		do j = 1, i - 1
			if ( activeMask(j) ) then
				strength = -relVortIn(j)*areaIn(j) / ( 4.0_kreal * PI * sphereRadius * &
					( sphereRadius * sphereRadius - x(i)*x(j) - y(i)*y(j) - z(i)*z(j) ))
				u(i) = u(i) + ( y(i)*z(j) - z(i)*y(j) ) * strength
				v(i) = v(i) + ( z(i)*x(j) - x(i)*z(j) ) * strength
				w(i) = w(i) + ( x(i)*y(j) - y(i)*x(j) ) * strength
			endif
		enddo
		!$acc loop vector
		do j = i + 1, size(x)
			if ( activeMask(j) ) then
				strength = -relVortIn(j)*areaIn(j) / ( 4.0_kreal * PI * sphereRadius * &
					( sphereRadius * sphereRadius - x(i)*x(j) - y(i)*y(j) - z(i)*z(j) ))
				u(i) = u(i) + ( y(i)*z(j) - z(i)*y(j) ) * strength
				v(i) = v(i) + ( z(i)*x(j) - x(i)*z(j) ) * strength
				w(i) = w(i) + ( x(i)*y(j) - y(i)*x(j) ) * strength
			endif
		enddo
	enddo
	
	do i = 0, numProcs - 1
		call MPI_BCAST( u(mpiParticles%indexStart(i):mpiParticles%indexEnd(i)), mpiParticles%messageLength(i), &
				MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)
		call MPI_BCAST( v(mpiParticles%indexStart(i):mpiParticles%indexEnd(i)), mpiParticles%messageLength(i), &
				MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)
		call MPI_BCAST( w(mpiParticles%indexStart(i):mpiParticles%indexEnd(i)), mpiParticles%messageLength(i), &
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

end module