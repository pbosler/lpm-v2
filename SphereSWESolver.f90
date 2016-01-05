module SphereSWESolverModule

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
use SphereSWEModule

implicit none

include 'mpif.h'

private
!
!----------------
! Module types and public declarations
!----------------
!
public SWESolver, New, Delete
public Timestep

type SWESolver
	real(kreal), allocatable, dimension(:) :: xStart
	real(kreal), allocatable, dimension(:) :: yStart
	real(kreal), allocatable, dimension(:) :: zStart
	real(kreal), allocatable, dimension(:) :: relVortStart
	real(kreal), allocatable, dimension(:) :: divStart
	real(kreal), allocatable, dimension(:) :: hStart
	real(kreal), allocatable, dimension(:) :: areaStart
	real(kreal), allocatable, dimension(:) :: u
	real(kreal), allocatable, dimension(:) :: v
	real(kreal), allocatable, dimension(:) :: w
	real(kreal), allocatable, dimension(:) :: potVort
	real(kreal), allocatable, dimension(:) :: doubleDot
	real(kreal), allocatable, dimension(:) :: lapSurf
	logical(klog), allocatable, dimension(:) :: mask
	
	real(kreal), allocatable, dimension(:) :: xIn
	real(kreal), allocatable, dimension(:) :: xStage1
	real(kreal), allocatable, dimension(:) :: xStage2
	real(kreal), allocatable, dimension(:) :: xStage3
	real(kreal), allocatable, dimension(:) :: xStage4
	real(kreal), allocatable, dimension(:) :: yIn
	real(kreal), allocatable, dimension(:) :: yStage1
	real(kreal), allocatable, dimension(:) :: yStage2
	real(kreal), allocatable, dimension(:) :: yStage3
	real(kreal), allocatable, dimension(:) :: yStage4
	real(kreal), allocatable, dimension(:) :: zIn
	real(kreal), allocatable, dimension(:) :: zStage1
	real(kreal), allocatable, dimension(:) :: zStage2
	real(kreal), allocatable, dimension(:) :: zStage3
	real(kreal), allocatable, dimension(:) :: zStage4
	real(kreal), allocatable, dimension(:) :: relVortIn
	real(kreal), allocatable, dimension(:) :: relVortStage1
	real(kreal), allocatable, dimension(:) :: relVortStage2
	real(kreal), allocatable, dimension(:) :: relVortStage3
	real(kreal), allocatable, dimension(:) :: relVortStage4
	real(kreal), allocatable, dimension(:) :: divIn
	real(kreal), allocatable, dimension(:) :: divStage1
	real(kreal), allocatable, dimension(:) :: divStage2
	real(kreal), allocatable, dimension(:) :: divStage3
	real(kreal), allocatable, dimension(:) :: divStage4
	real(kreal), allocatable, dimension(:) :: hIn
	real(kreal), allocatable, dimension(:) :: hStage1
	real(kreal), allocatable, dimension(:) :: hStage2
	real(kreal), allocatable, dimension(:) :: hStage3
	real(kreal), allocatable, dimension(:) :: hStage4
	real(kreal), allocatable, dimension(:) :: areaIn
	real(kreal), allocatable, dimension(:) :: areaStage1
	real(kreal), allocatable, dimension(:) :: areaStage2
	real(kreal), allocatable, dimension(:) :: areaStage3
	real(kreal), allocatable, dimension(:) :: areaStage4
	
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
character(len=28), save :: logKey = 'SWESphereK4'
integer(kint), parameter :: logLevel = DEBUG_LOGGING_LEVEL
character(len=MAX_STRING_LENGTH) :: logString

contains
!
!----------------
! public methods
!----------------
!
subroutine newPrivate( self, sphere, topoFn)
	type(SWESolver), intent(out) :: self
	type(SWEMesh), intent(inout) :: sphere
	procedure(scalarFnOf3dSpace) :: topoFn
	!
	integer(kint) :: i, nP
	
	if (.NOT. logInit ) call InitLogger(log, procRank)
	
	nP = sphere%mesh%particles%N

	allocate(self%xStart(1:nP))
	allocate(self%yStart(1:nP))
	allocate(self%zStart(1:nP))
	allocate(self%relVortStart(1:nP))
	allocate(self%divStart(1:nP))
	allocate(self%hStart(1:nP))
	allocate(self%areaStart(1:nP))
	allocate(self%u(1:nP))
	allocate(self%v(1:nP))
	allocate(self%w(1:nP))
	allocate(self%potVort(1:nP))
	allocate(self%doubleDot(1:nP))
	allocate(self%lapSurf(1:nP))
	allocate(self%mask(1:nP))
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
	allocate(self%zIn(1:nP))
	allocate(self%zStage1(1:nP))
	allocate(self%zStage2(1:nP))
	allocate(self%zStage3(1:nP))
	allocate(self%zStage4(1:nP))
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
	allocate(self%hIn(1:nP))
	allocate(self%hStage1(1:nP))
	allocate(self%hStage2(1:nP))
	allocate(self%hStage3(1:nP))
	allocate(self%hStage4(1:nP))
	allocate(self%areaIn(1:nP))
	allocate(self%areaStage1(1:nP))
	allocate(self%areaStage2(1:nP))
	allocate(self%areaStage3(1:nP))
	allocate(self%areaStage4(1:nP))
	
	self%xStart = sphere%mesh%particles%x(1:nP)
	self%yStart = sphere%mesh%particles%y(1:nP)
	self%zStart = sphere%mesh%particles%z(1:nP)
	self%areaStart = sphere%mesh%particles%area(1:nP)
	self%relVortStart = sphere%relVort%scalar(1:nP)
	self%divStart = sphere%divergence%scalar(1:nP)
	self%hStart = sphere%h%scalar(1:nP)
	self%mask = sphere%mesh%particles%isActive(1:nP)
	
	call SWESphereRHSIntegrals( self%u, self%v, self%w, self%doubleDot, self%lapSurf, self%xStart, self%yStart, self%zStart, &
				self%relVortStart, self%divStart, self%hStart, self%areaStart, topoFn, self%mask, sphere%radius, sphere%pseEps, &
				sphere%mpiParticles)
				
	! DEBUG
	if ( associated(sphere%tracers) .AND. sphere%tracers(1)%nDim == 1 ) then
		sphere%tracers(1)%N = nP
		do i = 1, nP
			sphere%velocity%xComp(i) = self%u(i)
			sphere%velocity%yComp(i) = self%v(i)
			sphere%velocity%zComp(i) = self%w(i)
			sphere%tracers(1)%scalar(i) = self%doubleDot(i)
		enddo
	endif
	
end subroutine

subroutine deletePrivate(self)
	type(SWESolver), intent(inout) :: self
	
	if ( allocated(self%xStart)) then
		deallocate(self%xStart)
		deallocate(self%yStart)
		deallocate(self%zStart)
		deallocate(self%relVortStart)
		deallocate(self%divStart)
		deallocate(self%hStart)
		deallocate(self%areaStart)
		deallocate(self%u)
		deallocate(self%v)
		deallocate(self%w)
		deallocate(self%potVort)
		deallocate(self%doubleDot)
		deallocate(self%lapSurf)
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
		deallocate(self%divIn)
		deallocate(self%divStage1)
		deallocate(self%divStage2)
		deallocate(self%divStage3)
		deallocate(self%divStage4)
		deallocate(self%hIn)
		deallocate(self%hStage1)
		deallocate(self%hStage2)
		deallocate(self%hStage3)
		deallocate(self%hStage4)
		deallocate(self%areaIn)
		deallocate(self%areaStage1)
		deallocate(self%areaStage2)
		deallocate(self%areaStage3)
		deallocate(self%areaStage4)
	endif
end subroutine

subroutine timestepPrivate( self, sphere, dt, topoFn)
	type(SWESolver), intent(inout) :: self
	type(SWEMesh), intent(inout) :: sphere
	real(kreal), intent(in) :: dt
	procedure(scalarFnOf3dSpace) :: topoFn
	!
	integer(kint) :: i, nP
	real(kreal) :: omg
	
	nP = sphere%mesh%particles%N
	omg = sphere%rotationRate
	
	!
	!	RK Stage 1
	!
	do i = 1, nP
		self%xStage1(i) = dt * self%u(i)
		self%yStage1(i) = dt * self%v(i)
		self%zStage1(i) = dt * self%w(i)
		self%relVortStage1(i) = dt * ( - ( self%relVortStart(i) + 2.0_kreal * omg * self%zStart(i) / &
			sphere%radius ) * self%divStart(i) - 2.0_kreal * omg * self%w(i) / sphere%radius )
		self%divStage1(i) = dt * ( - self%doubleDot(i) - 2.0_kreal * omg  * self%w(i) / sphere%radius * & 
			( self%xStart(i) * self%v(i) - self%yStart(i) * self%u(i) ) - sphere%g * self%lapSurf(i) )
		self%hStage1(i) = dt * ( - self%hStart(i) * self%divStart(i) )
		self%areaStage1(i) = dt * ( self%divStart(i) * self%areaStart(i) )
	enddo
	
	!
	!	RK Stage 2
	!
	do i = 1, nP
		self%xIn(i) = self%xStart(i) + 0.5_kreal * self%xStage1(i)
		self%yIn(i) = self%yStart(i) + 0.5_kreal * self%yStage1(i)
		self%zIn(i) = self%zStart(i) + 0.5_kreal * self%zStage1(i)
		self%relVortIn(i) = self%relVortStart(i) + 0.5_kreal * self%relVortStage1(i)
		self%divIn(i) = self%divStart(i) + 0.5_kreal * self%divStage1(i)
		self%hIn(i) = self%hStart(i) + 0.5_kreal * self%hStage1(i)
		self%areaIn(i) = self%areaStart(i) + 0.5_kreal * self%areaStage1(i)
	enddo
	
	call SWESphereRHSIntegrals( self%u, self%v, self%w, self%doubleDot, self%lapSurf, self%xIn, self%yIn, self%zIn, &
		self%relVortIn, self%divIn, self%hIn, self%areaIn, topoFn, self%mask, sphere%radius, sphere%pseEps, &
		sphere%mpiParticles)
		
	
	! DEBUG
	if ( associated(sphere%tracers) .AND. sphere%tracers(2)%nDim == 1 ) then
		sphere%tracers(2)%N = nP
		do i = 1, nP
			sphere%tracers(2)%scalar(i) = self%lapSurf(i)
		enddo
	endif
	
	do i = 1, nP
		self%xStage2(i) = dt * self%u(i)
		self%yStage2(i) = dt * self%v(i)
		self%zStage2(i) = dt * self%w(i)
		self%relVortStage2(i) = dt * ( - ( self%relVortIn(i) + 2.0_kreal * omg * self%zIn(i) / sphere%radius ) * &
			self%divIn(i) - 2.0_kreal * omg * self%w(i) / sphere%radius )
		self%divStage2(i) = dt * ( - self%doubleDot(i) - 2.0_kreal * omg * self%w(i) / sphere%radius * &
			( self%xIn(i) * self%v(i) - self%yIn(i) * self%u(i) ) - sphere%g * self%lapSurf(i) )
		self%hStage2(i) = dt * ( - self%hIn(i) * self%divIn(i) )
		self%areaStage2(i) = dt * ( self%areaIn(i) * self%divIn(i) )
	enddo
	
	!
	!	RK Stage 3
	!
	do i = 1, nP
		self%xIn(i) = self%xStart(i) + 0.5_kreal * self%xStage2(i)
		self%yIn(i) = self%yStart(i) + 0.5_kreal * self%yStage2(i)
		self%zIn(i) = self%zStart(i) + 0.5_kreal * self%zStage2(i)
		self%relVortIn(i) = self%relVortStart(i) + 0.5_kreal * self%relVortStage2(i)
		self%divIn(i) = self%divStart(i) + 0.5_kreal * self%divStage2(i)
		self%hIn(i) = self%hStart(i) + 0.5_kreal * self%hStage2(i)
		self%areaIn(i) = self%areaStart(i) + 0.5_kreal * self%areaStage2(i)
	enddo
	
	call SWESphereRHSIntegrals( self%u, self%v, self%w, self%doubleDot, self%lapSurf, self%xIn, self%yIn, self%zIn, &
		self%relVortIn, self%divIn, self%hIn, self%areaIn, topoFn, self%mask, sphere%radius, sphere%pseEps, &
		sphere%mpiParticles)
	
	! DEBUG
	if ( associated(sphere%tracers) .AND. sphere%tracers(3)%nDim == 1 ) then
		sphere%tracers(3)%N = nP
		do i = 1, nP
			sphere%tracers(3)%scalar(i) = self%doubleDot(i)
		enddo
	endif
	
	! DEBUG
	if ( associated(sphere%tracers) .AND. sphere%tracers(4)%nDim == 1 ) then
		sphere%tracers(4)%N = nP
		do i = 1, nP
			sphere%tracers(4)%scalar(i) = self%lapSurf(i)
		enddo
	endif
	
	do i = 1, nP
		self%xStage3(i) = dt * self%u(i)
		self%yStage3(i) = dt * self%v(i)
		self%zStage3(i) = dt * self%w(i)
		self%relVortStage3(i) = dt * ( - ( self%relVortIn(i) + 2.0_kreal * omg * self%zIn(i) / sphere%radius ) * &
			self%divIn(i) - 2.0_kreal * omg * self%w(i) / sphere%radius )
		self%divStage3(i) = dt * ( - self%doubleDot(i) - 2.0_kreal * omg * self%w(i) / sphere%radius * &
			( self%xIn(i) * self%v(i) - self%yIn(i) * self%u(i) ) - sphere%g * self%lapSurf(i) )
		self%hStage3(i) = dt * ( - self%hIn(i) * self%divIn(i) )
		self%areaStage3(i) = dt * ( self%areaIn(i) * self%divIn(i) )
	enddo
	
	!
	!	RK Stage 4
	!
	do i = 1, nP
		self%xIn(i) = self%xStart(i) + self%xStage3(i)
		self%yIn(i) = self%yStart(i) + self%yStage3(i)
		self%zIn(i) = self%zStart(i) + self%zStage3(i)
		self%relVortIn(i) = self%relVortStart(i) + self%relVortStage3(i)
		self%divIn(i) = self%divStart(i) + self%divStage3(i)
		self%hIn(i) = self%hStart(i) + self%hStage3(i)
		self%areaIn(i) = self%areaStart(i) + self%areaStage3(i)
	enddo
	
	call SWESphereRHSIntegrals( self%u, self%v, self%w, self%doubleDot, self%lapSurf, self%xIn, self%yIn, self%zIn, &
		self%relVortIn, self%divIn, self%hIn, self%areaIn, topoFn, self%mask, sphere%radius, sphere%pseEps, &
		sphere%mpiParticles)
		
	
	
	do i = 1, nP
		self%xStage4(i) = dt * self%u(i)
		self%yStage4(i) = dt * self%v(i)
		self%zStage4(i) = dt * self%w(i)
		self%relVortStage4(i) = dt * ( - ( self%relVortIn(i) + 2.0_kreal * omg * self%zIn(i) / sphere%radius ) * &
			self%divIn(i) - 2.0_kreal * omg * self%w(i) / sphere%radius )
		self%divStage4(i) = dt * ( - self%doubleDot(i) - 2.0_kreal * 0.0_kreal * self%w(i) / sphere%radius * &
			( self%xIn(i) * self%v(i) - self%yIn(i) * self%u(i) ) - sphere%g * self%lapSurf(i) )
		self%hStage4(i) = dt * ( - self%hIn(i) * self%divIn(i) )
		self%areaStage4(i) = dt * ( self%areaIn(i) * self%divIn(i) )
	enddo
	
	!
	!	RK update
	!
	do i = 1, nP
		self%xStart(i) = self%xStart(i) + self%xStage1(i) / 6.0_kreal + self%xStage2(i) / 3.0_kreal + &
				self%xStage3(i) / 3.0_kreal + self%xStage4(i) / 6.0_kreal
		self%yStart(i) = self%yStart(i) + self%yStage1(i) / 6.0_kreal + self%yStage2(i) / 3.0_kreal + &
				self%yStage3(i) / 3.0_kreal + self%yStage4(i) / 6.0_kreal		
		self%zStart(i) = self%zStart(i) + self%zStage1(i) / 6.0_kreal + self%zStage2(i) / 3.0_kreal + &
				self%zStage3(i) / 3.0_kreal + self%zStage4(i) / 6.0_kreal
		self%relVortStart(i) = self%relVortStart(i) + self%relVortStage1(i) / 6.0_kreal + self%relVortStage2(i) / 3.0_kreal + &
				self%relVortStage3(i) / 3.0_kreal + self%relVortStage4(i) / 6.0_kreal			
		self%divStart(i) = self%divStart(i) + self%divStage1(i) / 6.0_kreal + self%divStage2(i) / 3.0_kreal + &
				self%divStage3(i) / 3.0_kreal + self%divStage4(i) / 6.0_kreal	
		self%hStart(i) = self%hStart(i) + self%hStage1(i) / 6.0_kreal + self%hStage2(i) / 3.0_kreal + &
				self%hStage3(i) / 3.0_kreal + self%hStage4(i) / 6.0_kreal
		self%areaStart(i) = self%areaStart(i) + self%areaStage1(i) / 6.0_kreal + self%areaStage2(i) / 3.0_kreal + &
				self%areaStage3(i) / 3.0_kreal + self%areaStage4(i) / 6.0_kreal		
	enddo
	
	call SWESphereRHSIntegrals( self%u, self%v, self%w, self%doubleDot, self%lapSurf, self%xStart, self%yStart, &
		self%zStart, self%relVortStart, self%divStart, self%hStart, self%areaStart, topoFn, self%mask, &
		sphere%radius, sphere%pseEps, sphere%mpiParticles)
	
	do i = 1, nP
		sphere%mesh%particles%x(i) = self%xStart(i)
		sphere%mesh%particles%y(i) = self%yStart(i)
		sphere%mesh%particles%z(i) = self%zStart(i)
		sphere%mesh%particles%area(i) = self%areaStart(i)
		sphere%relVort%scalar(i) = self%relVortStart(i)
		sphere%divergence%scalar(i) = self%divStart(i)
		sphere%h%scalar(i) = self%hStart(i)
		sphere%velocity%xComp(i) = self%u(i)
		sphere%velocity%yComp(i) = self%v(i)
		sphere%velocity%zComp(i) = self%w(i)
	enddo
	call SetBottomHeightOnMesh( sphere, topoFn )
end subroutine

!
!----------------
! private methods
!----------------
!
subroutine SWESphereRHSIntegrals( u, v, w, doubleDot, lapSurf, x, y, z, relVort, div, h, area, topoFn, mask, &
	radius, pseEps, mpiParticles)
	real(kreal), dimension(:), intent(inout) :: u
	real(kreal), dimension(:), intent(inout) :: v
	real(kreal), dimension(:), intent(inout) :: w
	real(kreal), dimension(:), intent(inout) :: doubleDot
	real(kreal), dimension(:), intent(inout) :: lapSurf
	real(kreal), dimension(:), intent(in) :: x
	real(kreal), dimension(:), intent(in) :: y
	real(kreal), dimension(:), intent(in) :: z
	real(kreal), dimension(:), intent(in) :: relVort
	real(kreal), dimension(:), intent(in) :: div
	real(kreal), dimension(:), intent(in) :: h
	real(kreal), dimension(:), intent(in) :: area
	procedure(scalarFnOf3dSpace) :: topoFn
	logical(klog), dimension(:), intent(in) :: mask
	real(kreal), intent(in) :: radius
	real(kreal), intent(in) :: pseEps
	type(MPISetup), intent(in) :: mpiParticles
	!
	integer(kint) :: i, j, mpiErrCode
	real(kreal) :: rotStrength, potStrength, denom, denom2, pseKin, lapKernel
	real(kreal) :: surfHeightI, surfHeightJ
	real(kreal) :: fourPiRsq, dotProd, u0
	real(kreal), dimension(3) :: px, py, pz, gradU, gradV, gradW
		
	fourPiRsq = 4.0_kreal * PI * radius * radius
	u0 = 2.0_kreal * PI / 12.0_kreal
	
	do i = mpiParticles%indexStart(procRank), mpiParticles%indexEnd(procRank)
		u(i) = 0.0_kreal
		v(i) = 0.0_kreal
		w(i) = 0.0_kreal
		lapSurf(i) = 0.0_kreal
		surfHeightI = h(i) + topoFn( x(i), y(i), z(i) )
		gradU = 0.0_kreal
		gradV = 0.0_kreal 
		gradW = 0.0_kreal
		px(1) = 1.0_kreal - x(i) * x(i)
		px(2) = - x(i) * y(i)
		px(3) = - x(i) * z(i)
		py(1) = - x(i) * y(i)
		py(2) = 1.0_kreal - y(i) * y(i)
		py(3) = - y(i) * z(i)
		pz(1) = - x(i) * z(i)
		pz(2) = - y(i) * z(i)
		pz(3) = 1.0_kreal - z(i) * z(i)
		do j = 1, i - 1
			if ( mask(j) ) then
				surfHeightJ = h(j) + topoFn(x(j), y(j), z(j))
				denom = radius * radius - (x(i) * x(j) + y(i) * y(j) + z(i) * z(j))
				rotStrength =  relVort(j) * area(j) / denom / fourPiRsq
				potStrength = div(j) * area(j) / denom / fourPiRsq
				
				u(i) = u(i) - ( y(i) * z(j) - z(i) * y(j) ) * rotStrength - radius * x(j) * potStrength
				v(i) = v(i) - ( z(i) * x(j) - x(i) * z(j) ) * rotStrength - radius * y(j) * potStrength
				w(i) = w(i) - ( x(i) * y(j) - y(i) * x(j) ) * rotStrength - radius * z(j) * potStrength
				
				gradU(1) = gradU(1) + x(j) * ( -( y(i) * z(j) - z(i) * y(j) ) * rotStrength - radius * x(j) * potStrength ) / denom
				gradU(2) = gradU(2) - z(j) * rotStrength + y(j) * ( -( y(i) * z(j) - z(i) * y(j)) * rotStrength - &
					radius * x(j) * potStrength) / denom
				gradU(3) = gradU(3) + y(j) * rotStrength + z(j) * ( -( y(i) * z(j) - z(i) * y(j)) * rotStrength - &
					radius * x(j) * potStrength) / denom
			
				gradV(1) = gradV(1) + z(j) * rotStrength + x(j) * ( -( z(i) * x(j) - x(i) * z(j) ) * rotStrength - &
					radius * y(j) * potStrength ) / denom 
			    gradV(2) = gradV(2) + y(j) * ( -(z(i) * x(j) - x(i) * z(j)) * rotStrength - radius * y(j) * potStrength ) / denom
				gradV(3) = gradV(3) - x(j) * rotStrength + z(j) * ( -( z(i) * x(j) - x(i) * z(j) ) * rotStrength - &
					radius * y(j) * potStrength ) /denom
			    
			    gradW(1) = gradW(1) - y(j) * rotStrength + x(j) * ( -(x(i) * y(j) - y(i) * x(j)) * rotStrength - & 
			    	radius * z(j) * potStrength ) / denom
			    gradW(2) = gradW(2) + x(j) * rotStrength + y(j) * ( -(x(i) * y(j) - y(i) * x(j)) * rotStrength - &
			    	radius * z(j) * potStrength ) / denom
			    gradW(3) = gradW(3) + z(j) * ( -(x(i) * y(j) - y(i) * x(j) ) * rotStrength - radius * z(j) * potStrength) / denom
			    		  
			    pseKin = SphereDistance( [x(i), y(i), z(i)], [x(j), y(j), z(j)] ) / pseEps
			    lapKernel = bivariateLaplacianKernel8( pseKin ) / pseEps**2 
				lapSurf(i) = lapSurf(i) + lapKernel * (surfHeightJ - surfHeightI) * area(j)
			endif
		enddo
		do j = i + 1, size(x)
			if ( mask(j) ) then
				surfHeightJ = h(j) + topoFn(x(j), y(j), z(j))
				denom = radius * radius - (x(i) * x(j) + y(i) * y(j) + z(i) * z(j))
				rotStrength =  relVort(j) * area(j) / denom / fourPiRsq
				potStrength = div(j) * area(j) / denom / fourPiRsq
				
				u(i) = u(i) - ( y(i) * z(j) - z(i) * y(j) ) * rotStrength - radius * x(j) * potStrength
				v(i) = v(i) - ( z(i) * x(j) - x(i) * z(j) ) * rotStrength - radius * y(j) * potStrength
				w(i) = w(i) - ( x(i) * y(j) - y(i) * x(j) ) * rotStrength - radius * z(j) * potStrength
				
				gradU(1) = gradU(1) + x(j) * ( -( y(i) * z(j) - z(i) * y(j) ) * rotStrength - radius * x(j) * potStrength ) / denom
				gradU(2) = gradU(2) - z(j) * rotStrength + y(j) * ( -( y(i) * z(j) - z(i) * y(j)) * rotStrength - &
					radius * x(j) * potStrength) / denom
				gradU(3) = gradU(3) + y(j) * rotStrength + z(j) * ( -( y(i) * z(j) - z(i) * y(j)) * rotStrength - &
					radius * x(j) * potStrength) / denom
			
				gradV(1) = gradV(1) + z(j) * rotStrength + x(j) * ( -( z(i) * x(j) - x(i) * z(j) ) * rotStrength - &
					radius * y(j) * potStrength ) / denom 
			    gradV(2) = gradV(2) + y(j) * ( -(z(i) * x(j) - x(i) * z(j)) * rotStrength - radius * y(j) * potStrength ) / denom
				gradV(3) = gradV(3) - x(j) * rotStrength + z(j) * ( -( z(i) * x(j) - x(i) * z(j) ) * rotStrength - &
					radius * y(j) * potStrength ) /denom
			    
			    gradW(1) = gradW(1) - y(j) * rotStrength + x(j) * ( -(x(i) * y(j) - y(i) * x(j)) * rotStrength - & 
			    	radius * z(j) * potStrength ) / denom
			    gradW(2) = gradW(2) + x(j) * rotStrength + y(j) * ( -(x(i) * y(j) - y(i) * x(j)) * rotStrength - &
			    	radius * z(j) * potStrength ) / denom
			    gradW(3) = gradW(3) + z(j) * ( -(x(i) * y(j) - y(i) * x(j) ) * rotStrength - radius * z(j) * potStrength) / denom
			    		  
			    pseKin = SphereDistance( [x(i), y(i), z(i)], [x(j), y(j), z(j)] ) / pseEps
			    lapKernel = bivariateLaplacianKernel8( pseKin ) / pseEps**2 
				lapSurf(i) = lapSurf(i) + lapKernel * (surfHeightJ - surfHeightI) * area(j)
			endif
		enddo

		u(i) = - u0 * y(i)
		v(i) = u0 * x(i)
		w(i) = 0.0_kreal

		lapSurf(i) = lapSurf(i) / pseEps**2
			
!		doubleDot(i) = sum( px * gradU )**2 + sum( py * gradV )**2 + sum( pz * gradW )**2 + 2.0_kreal * ( &
!			sum( py * gradU ) * sum( px * gradV ) + sum( pz * gradU ) * sum( px * gradW ) + &
!			sum( py * gradW ) * sum( pz * gradV ) )

		doubleDot(i) = 0.0_kreal !2.0_kreal * u0 * u0 * ( - z(i) * z(i) )

		if ( procRank == 0 ) then
!			print *, "i = ", i, " p_x = ", px
!			print *, "i = ", i, " p_y = ", py
!			print *, "i = ", i, " p_z = ", pz
!			print *, "i = ", i, " gradU = ", gradU
!			print *, "i = ", i, " gradV = ", gradV
!			print *, "i = ", i, " gradW = ", gradW
!			print *, "i = ", i, " dDot(i) = ", doubleDot(i)
		endif
	enddo
	
	do i = 0, numProcs - 1
		call MPI_BCAST(u(mpiParticles%indexStart(i):mpiParticles%indexEnd(i)), mpiParticles%messageLength(i), &
				MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)
		call MPI_BCAST(v(mpiParticles%indexStart(i):mpiParticles%indexEnd(i)), mpiParticles%messageLength(i), &
				MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)
		call MPI_BCAST(w(mpiParticles%indexStart(i):mpiParticles%indexEnd(i)), mpiParticles%messageLength(i), &
				MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)
		call MPI_BCAST(doubleDot(mpiParticles%indexStart(i):mpiParticles%indexEnd(i)), mpiParticles%messageLength(i), &
				MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)				
		call MPI_BCAST(lapSurf(mpiParticles%indexStart(i):mpiParticles%indexEnd(i)), mpiParticles%messageLength(i), &
				MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)
	enddo
	
	call LogMessage(log, DEBUG_LOGGING_LEVEL, trim(logKey)//" max u = ", maxval(u) )
	call LogMessage(log, DEBUG_LOGGING_LEVEL, trim(logKey)//" min u = ", minval(u) )
	call LogMessage(log, DEBUG_LOGGING_LEVEL, trim(logKey)//" max v = ", maxval(v) )
	call LogMessage(log, DEBUG_LOGGING_LEVEL, trim(logKey)//" min v = ", minval(v) )
	call LogMessage(log, DEBUG_LOGGING_LEVEL, trim(logKey)//" max w = ", maxval(w) )
	call LogMessage(log, DEBUG_LOGGING_LEVEL, trim(logKey)//" min w = ", minval(w) )
	call LogMessage(log, DEBUG_LOGGING_LEVEL, trim(logKey)//" max(ddot) = ", maxval(doubleDot) )
	call LogMessage(log, DEBUG_LOGGING_LEVEL, trim(logKey)//" min(ddot) = ", minval(doubleDot) )
	call LogMessage(log, DEBUG_LOGGING_LEVEL, trim(logKey)//" max(lapS) = ", maxval(lapSurf) )
	call LogMessage(log, DEBUG_LOGGING_LEVEL, trim(logKey)//" min(lapS) = ", minval(lapSurf) )
end subroutine

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
