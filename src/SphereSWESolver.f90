module SphereSWESolverModule
!> @file SphereSWESolver.f90
!> Data structure for solving the Shallow Water Equations (SWE) on a rotating sphere.
!> @author Peter Bosler, Sandia National Laboratories Center for Computing Research
!> 
!>
!> @defgroup SphereSWESolver SphereSWESolver
!> Data structure for solving the Shallow Water Equations (SWE) on a rotating sphere.
!> Used with the @ref SphereWE module.
!>
!> The shallow water equations are presented in the detailed description of the @ref SphereSWE module.  @n
!> 
!> Integrals are used to compute the velocity @f$\vec{u}@f$, the double dot product @f$ \nabla\vec{u}:\nabla\vec{u} @f$, and a PSE integral (see @ref PSEDirectSum) 
!> approximates the Laplacian of the fluid surface.@n
!> As with the other integrals in LPM, these integrals are computed as a parallel direct summation across all MPI ranks.@n
!> 
!> All other terms are computed as ODEs along each particle trajectory.
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
	real(kreal), allocatable, dimension(:) :: xStart !< starting x-coordinate of each particle
	real(kreal), allocatable, dimension(:) :: yStart !< starting y-coordinate of each particle
	real(kreal), allocatable, dimension(:) :: zStart !< starting z-coordinate of each particle
	real(kreal), allocatable, dimension(:) :: relVortStart !< starting relative vorticity of each particle
	real(kreal), allocatable, dimension(:) :: divStart !< starting divergence of each particle
	real(kreal), allocatable, dimension(:) :: hStart !< starting fluid depth at each particle
	real(kreal), allocatable, dimension(:) :: areaStart !< starting area of each particle
	real(kreal), allocatable, dimension(:) :: u !< x-component of fluid velocity
	real(kreal), allocatable, dimension(:) :: v !< y-component of fluid velocity
	real(kreal), allocatable, dimension(:) :: w !< z-component of fluid velocity
	real(kreal), allocatable, dimension(:) :: potVort !< potential vorticity carried by each particle
	real(kreal), allocatable, dimension(:) :: doubleDot !< @f$ \nabla\vec{u}:\nabla\vec{u} @f$ at each particle
	real(kreal), allocatable, dimension(:) :: lapSurf !< Laplacian of the fluid surface at each particle
	logical(klog), allocatable, dimension(:) :: mask !< mask(i) = true if particle i is active
	
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
	real(kreal), allocatable :: zIn(:) !< z-coordinate input to RK4
	real(kreal), allocatable :: zStage1(:) !< z-coordinates of each particle at RK4 stage 1
	real(kreal), allocatable :: zStage2(:) !< z-coordinates of each particle at RK4 stage 2
	real(kreal), allocatable :: zStage3(:) !< z-coordinates of each particle at RK4 stage 3
	real(kreal), allocatable :: zStage4(:) !< z-coordinates of each particle at RK4 stage 4
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
	real(kreal), allocatable :: areaIn(:) !< area input to RK4
	real(kreal), allocatable :: areaStage1(:) !< area of each particle at RK4 stage 1
	real(kreal), allocatable :: areaStage2(:) !< area of each particle at RK4 stage 2
	real(kreal), allocatable :: areaStage3(:) !< area of each particle at RK4 stage 3
	real(kreal), allocatable :: areaStage4(:) !< area of each particle at RK4 stage 4
	
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
	real(kreal) :: fourPiRsq, dotProd
	real(kreal) :: ux, uy, uz, vx, vy, vz, wx, wy, wz

	
	fourPiRsq = 4.0_kreal * PI * radius * radius
		
	do i = mpiParticles%indexStart(procRank), mpiParticles%indexEnd(procRank)
		u(i) = 0.0_kreal
		v(i) = 0.0_kreal
		w(i) = 0.0_kreal
		doubleDot(i) = 0.0_kreal
		lapSurf(i) = 0.0_kreal
		surfHeightI = h(i) + topoFn( x(i), y(i), z(i) )
		ux = 0.0_kreal
		uy = 0.0_kreal
		uz = 0.0_kreal
		vx = 0.0_kreal
		vy = 0.0_kreal
		vz = 0.0_kreal
		wx = 0.0_kreal
		wy = 0.0_kreal
		wz = 0.0_kreal
		do j = 1, size(x)
			if ( mask(j) ) then
				!
				! regularized kernels
				!
				surfHeightJ = h(j) + topoFn( x(j), y(j), z(j) )
				pseKin = SphereDistance( [x(i), y(i), z(i)], [x(j), y(j), z(j)] ) / pseEps
				lapKernel = bivariateLaplacianKernel8( pseKin ) / pseEps**2
				lapSurf(i) = lapSurf(i) + lapKernel * (surfHeightJ - surfHeightI) * area(j)
				
				if ( i == j ) cycle
				
				dotProd = x(i) * x(j) + y(i) * y(j) + z(i) * z(j)
				denom = fourPiRsq * ( radius * radius - dotProd ) 
				denom2 = fourPiRsq * ( radius * radius - dotProd )**2
				
				rotStrength = relVort(j) * area(j) / denom
				potStrength = radius * div(j) * area(j) / denom
				
				u(i) = u(i) - ( y(i) * z(j) - z(i) * y(j) ) * rotStrength - x(j) * potStrength
				v(i) = v(i) - ( z(i) * x(j) - x(i) * z(j) ) * rotStrength - y(j) * potStrength
				w(i) = w(i) - ( x(i) * y(j) - y(i) * x(j) ) * rotStrength - z(j) * potStrength
				
				ux = ux - x(j) * ( ( y(i) * z(j) - z(i) * y(j)  ) * rotStrength - x(j) * potStrength)
			endif
		enddo
	enddo
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

!> @}
end module
