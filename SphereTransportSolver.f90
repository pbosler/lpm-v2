module SphereTransportSolverModule
!> @file SphereTransportSolver.f90
!> Data structure for representing solutions of the Barotropic Vorticity Equation (BVE) on the surface of a rotating sphere
!> @author Peter Bosler, Sandia National Laboratories Center for Computing Research
!> 
!>
!> @defgroup SphereTransportSolver SphereTransportSolver
!> @brief Data structure for representing solutions of the advection equation on the surface of a rotating sphere.
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
use SphereGeomModule
use SphereTransportModule

implicit none

include 'mpif.h'

private
!
!----------------
! Module types and public declarations
!----------------
!
public TransportSolver, New, Delete
public Timestep
public LauritzenEtalDeformationalVelocity


type TransportSolver
	! Input state : copy from host to device
	real(kreal), allocatable, dimension(:) :: xStart
	real(kreal), allocatable, dimension(:) :: yStart
	real(kreal), allocatable, dimension(:) :: zStart
	real(kreal), allocatable, dimension(:) :: areaStart
	real(kreal), allocatable, dimension(:) :: densityStart
	logical(klog), allocatable, dimension(:) :: mask

	! DEVICE	
	! Knowns
	real(kreal), allocatable, dimension(:) :: u
	real(kreal), allocatable, dimension(:) :: v
	real(kreal), allocatable, dimension(:) :: w
	real(kreal), allocatable, dimension(:) :: div

	! RK4
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
	real(kreal), allocatable, dimension(:) :: areaIn
	real(kreal), allocatable, dimension(:) :: areaStage1
	real(kreal), allocatable, dimension(:) :: areaStage2
	real(kreal), allocatable, dimension(:) :: areaStage3
	real(kreal), allocatable, dimension(:) :: areaStage4
	real(kreal), allocatable, dimension(:) :: densityIn
	real(kreal), allocatable, dimension(:) :: densityStage1
	real(kreal), allocatable, dimension(:) :: densityStage2
	real(kreal), allocatable, dimension(:) :: densityStage3
	real(kreal), allocatable, dimension(:) :: densityStage4
	
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
character(len=28), save :: logKey = 'TransportSolver'
integer(kint), parameter :: logLevel = DEBUG_LOGGING_LEVEL

contains
!
!----------------
! public methods
!----------------
!
function LauritzenEtalDeformationalVelocity( x, y, z, t )
	real(kreal), dimension(3) :: LauritzenEtalDeformationalVelocity
	real(kreal), intent(in) :: x, y, z, t
	real(kreal), parameter :: RR = 1.0_kreal
	real(kreal), parameter :: TT = 5.0_kreal
	real(kreal) :: u, v, lat, lon
	
	lat = Latitude(x, y, z)
	lon = Longitude(x, y, z)
	u = 10.0_kreal * RR / TT * sin(lon - 2.0_kreal * PI * t / TT) * sin(lon - 2.0_kreal * PI * t / TT ) * &
		sin(2.0_kreal * lat ) * cos( PI * t / TT) + 2.0_kreal * PI * RR / TT * cos(lat)
	v = 10.0_kreal * RR / TT * sin(2.0_kreal * (lon - 2.0_kreal * PI * t / TT)) * cos(lat) * cos(PI * t/TT )
	
	LauritzenEtalDeformationalVelocity(1) = -u * sin(lon) - v * sin(lat) * cos(lon)
	LauritzenEtalDeformationalVelocity(2) =  u * cos(lon) - v * sin(lat) * sin(lon)
	LauritzenEtalDeformationalVelocity(3) =  v * cos(lat)
end function

function LauritzenEtalDivergentFlowVelocity( x, y, z, t )
	real(kreal), dimension(3) :: LauritzenEtalDivergentFlowVelocity
	real(kreal), intent(in) :: x, y, z, t
	real(kreal), parameter :: RR = 1.0_kreal
	real(kreal), parameter :: TT = 5.0_kreal
	real(kreal), parameter :: relPeriod = RR / TT
	real(kreal) :: u, v, lat, lon, lamPrime
	
	lat = Latitude(x, y, z)
	lon = Longitude(x, y, z)
	lamPrime = lon - 2.0_kreal * PI * t / TT
	
	u = - 5.0_kreal * relPeriod * sin( 0.5_kreal * lamPrime) * sin( 0.5_kreal * lamPrime ) * sin( 2.0_kreal * lat ) * &
		cos(lat) * cos(lat) * cos( PI * t / TT ) + 2.0_kreal * PI * relPeriod * cos(lat)
	v = 2.5_kreal * relPeriod * sin(lamPrime) * cos(lat) * cos(lat) * cos(lat) * cos( PI * t / TT )
	
	LauritzenEtalDivergentFlowVelocity(1) = -u * sin(lon) - v * sin(lat) * cos(lon)
	LauritzenEtalDivergentFlowVelocity(2) =  u * cos(lon) - v * sin(lat) * sin(lon)
	LauritzenEtalDivergentFlowVelocity(3) =  v * cos(lat)
end function

function LauritzenEtalDivergentFlowDivergence( x, y, z, t )
	real(kreal) :: LauritzenEtalDivergentFlowDivergence
	real(kreal), intent(in) :: x, y, z, t
	real(kreal), parameter :: TT = 5.0_kreal
	real(kreal) :: u, v, lat, lon, lamPrime
	
	lat = Latitude(x, y, z)
	lon = Longitude(x, y, z)
	lamPrime = lon - 2.0_kreal * PI * t / TT
	
	LauritzenEtalDivergentFlowDivergence = - 15.0_kreal * cos(lat) * cos(lat) * cos( PI * t / TT ) * sin(lat) * &
		sin(lamPrime) / TT
end function 

function RH4Velocity( x, y, z, t )
	real(kreal), dimension(3) :: RH4Velocity
	real(kreal), intent(in) :: x, y, z, t
	real(kreal) :: u, v, lat, lon
	real(kreal), parameter :: u0 = 0.0_kreal
	real(kreal), parameter :: amp = 1.0_kreal
	
	lat = Latitude(x, y, z)
	lon = Longitude(x, y, z)
	
	u = cos(lat) * ( u0 - 0.5_kreal * amp * cos( 4.0_kreal * lon ) * cos(lat) * cos(lat) * &
		(-3.0_kreal + 5.0_kreal * cos(2.0_kreal * lat )) )
	v = -4.0_kreal * amp * cos(lat) * cos(lat) * cos(lat) * sin(4.0_kreal * lon ) * sin(lat)
	
	RH4Velocity(1) = -u * sin(lon) - v * sin(lat) * cos(lon)
	RH4Velocity(2) =  u * cos(lon) - v * sin(lat) * sin(lon)
	RH4Velocity(3) =  v * cos(lat)
end function	

!
!----------------
! private methods
!----------------
!
subroutine newPrivate(self, sphereTrans )
	type(TransportSolver), intent(out) :: self
	type(TransportMesh), intent(in) :: sphereTrans
	integer(kint) :: nP
	
	if ( .NOT. logInit ) call InitLogger(log, procRank)
	
	nP = sphereTrans%mesh%particles%N
	
	allocate(self%xStart(nP))
	allocate(self%yStart(nP))
	allocate(self%zStart(nP))
	allocate(self%mask(nP))
	
	if ( sphereTrans%nonzeroDivergence ) then
		allocate(self%areaStart(nP))
		allocate(self%densityStart(nP))
		allocate(self%densityIn(nP))
		allocate(self%densityStage1(nP))
		allocate(self%densityStage2(nP))
		allocate(self%densityStage3(nP))
		allocate(self%densityStage4(nP))
		allocate(self%areaIn(nP))
		allocate(self%areaStage1(nP))
		allocate(self%areaStage2(nP))
		allocate(self%areaStage3(nP))
		allocate(self%areaStage4(nP))
	endif
	
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
	
	self%mask = sphereTrans%mesh%particles%isActive(1:nP)
end subroutine

subroutine deletePrivate(self)
	type(TransportSolver), intent(inout) :: self
	
	if ( allocated(self%xStart) ) then
		deallocate(self%xStart)
		deallocate(self%yStart)
		deallocate(self%zStart)
		deallocate(self%mask)
	
		if ( allocated(self%areaStart) ) then
			deallocate(self%areaStart)
			deallocate(self%densityStart)
			deallocate(self%densityIn)
			deallocate(self%densityStage1)
			deallocate(self%densityStage2)
			deallocate(self%densityStage3)
			deallocate(self%densityStage4)
			deallocate(self%areaIn)
			deallocate(self%areaStage1)
			deallocate(self%areaStage2)
			deallocate(self%areaStage3)
			deallocate(self%areaStage4)
		endif
	
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
	endif
end subroutine

subroutine timestepPrivate( self, sphere, t, dt, velFn, divFn )
	type(TransportSolver), intent(inout) :: self
	type(TransportMesh), intent(inout) :: sphere
	real(kreal), intent(in) :: t
	real(kreal), intent(in) :: dt
	procedure(vectorFnOf3DSpaceAndTime) :: velFn
	procedure(scalarFnOf3DSpaceAndTime), optional :: divFn
	integer(kint) :: i, nP
	real(kreal), dimension(3) :: vec
	real(kreal) :: div, tIn
	integer(kint) :: mpiErrCode
	
	nP = sphere%mesh%particles%N
	
	!
	!	setup starting point (Later: copy-to-device )
	!
	do i = sphere%mpiParticles%indexStart(procRank), sphere%mpiParticles%indexEnd(procRank)
		self%xStart(i) = sphere%mesh%particles%x(i)
		self%yStart(i) = sphere%mesh%particles%y(i)
		self%zStart(i) = sphere%mesh%particles%z(i)
	enddo
	if ( sphere%nonzeroDivergence ) then
		do i = sphere%mpiParticles%indexStart(procRank), sphere%mpiParticles%indexEnd(procRank)
			self%densityStart(i) = sphere%density%scalar(i)
			if ( self%mask(i) ) then
				 self%areaStart(i) = sphere%mesh%particles%area(i)
			else
				self%areaStart(i) = 0.0_kreal
			endif
		enddo
	endif
!
!	Later : on-device
!	
	!
	!	RK Stage 1
	!
	do i = sphere%mpiParticles%indexStart(procRank), sphere%mpiParticles%indexEnd(procRank)
		vec = velFn( self%xStart(i), self%yStart(i), self%zStart(i), t )
		self%xStage1(i) = dt * vec(1)
		self%yStage1(i) = dt * vec(2)
		self%zStage1(i) = dt * vec(3)
	enddo
	if ( sphere%nonzeroDivergence ) then
		do i = sphere%mpiParticles%indexStart(procRank), sphere%mpiParticles%indexEnd(procRank)
			div = divFn( self%xStart(i), self%yStart(i), self%zStart(i), t )
			self%densityStage1(i) = - dt * div * self%densityStart(i)
			if ( self%mask(i) )  then 
				self%areaStage1(i) = dt * div * self%areaStart(i)			
			else
				self%areaStage1(i) = 0.0_kreal
			endif
		enddo
	endif
	
	!
	!	RK Stage 2
	!
	tIn = t + 0.5_kreal * dt
	do i = sphere%mpiParticles%indexStart(procRank), sphere%mpiParticles%indexEnd(procRank)
		self%xIn(i) = self%xStart(i) + 0.5_kreal * self%xStage1(i)
		self%yIn(i) = self%yStart(i) + 0.5_kreal * self%yStage1(i)
		self%zIn(i) = self%zStart(i) + 0.5_kreal * self%zStage1(i)
	enddo
	if ( sphere%nonzeroDivergence) then
		do i = sphere%mpiParticles%indexStart(procRank), sphere%mpiParticles%indexEnd(procRank)
			self%densityIn(i) = self%densityStart(i) + 0.5_kreal * self%densityStage1(i)
			if ( self%mask(i) ) then
				self%areaIn(i) = self%areaStart(i) + 0.5_kreal * self%areaStage1(i)
			else
				self%areaIn(i) = 0.0_kreal
			endif
		enddo
	endif

	do i = sphere%mpiParticles%indexStart(procRank), sphere%mpiParticles%indexEnd(procRank)
		vec = velFn( self%xIn(i), self%yIn(i), self%zIn(i), tIn )
		self%xStage2(i) = dt * vec(1)
		self%yStage2(i) = dt * vec(2)
		self%zStage2(i) = dt * vec(3)
	enddo	
	if ( sphere%nonzeroDivergence ) then
		do i = sphere%mpiParticles%indexStart(procRank), sphere%mpiParticles%indexEnd(procRank)
			div = divFn( self%xIn(i), self%yIn(i), self%zIn(i), tIn )
			self%densityStage2(i) = - dt * div * self%densityIn(i)
			if ( self%mask(i) ) then 
				self%areaStage2(i) = dt * div * self%areaIn(i)
			else
				self%areaStage2(i) = 0.0_kreal
			endif
		enddo
	endif	
		
	!
	!	RK Stage 3
	!
	tIn = t + 0.5_kreal * dt
	do i = sphere%mpiParticles%indexStart(procRank), sphere%mpiParticles%indexEnd(procRank)
		self%xIn(i) = self%xStart(i) + 0.5_kreal * self%xStage2(i)
		self%yIn(i) = self%yStart(i) + 0.5_kreal * self%yStage2(i)
		self%zIn(i) = self%zStart(i) + 0.5_kreal * self%zStage2(i)
	enddo
	if ( sphere%nonzeroDivergence) then
		do i = sphere%mpiParticles%indexStart(procRank), sphere%mpiParticles%indexEnd(procRank)
			self%densityIn(i) = self%densityStart(i) + 0.5_kreal * self%densityStage2(i)
			if ( self%mask(i) ) then
				 self%areaIn(i) = self%areaStart(i) + 0.5_kreal * self%areaStage2(i)
			else
				self%areaIn(i) = 0.0_kreal
			endif
		enddo
	endif

	do i = sphere%mpiParticles%indexStart(procRank), sphere%mpiParticles%indexEnd(procRank)
		vec = velFn( self%xIn(i), self%yIn(i), self%zIn(i), tIn )
		self%xStage3(i) = dt * vec(1)
		self%yStage3(i) = dt * vec(2)
		self%zStage3(i) = dt * vec(3)
	enddo	
	if ( sphere%nonzeroDivergence ) then
		do i = sphere%mpiParticles%indexStart(procRank), sphere%mpiParticles%indexEnd(procRank)
			div = divFn( self%xIn(i), self%yIn(i), self%zIn(i), tIn )	
			self%densityStage3(i) = - dt * div * self%densityIn(i)
			if ( self%mask(i) )  then
				self%areaStage3(i) = dt * div * self%areaIn(i)
			else
				self%areaStage3(i) = 0.0_kreal
			endif
		enddo
	endif	
	
	!
	!	RK Stage 4
	!
	tIn = t + dt
	do i = sphere%mpiParticles%indexStart(procRank), sphere%mpiParticles%indexEnd(procRank)
		self%xIn(i) = self%xStart(i) + self%xStage3(i)
		self%yIn(i) = self%yStart(i) + self%yStage3(i)
		self%zIn(i) = self%zStart(i) + self%zStage3(i)
	enddo
	if ( sphere%nonzeroDivergence) then
		do i = sphere%mpiParticles%indexStart(procRank), sphere%mpiParticles%indexEnd(procRank)
			self%densityIn(i) = self%densityStart(i) + self%densityStage3(i)
			if ( self%mask(i) ) then
				self%areaIn(i) = self%areaStart(i) + self%areaStage3(i)
			else
				self%areaIn(i) = 0.0_kreal
			endif
		enddo
	endif

	do i = sphere%mpiParticles%indexStart(procRank), sphere%mpiParticles%indexEnd(procRank)
		vec = velFn( self%xIn(i), self%yIn(i), self%zIn(i), tIn )
		self%xStage4(i) = dt * vec(1)
		self%yStage4(i) = dt * vec(2)
		self%zStage4(i) = dt * vec(3)
	enddo	
	if ( sphere%nonzeroDivergence ) then
		do i = sphere%mpiParticles%indexStart(procRank), sphere%mpiParticles%indexEnd(procRank)
			div = divFn( self%xIn(i), self%yIn(i), self%zIn(i), tIn )
			self%densityStage4(i) = - dt * div * self%densityIn(i)
			if ( self%mask(i) ) then
				self%areaStage4(i) = dt * div * self%areaIn(i)
			else
				self%areaStage4(i) = 0.0_kreal
			endif
		enddo
	endif
	
	!
	!	RK update
	!
	do i = sphere%mpiParticles%indexStart(procRank), sphere%mpiParticles%indexEnd(procRank)
		sphere%mesh%particles%x(i) = self%xStart(i) + self%xStage1(i) / 6.0_kreal + self%xStage2(i) / 3.0_kreal + &
				self%xStage3(i) / 3.0_kreal + self%xStage4(i) / 6.0_kreal
		sphere%mesh%particles%y(i) = self%yStart(i) + self%yStage1(i) / 6.0_kreal + self%yStage2(i) / 3.0_kreal + &
				self%yStage3(i) / 3.0_kreal + self%yStage4(i) / 6.0_kreal
		sphere%mesh%particles%z(i) = self%zStart(i) + self%zStage1(i) / 6.0_kreal + self%zStage2(i) / 3.0_kreal + &
				self%zStage3(i) / 3.0_kreal + self%zStage4(i) / 6.0_kreal
	enddo
	if ( sphere%nonzeroDivergence ) then
		do i = sphere%mpiParticles%indexStart(procRank), sphere%mpiParticles%indexEnd(procRank)
			sphere%mesh%particles%area(i) = self%areaStart(i) + self%areaStage1(i) / 6.0_kreal + & 
					self%areaStage2(i) / 3.0_kreal + self%areaStage3(i) / 3.0_kreal + self%areaStage4(i) / 6.0_kreal
			sphere%density%scalar(i) = self%densityStart(i) + self%densityStage1(i) / 6.0_kreal + &
				self%densityStage2(i) / 3.0_kreal + self%densityStage3(i) / 3.0_kreal + self%densityStage4(i)/6.0_kreal
		enddo
	endif
	
	do i = 0, numProcs - 1
		call MPI_BCAST( sphere%mesh%particles%x( sphere%mpiParticles%indexStart(i):sphere%mpiParticles%indexEnd(i)), &
				sphere%mpiParticles%messageLength(i), MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)
		call MPI_BCAST( sphere%mesh%particles%y( sphere%mpiParticles%indexStart(i):sphere%mpiParticles%indexEnd(i)), &
				sphere%mpiParticles%messageLength(i), MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)
		call MPI_BCAST( sphere%mesh%particles%z( sphere%mpiParticles%indexStart(i):sphere%mpiParticles%indexEnd(i)), &
				sphere%mpiParticles%messageLength(i), MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)
		if ( sphere%nonzeroDivergence ) then
			call MPI_BCAST( &
					sphere%mesh%particles%area(sphere%mpiParticles%indexStart(i):sphere%mpiParticles%indexEnd(i)), &
					sphere%mpiParticles%messageLength(i), MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)
			call MPI_BCAST( sphere%density%scalar(sphere%mpiParticles%indexStart(i):sphere%mpiParticles%indexEnd(i)), &
					sphere%mpiParticles%messageLength(i), MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)
		endif
	enddo
	
	do i = 1, nP
		vec = velFn( sphere%mesh%particles%x(i),  sphere%mesh%particles%y(i), sphere%mesh%particles%z(i), t + dt)
		sphere%velocity%xComp(i) = vec(1)
		sphere%velocity%yComp(i) = vec(2)
		sphere%velocity%zComp(i) = vec(3)
	enddo
end subroutine



!> @brief Initializes a logger for the SphereTransport solver module
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


!> }
end module
