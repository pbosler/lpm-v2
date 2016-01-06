module BIVARInterfaceModule
!> @file BIVARParticlesInterface.f90
!> Interface and workspace for interpolation of LPM data by the BIVAR package.
!> @author Peter Bosler, Sandia National Laboratories Center for Computing Research
!> 
!>
!> @defgroup BIVARInterface BIVARInterface
!> Interface and workspace for interpolation of LPM data by the BIVAR package.
!> 
!> Bivar performs its interpolation using quintic Hermite polynomials on triangles.  
!> Derivatives are estimated as described in reference 2.  
!> It must either build its own Delaunay triangulation of all LPM particles, or, if used with a triangular
!> @ref PolyMesh2d object, the native LPM triangulation may be used.
!>
!> The BIVAR package is due to H. Akima and can be downloaded from the ACM TOMS web page.  
!> LPM uses the updated Fortran 90 version provided by [Jeff Burkhardt](http://people.sc.fsu.edu/~jburkhardt). @n
!> @n
!> The complete description of the BIVAR package may be found in 
!> 1. H. Akima, _A method of bivariate interpolation and smooth surface fitting
!> for values given at irregularly distributed data points_, Office of Telecommunications Report OT 75:70, U. S. Dept. of Commerce, 1975.
!> 2. H. Akima, A method of bivariate interpolation and smooth surface fitting
!> for values given at irregularly distributed data points, _ACM TOMS_, 4:2, 1978.
!> 3. H. Akima, On estimating partial derivatives for bivariate interpolation of scattered data, _Rocky Mtn. J. Math._, 14, 1984. @n
!>
!> @{
use NumberKindsModule
use LoggerModule
use ParticlesModule
use FieldModule
use BIVARModule

implicit none

private
public BIVARInterface, New, Delete, SetBIVARMD
public InterpolateScalar, InterpolateVector, InterpolateLagParam

!
!----------------
! Types and module constants
!----------------
!
type BIVARInterface
	integer(kint), dimension(:), allocatable :: intWork !< Integer variable workspace, see bivar.f90 for requirements
	real(kreal), dimension(:), allocatable :: realWork !< Real variable workspace, see bivar.f90 for requirements
	integer(kint) :: md = 1 !< Flag to determine state of Delaunay triangulation, see bivar.f90 for definition
	
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

!
!----------------
! Logging
!----------------
!
logical(klog), save :: logInit = .FALSE.
type(Logger) :: log
character(len=28), save :: logKey = 'Bivar'
integer(kint), parameter :: logLevel = DEBUG_LOGGING_LEVEL
character(len=MAX_STRING_LENGTH) :: logString

contains

!> @brief Allocates memory for interpolation of LPM data using the Bivar package.
!> 
!> @param[out] self Target interpolation interface
!> @param[in] sourceParticles @ref Particles object
subroutine newPrivate( self, sourceParticles )
	type(BIVARInterface), intent(out) :: self
	type(Particles), intent(in) :: sourceParticles
	
	if ( .NOT. logInit ) call InitLogger(log, procRank)
	
	allocate(self%intWork( 31 * sourceParticles%N + sourceParticles%N_Max))
	allocate(self%realWork( 8 * sourceParticles%N ))
end subroutine

!> @brief Deletes and frees memory associated with a BIVARInterface object.
!> @param[inout] self Target interpolation interface
subroutine deletePrivate( self )
	type(BIVARInterface), intent(inout)  :: self
	
	if ( allocated( self%intWork )) then
		deallocate(self%intWork)
		deallocate(self%realWork)
	endif
end subroutine

!> @brief Resets the value of the Bivar "md" parameter.  
!> See bivar.f90 for its definition.
!> 
!> @param[inout] self Target interpolation interface
!> @param[in] md New parameter value
subroutine SetBIVARMD(self, md)
	type(BIVARInterface), intent(inout) :: self
	integer(kint), intent(in) :: md
	
	if ( md > 3 .OR. md < 1 ) then
		call LogMessage(log, WARNING_LOGGING_LEVEL,trim(logKey)//" SetBIVARMD WARNING : ", "invalid md value")
	else
		!self%md = md
	endif
end subroutine

!> @brief Interpolates a scalar @ref Field from a set of @ref Particles to a set of destination points.
!> 
!> @param[inout] scalarOut Interpolated scalar values, array must have been allocated prior to calling this routine.
!> @param[in] xOut x-coordinates of destination points
!> @param[in] yOut y-coordinates of destination points
!> @param[inout] self Target interpolation interface
!> @param[in] sourceParticles Existing @ref Particles
!> @param[in] sourceField Existing scalar @ref Field
subroutine InterpolateScalar( scalarOut, xOut, yOut, self, sourceParticles, sourceField )
	real(kreal), dimension(:), intent(inout) :: scalarOut
	real(kreal), dimension(:), intent(in) :: xOut
	real(kreal), dimension(:), intent(in) :: yOut
	type(BIVARInterface), intent(inout) :: self
	type(Particles), intent(in) :: sourceParticles
	type(Field), intent(in) :: sourceField
	!
	integer(kint) :: nParticles
	integer(kint) :: nOut
	
	nParticles = sourceParticles%N
	nOut = size(xOut)
	call IDBVIP( self%md, sourceParticles%N, sourceParticles%x(1:nParticles), sourceParticles%y(1:nParticles), &
				 sourceField%scalar(1:nParticles), nOut, xOut, yOut, scalarOut, self%intWork, self%realWork)
end subroutine	

!> @brief Interpolates the Lagrangian coordinates from a set of @ref Particles to a set of destination points.
!> 
!> @param[inout] lagXOut Interpolated Lagrangian x-coordinate values, array must have been allocated prior to calling this routine.
!> @param[inout] lagYOut Interpolated Lagrangian y-coordinate values, array must have been allocated prior to calling this routine.
!> @param[in] xOut x-coordinates of destination points
!> @param[in] yOut y-coordinates of destination points
!> @param[inout] self Target interpolation interface
!> @param[in] sourceParticles Existing @ref Particles
subroutine InterpolateLagParam( lagXOut, lagYOut, xOut, yOut, self, sourceParticles )
	real(kreal), dimension(:), intent(inout) :: lagXOut
	real(kreal), dimension(:), intent(inout) :: lagYOut
	real(kreal), dimension(:), intent(in) :: xOut
	real(kreal), dimension(:), intent(in) :: yOut
	type(BIVARInterface), intent(inout) :: self
	type(Particles), intent(in) :: sourceParticles
	!
	integer(kint) :: nParticles
	integer(kint) :: nOut
	integer(kint) :: tempMD
	
	nParticles = sourceParticles%N
	nOut = size(xOut)
	!call LogMessage(log, DEBUG_LOGGING_LEVEL, trim(logKey)//" nOut = ", nOut)
	!print *, "nOut = ", nOut
	
	tempMD = self%md
	call IDBVIP( self%md, nParticles, sourceParticles%x(1:nParticles), sourceParticles%y(1:nParticles), &
				 sourceParticles%x0(1:nParticles), nOut, xOut, yOut, lagXOut, self%intWork, self%realWork)
	self%md = 3
	call IDBVIP( self%md, nParticles, sourceParticles%x(1:nParticles), sourceParticles%y(1:nParticles), &
				 sourceParticles%y0(1:nParticles), nOut, xOut, yOut, lagYOut, self%intWork, self%realWork)	
	self%md = tempMD			 			 
end subroutine

!> @brief Interpolates a vector @ref Field from a set of @ref Particles to a set of destination points.
!> 
!> @param[inout] vecXOut Interpolated values of the vector field's x component, array must have been allocated prior to calling this routine.
!> @param[inout] vecYOut Interpolated values of the vector field's y component, array must have been allocated prior to calling this routine.
!> @param[in] xOut x-coordinates of destination points
!> @param[in] yOut y-coordinates of destination points
!> @param[inout] self Target interpolation interface
!> @param[in] sourceParticles Existing @ref Particles
!> @param[in] sourceField Existing vector @ref Field
subroutine InterpolateVector( vecXOut, vecYOut, xOut, yOUt, self, sourceParticles, sourceField )
	real(kreal), dimension(:), intent(inout) :: vecXOut
	real(kreal), dimension(:), intent(inout) :: vecYOut
	real(kreal), dimension(:), intent(in) :: xOut
	real(kreal), dimension(:), intent(in) :: yOut
	type(BIVARInterface), intent(inout) :: self
	type(Particles), intent(in) :: sourceParticles
	type(Field), intent(in) :: sourceField
	!
	integer(kint) :: nParticles
	integer(kint) :: nOut
	integer(kint) :: tempMD
	
	nParticles = sourceParticles%n
	nOut = size(xOut)
	tempMD = self%md
	
	call IDBVIP( self%md, nParticles, sourceParticles%x(1:nParticles), sourceParticles%y(1:nParticles), &
				 sourceField%xComp(1:nParticles), nout, xOut, yOut, vecXOut, self%intWork, self%realWork)
	self%md = 3
	call IDBVIP( self%md, nParticles, sourceParticles%x(1:nParticles), sourceParticles%y(1:nParticles), &
				 sourceField%yComp(1:nParticles), nOut, xOut, yOut, vecYOut, self%intWork, self%realWork)
	self%md = tempMD
end subroutine	

!
!----------------
! private methods
!----------------
!

!> @brief Initializes a logger for the BIVARInterface module
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
