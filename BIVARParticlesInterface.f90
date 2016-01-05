module BIVARInterfaceModule

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
	integer(kint), dimension(:), allocatable :: intWork 
	real(kreal), dimension(:), allocatable :: realWork 
	integer(kint) :: md = 1
	
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

subroutine newPrivate( self, sourceParticles )
	type(BIVARInterface), intent(out) :: self
	type(Particles), intent(in) :: sourceParticles
	
	if ( .NOT. logInit ) call InitLogger(log, procRank)
	
	allocate(self%intWork( 31 * sourceParticles%N + sourceParticles%N_Max))
	allocate(self%realWork( 8 * sourceParticles%N ))
end subroutine

subroutine deletePrivate( self )
	type(BIVARInterface), intent(inout)  :: self
	
	if ( allocated( self%intWork )) then
		deallocate(self%intWork)
		deallocate(self%realWork)
	endif
end subroutine

subroutine SetBIVARMD(self, md)
	type(BIVARInterface), intent(inout) :: self
	integer(kint), intent(in) :: md
	
	if ( md > 3 .OR. md < 1 ) then
		call LogMessage(log, WARNING_LOGGING_LEVEL,trim(logKey)//" SetBIVARMD WARNING : ", "invalid md value")
	else
		!self%md = md
	endif
end subroutine

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
