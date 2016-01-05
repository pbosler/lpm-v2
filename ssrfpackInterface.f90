module SSRFPACKInterfaceModule

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
public SSRFPACKInterface, DelaunayTriangulation, New, Delete
public SetScalarSourceData, SetVectorSourceData, SetSourceLagrangianParameter
public InterpolateScalar, InterpolateVector, InterpolateLagParam
public InterpolateScalarToUnifLatLonGrid, InterpolateVectorToUnifLatLonGrid
public SetSigmaFlag

type DelaunayTriangulation
	integer(kint), dimension(:), allocatable :: list 
	integer(kint), dimension(:), allocatable :: lptr 
	integer(kint), dimension(:), allocatable :: lend 
	
	contains
	
		final :: deleteDelTri
end type

type SSRFPACKInterface
	real(kreal), dimension(:,:), allocatable :: grad1 
	real(kreal), dimension(:,:), allocatable :: grad2 
	real(kreal), dimension(:,:), allocatable :: grad3 
	real(kreal), dimension(:), allocatable :: sigma1 
	real(kreal), dimension(:), allocatable :: sigma2 
	real(kreal), dimension(:), allocatable :: sigma3 
	
	contains
	
		final :: deletePrivate
end type

integer(kint), save :: SIGMA_FLAG = 0
integer(kint), parameter :: GRAD_FLAG = 1
integer(kint), save :: startTriangle = 1
real(kreal), save :: sigmaTol = 0.01_kreal

!
!----------------
! Module interfaces
!----------------
!
interface New
	module procedure newPrivate
	module procedure newDelTri
end interface

interface Delete
	module procedure deletePrivate
	module procedure deleteDelTri
end interface


!
!----------------
! Logging
!----------------
!
logical(klog), save :: logInit = .FALSE.
type(Logger) :: log
character(len=28), save :: logKey = 'SSRFPACK'
integer(kint), parameter :: logLevel = DEBUG_LOGGING_LEVEL
character(len=MAX_STRING_LENGTH) :: logString

contains
!
!----------------
! public methods
!----------------
!
subroutine newPrivate(self, aMesh, nDim )
	type(SSRFPACKInterface), intent(out) :: self
	type(PolyMesh2d), intent(inout) :: aMesh
	integer(kint), intent(in) :: nDim 

	if ( .NOT. logInit) call InitLogger(log, procRank)
	
	call LogMessage(log, DEBUG_LOGGING_LEVEL, trim(logKey)//" newSSRFPACKInterface : ", "entering.")
	
	allocate(self%grad1(3, aMesh%particles%N))
	allocate(self%sigma1( 6 * aMesh%particles%N - 12 ))
	self%sigma1 = 0.0_kreal
	
	if ( nDim == 3 ) then
		allocate(self%grad2(3, aMesh%particles%N))
		allocate(self%grad3(3, aMesh%particles%N))
		allocate(self%sigma2( 6 * aMesh%particles%N - 12 ))
		allocate(self%sigma3( 6 * aMesh%particles%N - 12 ))
		self%sigma2 = 0.0_kreal
		self%sigma3 = 0.0_kreal
	endif
	
	call LogMessage(log, DEBUG_LOGGING_LEVEL, trim(logKey)//" newSSRFPACKInterface : ", "returning.")
end subroutine

subroutine newDelTri(self, aMesh)
	type(DelaunayTriangulation), intent(out) :: self
	type(PolyMesh2d), intent(inout) :: aMesh
	
	if ( .NOT. logInit) call InitLogger(log, procRank)
	
	call ProjectParticlesToSphere(aMesh, 1.0_kreal)
	
	allocate(self%list( 6 * aMesh%particles%N - 12 ))
	allocate(self%lptr( 6 * amesh%particles%N - 12 ))
	allocate(self%lend( aMesh%particles%N ))
	
	call BuildDelaunayTriangulation(self, aMesh)
end subroutine

subroutine deleteDelTri(self)
	type(DelaunayTriangulation), intent(inout) :: self
	if ( allocated(self%list)) then
		deallocate(self%list)
		deallocate(self%lptr)
		deallocate(self%lend)
	endif
end subroutine

subroutine deletePrivate( self )
	type(SSRFPACKInterface), intent(inout) :: self
	if ( allocated(self%grad1) ) then
		deallocate(self%grad1)
		deallocate(self%sigma1)
	endif
	if ( allocated(self%grad2) ) then
		deallocate(self%grad2)
		deallocate(self%grad3)
		deallocate(self%sigma2)
		deallocate(self%sigma3)
	endif
end subroutine



subroutine SetVectorSourceData(self, aMesh, delTri, vectorField)
	type(SSRFPACKInterface), intent(inout) :: self
	type(PolyMesh2d), intent(in) :: aMesh
	type(DelaunayTriangulation), intent(in) :: delTri
	type(Field), intent(in) :: vectorField
	integer(kint) :: i, errCode
	real(kreal) :: dSig
	
	do i = 1, aMesh%particles%N
		call GRADL( aMesh%particles%N, i, aMesh%particles%x, aMesh%particles%y, aMesh%particles%z, &
				vectorField%xComp, delTri%list, delTri%lptr, delTri%lend, self%grad1(:,i), errCode)
		call GRADL( aMesh%particles%N, i, aMesh%particles%x, aMesh%particles%y, aMesh%particles%z, &
				vectorField%yComp, delTri%list, delTri%lptr, delTri%lend, self%grad2(:,i), errCode)
		call GRADL( aMesh%particles%N, i, aMesh%particles%x, aMesh%particles%y, aMesh%particles%z, &
				vectorField%zComp, delTri%list, delTri%lptr, delTri%lend, self%grad3(:,i), errCode)
	enddo
	
	if ( SIGMA_FLAG > 0 ) then
	call GETSIG(amesh%particles%N, aMesh%particles%x, aMesh%particles%y, aMesh%particles%z, &
			vectorField%xComp, delTri%list, delTri%lptr, delTri%lend, self%grad1, sigmaTol, self%sigma1, dSig, errCode)
	call GETSIG(amesh%particles%N, aMesh%particles%x, aMesh%particles%y, aMesh%particles%z, &
			vectorField%yComp, delTri%list, delTri%lptr, delTri%lend, self%grad2, sigmaTol, self%sigma2, dSig, errCode)
	call GETSIG(amesh%particles%N, aMesh%particles%x, aMesh%particles%y, aMesh%particles%z, &
			vectorField%zComp, delTri%list, delTri%lptr, delTri%lend, self%grad3, sigmaTol, self%sigma3, dSig, errCode)						
	endif
end subroutine

subroutine SetScalarSourceData(self, aMesh, delTri, scalarField)
	type(SSRFPACKInterface), intent(inout) :: self
	type(PolyMesh2d), intent(in) :: aMesh
	type(DelaunayTriangulation), intent(in) :: delTri
	type(Field), intent(in) :: scalarField
	integer(kint) :: i, errCode
	real(kreal) :: dSig
	
	call LogMessage(log, DEBUG_LOGGING_LEVEL, trim(logKey)//" SetScalarSourceData : ", "entering.")
	
	do i = 1, aMesh%particles%N
		call GRADL( aMesh%particles%n, i, aMesh%particles%x, aMesh%particles%y, aMesh%particles%z, &
					scalarField%scalar, delTri%list, delTri%lptr, delTri%lend, self%grad1(:,i), errCode)
	enddo
	call LogMessage(log, DEBUG_LOGGING_LEVEL, trim(logKey)//" SetScalarSourceData : ", "gradient estimates done.")
	if ( SIGMA_FLAG > 0 ) then
	call GETSIG(aMesh%particles%N, aMesh%particles%x, aMesh%particles%y, aMesh%particles%z, &
			scalarField%scalar, delTri%list, delTri%lptr, delTri%lend, self%grad1, sigmaTol, self%sigma1, dSig, errCode)
	endif
	call LogMessage(log, DEBUG_LOGGING_LEVEL, trim(logKey)//" SetScalarSourceData : ", "returning.")
end subroutine

subroutine SetSourceLagrangianParameter( self, aMesh, delTri )
	type(SSRFPACKInterface), intent(inout) :: self
	type(PolyMesh2d), intent(in) :: aMesh
	type(DelaunayTriangulation), intent(in) :: delTri
	integer(kint) :: i, errCode
	real(kreal) :: dSig
	
	do i = 1, aMesh%particles%N
		call GRADL( aMesh%particles%N, i, aMesh%particles%x, aMesh%particles%y, aMesh%particles%z, &
				aMesh%particles%x0, delTri%list, delTri%lptr, delTri%lend, self%grad1(:,i), errCode)
		call GRADL( aMesh%particles%N, i, aMesh%particles%x, aMesh%particles%y, aMesh%particles%z, &
				aMesh%particles%y0, delTri%list, delTri%lptr, delTri%lend, self%grad2(:,i), errCode)
		call GRADL( aMesh%particles%N, i, aMesh%particles%x, aMesh%particles%y, aMesh%particles%z, &
				aMesh%particles%z0, delTri%list, delTri%lptr, delTri%lend, self%grad3(:,i), errCode)
	enddo
	
	if ( SIGMA_FLAG > 0 ) then
	call GETSIG(amesh%particles%N, aMesh%particles%x, aMesh%particles%y, aMesh%particles%z, &
			aMesh%particles%x0, delTri%list, delTri%lptr, delTri%lend, self%grad1, sigmaTol, self%sigma1, dSig, errCode)
	call GETSIG(amesh%particles%N, aMesh%particles%x, aMesh%particles%y, aMesh%particles%z, &
			aMesh%particles%y0, delTri%list, delTri%lptr, delTri%lend, self%grad2, sigmaTol, self%sigma2, dSig, errCode)
	call GETSIG(amesh%particles%N, aMesh%particles%x, aMesh%particles%y, aMesh%particles%z, &
			aMesh%particles%z0, delTri%list, delTri%lptr, delTri%lend, self%grad3, sigmaTol, self%sigma3, dSig, errCode)
	endif
end subroutine

function InterpolateScalar( lon, lat, self, aMesh, delTri, scalarField)
	real(kreal) :: InterpolateScalar
	type(SSRFPACKInterface), intent(in) :: self
	type(PolyMesh2d), intent(in) :: aMesh
	type(DelaunayTriangulation), intent(in) :: delTri
	type(Field), intent(in) :: scalarField
	real(kreal), intent(in) :: lon
	real(kreal), intent(in) :: lat
	integer(kint) :: errCode
	
	call INTRC1( aMesh%particles%N, lat, lon, amesh%particles%x, amesh%particles%y, amesh%particles%z, &
			 	 scalarField%scalar, delTri%list, delTri%lptr, delTri%lend, SIGMA_FLAG, self%sigma1, GRAD_FLAG, &
			 	 self%grad1, startTriangle, InterpolateScalar, errCode)
end function

subroutine InterpolateScalarToUnifLatLonGrid( self, aMesh, delTri, scalarField, lons, lats, interpOut)
	type(SSRFPACKInterface), intent(inout) :: self
	type(PolyMesh2d), intent(in) :: aMesh
	type(DelaunayTriangulation), intent(in) :: delTri
	type(Field), intent(in) :: scalarField
	real(kreal), dimension(:), intent(in) :: lons
	real(kreal), dimension(:), intent(in) :: lats
	real(kreal), dimension(:,:), intent(out) :: interpOut
	!
	type(MPISetup) :: mpiLons
	integer(kint) :: i, j, nLat, nLon, errCode
	
	nLon = size(lons)
	nLat = size(lats)
	call New(mpiLons, nLon, numProcs)
	
	do j = mpiLons%indexStart(procRank), mpiLons%indexEnd(procrank)
		do i = 1, nLat
			call INTRC1( aMesh%particles%N, lats(i), lons(j), aMesh%particles%x, aMesh%particles%y, aMesh%particles%z, &
						 scalarField%scalar, delTri%list, delTri%lptr, delTri%lend, SIGMA_FLAG, self%sigma1, GRAD_FLAG, &
						 self%grad1, startTriangle, interpOut(i,j), errCode)
		enddo
	enddo
	
	do i = 0, numProcs - 1
		call MPI_BCAST( interpOut(:,mpiLons%indexStart(i):mpiLons%indexEnd(i)), nLat * mpiLons%messageLength(i), &
				MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, errCode)
	enddo
	
	call Delete(mpiLons)
end subroutine

function InterpolateVector( lon, lat, self, aMesh, delTri, vectorField)
	real(kreal), dimension(3) :: InterpolateVector
	type(SSRFPACKInterface), intent(in) :: self
	type(PolyMesh2d), intent(in) :: aMesh
	type(DelaunayTriangulation), intent(in) :: delTri
	type(Field), intent(in) :: vectorField
	real(kreal), intent(in) :: lon
	real(kreal), intent(in) :: lat
	integer(kint) :: errCode
	
	call INTRC1( aMesh%particles%n, lat, lon, aMesh%particles%x, aMesh%particles%y, aMesh%particles%z, &
				 vectorField%xComp, delTri%list, delTri%lptr, delTri%lend, SIGMA_FLAG, self%sigma1, GRAD_FLAG, &
				 self%grad1, startTriangle, InterpolateVector(1), errCode)
	call INTRC1( aMesh%particles%n, lat, lon, aMesh%particles%x, aMesh%particles%y, aMesh%particles%z, &
				 vectorField%yComp, delTri%list, delTri%lptr, delTri%lend, SIGMA_FLAG, self%sigma2, GRAD_FLAG, &
				 self%grad2, startTriangle, InterpolateVector(2), errCode)
	call INTRC1( aMesh%particles%n, lat, lon, aMesh%particles%x, aMesh%particles%y, aMesh%particles%z, &
				 vectorField%zComp, delTri%list, delTri%lptr, delTri%lend, SIGMA_FLAG, self%sigma3, GRAD_FLAG, &
				 self%grad3, startTriangle, InterpolateVector(3), errCode)				 				 
end function

subroutine InterpolateVectorToUnifLatLonGrid( self, aMesh, delTri, vectorField, lons, lats, interpX, interpY, interpZ)
	type(SSRFPACKInterface), intent(inout) :: self
	type(PolyMesh2d), intent(in) :: aMesh
	type(DelaunayTriangulation), intent(in) :: delTri
	type(Field), intent(in) :: vectorField
	real(kreal), dimension(:), intent(in) :: lons
	real(kreal), dimension(:), intent(in) :: lats
	real(kreal), dimension(:,:), intent(out) :: interpX
	real(kreal), dimension(:,:), intent(out) :: interpY
	real(kreal), dimension(:,:), intent(out) :: interpZ
	!
	integer(kint) :: i, j, errCode, nLat, nLon
	type(MPISetup) :: mpiLons
	
	nLat = size(lats)
	nLon = size(lons)
	call New(mpiLons, nLon, numProcs)
	
	do j = mpiLons%indexStart(procRank), mpiLons%indexEnd(procRank)
		do i = 1, nLat
			call INTRC1( aMesh%particles%n, lats(i), lons(j), aMesh%particles%x, aMesh%particles%y, aMesh%particles%z, &
					 vectorField%xComp, delTri%list, delTri%lptr, delTri%lend, SIGMA_FLAG, self%sigma1, GRAD_FLAG, &
					 self%grad1, startTriangle, interpX(i,j), errCode)
			call INTRC1( aMesh%particles%n, lats(i), lons(j), aMesh%particles%x, aMesh%particles%y, aMesh%particles%z, &
					 vectorField%yComp, delTri%list, delTri%lptr, delTri%lend, SIGMA_FLAG, self%sigma2, GRAD_FLAG, &
					 self%grad2, startTriangle, interpY(i,j), errCode)
			call INTRC1( aMesh%particles%n, lats(i), lons(j), aMesh%particles%x, aMesh%particles%y, aMesh%particles%z, &
					 vectorField%zComp, delTri%list, delTri%lptr, delTri%lend, SIGMA_FLAG, self%sigma3, GRAD_FLAG, &
					 self%grad3, startTriangle, interpZ(i,j), errCode)	
		enddo
	enddo
	
	do i = 0, numProcs - 1
		call MPI_BCAST( interpX(:,mpiLons%indexstart(i):mpiLons%indexEnd(i)), nLat * mpiLons%messageLength(i), &
				MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, errCode)
		call MPI_BCAST( interpY(:,mpiLons%indexstart(i):mpiLons%indexEnd(i)), nLat * mpiLons%messageLength(i), &
				MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, errCode)
		call MPI_BCAST( interpZ(:,mpiLons%indexstart(i):mpiLons%indexEnd(i)), nLat * mpiLons%messageLength(i), &
				MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, errCode)								
	enddo
	
	call Delete(mpiLons)
end subroutine

function InterpolateLagParam(lon, lat, self, aMesh, delTri )
	real(kreal), dimension(3) :: InterpolateLagParam
	type(SSRFPACKInterface), intent(in) :: self
	type(PolyMesh2d), intent(in) :: aMesh
	type(DelaunayTriangulation), intent(in) :: delTri
	real(kreal), intent(in) :: lon
	real(kreal), intent(in) :: lat
	integer(kint) :: errCode
	
	call INTRC1( aMesh%particles%N, lat, lon, aMesh%particles%x, aMesh%particles%y, aMesh%particles%z, &
				 aMesh%particles%x0, delTri%list, delTri%lptr, delTri%lend, SIGMA_FLAG, self%sigma1, GRAD_FLAG, &
				 self%grad1, startTriangle, InterpolateLagParam(1), errCode)
	call INTRC1( aMesh%particles%N, lat, lon, aMesh%particles%x, aMesh%particles%y, aMesh%particles%z, &
				 aMesh%particles%y0, delTri%list, delTri%lptr, delTri%lend, SIGMA_FLAG, self%sigma2, GRAD_FLAG, &
				 self%grad2, startTriangle, InterpolateLagParam(2), errCode)
	call INTRC1( aMesh%particles%N, lat, lon, aMesh%particles%x, aMesh%particles%y, aMesh%particles%z, &
				 aMesh%particles%z0, delTri%list, delTri%lptr, delTri%lend, SIGMA_FLAG, self%sigma3, GRAD_FLAG, &
				 self%grad3, startTriangle, InterpolateLagParam(3), errCode)				 
end function

subroutine SetSigmaFlag( newVal )
	integer(kint), intent(in) :: newVal
	SIGMA_FLAG = newVal
end subroutine

!
!----------------
! private methods
!----------------
!
subroutine BuildDelaunayTriangulation(self, aMesh)
	type(DelaunayTriangulation), intent(inout) :: self
	type(PolyMesh2d), intent(in) :: aMesh
	!
	real(kreal), allocatable, dimension(:) :: dist
	integer(kint), allocatable, dimension(:) :: near
	integer(kint), allocatable, dimension(:) :: next
	integer(kint) :: lnew, errCode
	
	allocate(dist(amesh%particles%n))
	allocate(near(amesh%particles%n))
	allocate(next(amesh%particles%n))
	
	call TRMESH(aMesh%particles%N, aMesh%particles%x, aMesh%particles%y, aMesh%particles%z, &
				self%list, self%lptr, self%lend, lnew, near, next, dist, errCode )

	if ( errCode == -1 ) then
		call LogMessage(log,ERROR_LOGGING_LEVEL, trim(logKey)//' TRMESH ERROR :',' found n < 3 points.')
	elseif (errCode == -2 ) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,trim(logKey)//' TRMESH ERROR :',' found first three nodes to be colinear.')
	elseif ( errCode > 0 ) then
		write(logString,'(A,I8,A)') ' node ',errCode,' is a duplicate.'
		call LogMessage(log,ERROR_LOGGING_LEVEL,trim(logKey)//' TRMESH ERROR :', trim(logString))
	endif				
	
	deallocate(dist)
	deallocate(near)
	deallocate(next)
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
