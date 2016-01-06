module SSRFPACKInterfaceModule
!> @file ssrfpackInterface.f90
!> @author Peter Bosler, Sandia National Laboratories, Center for Computing Research
!> 
!> @defgroup ssrfpackInterface  ssrfpackInterface
!> @brief Interface and workspace for interpolation of LPM data in spherical domains using the STRIPACK and SSRFPACK libraries.
!> 
!> STRIPACK produces a Delaunay triangulation of scattered data points on the surface of a unit sphere. @n
!> SSRFPACK performs interpolation on that triangulation using a linear combination of cubic Hermite splines with 
!> exponential tension factors.  
!> Derivatives for the Hermite interpolating polynomials are estimated using a least-squares projection onto local quadratic polynomials.
!> 
!> For more information, see
!> 1. R. Renka, Algorithm 772: STRIPACK: Delaunay triangulation and Voronoi diagram on the surface of a sphere, _ACM TOMS_ 23, 1997. (stripack.f)
!> 2. R. Renka, Algorithm 773: SSRFPACK: Interpolation scattered data on the surface of a sphere with a surface under tension, _ACM TOMS_ 23, 1997. (ssrfpack.f)
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
	integer(kint), dimension(:), allocatable :: list !< Stripack data structure, see stripack.f for more details
	integer(kint), dimension(:), allocatable :: lptr !< Stripack data structure, see stripack.f for more details
	integer(kint), dimension(:), allocatable :: lend !< Stripack data structure, see stripack.f for more details
	
	contains
	
		final :: deleteDelTri
end type

type SSRFPACKInterface
	real(kreal), dimension(:,:), allocatable :: grad1 !< Storage for estimated gradient vectors
	real(kreal), dimension(:,:), allocatable :: grad2 !< Storage for estimated gradient vectors
	real(kreal), dimension(:,:), allocatable :: grad3 !< Storage for estimated gradient vectors
	real(kreal), dimension(:), allocatable :: sigma1 !< Storage for smoothing factors
	real(kreal), dimension(:), allocatable :: sigma2 !< Storage for smoothing factors
	real(kreal), dimension(:), allocatable :: sigma3 !< Storage for smoothing factors
	
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

!> @brief Allocates memory for an SsrfpackInterface object.
!> 
!> @param[out] self Target SSRFPACKInterface
!> @param[inout] aMesh @ref PolyMesh2d spherical mesh
!> @param[in] nDim Dimension of data to be interpolated
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

!> @brief Allocates memory for STRIPACK's Delaunay triangulation data structure, constructs Delaunay triangulation.
!> 
!> @param[out] self Target Delaunay triangulation of all particles (both centers and vertices)
!> @param[in] aMesh @ref PolyMesh2d spherical mesh
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

!> @brief Deletes and frees memory associated with a Delaunay triangulation
!> @param[inout] self Target Del. tri.
subroutine deleteDelTri(self)
	type(DelaunayTriangulation), intent(inout) :: self
	if ( allocated(self%list)) then
		deallocate(self%list)
		deallocate(self%lptr)
		deallocate(self%lend)
	endif
end subroutine

!> @brief Deletes and frees memory associated with an SSRFPACKInterface object.
!> @param[inout] Target SSRFPACKInterface
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

!> @brief Sets up an SSRFPACKInterface object for interpolation.  
!>
!> * Estimates scalar gradients of a scalar field, or scalar gradients of each component of a vector field
!> * Computes optimal smoothing factors, see ssrfpack.f for more detail
!> 
!> @param[inout] self target SSRFPACKInterface, on output, target is ready to interpolate data
!> @param[in] aMesh @ref PolyMesh2d
!> @param[in] delTri Delaunay triangulation
!> @param[in] vectorField vector @ref Field associated with aMesh
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

!> @brief Sets up an SSRFPACKInterface object for interpolation.  
!>
!> * Estimates scalar gradients of a scalar field, or scalar gradients of each component of a vector field
!> * Computes optimal smoothing factors, see ssrfpack.f for more detail
!> 
!> @param[inout] self target SSRFPACKInterface, on output, target is ready to interpolate data
!> @param[in] aMesh @ref PolyMesh2d
!> @param[in] delTri Delaunay triangulation
!> @param[in] scalarField scalar @ref Field associated with aMesh
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

!> @brief Sets up an SSRFPACKInterface object for interpolation of the Lagrangian parameter.  
!>
!> * Estimates scalar gradients of each component of the Lagrangian coordinate vector at each particle
!> * Computes optimal smoothing factors, see ssrfpack.f for more detail
!> 
!> @param[inout] self target SSRFPACKInterface, on output, target is ready to interpolate data
!> @param[in] aMesh @ref PolyMesh2d
!> @param[in] delTri Delaunay triangulation
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

!> @brief Performs interpolation of a scalar field at a single point.
!>
!> * Locates (lat, lon) in Delaunay triangulation
!> * Computes interpolated value using cubic Hermite polynomial on triangle containing (lat,lon) 
!>
!> @param[in] lon Longitude of interpolation output 
!> @param[in] lat Latitude of interpolation output
!> @param[in] self SSRFPACKInterface ready for interpolation (a SetSource* subroutine has already been called)
!> @param[in] aMesh @ref PolyMesh2d spherical mesh
!> @param[in] delTri Delaunay triangulation
!> @param[in] scalarField source data
!> @return interpolated scalar value
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

!> @brief Performs parallel interpolation of a scalar field from an LPM particle set to a uniform latitude-longitude grid.
!> 
!> Distributes work of interpolation across MPI processes by assigning each process a subset of the output points (divided by longitude).
!> Then broadcasts each processes work to the other processes.
!> 
!> @param[inout] self SSRFPACKInterface
!> @param[in] aMesh @ref PolyMesh2d source mesh
!> @param[in] delTri Delaunay triangulation
!> @param[in] scalarField source data
!> @param[in] lons vector of longitudes associated with uniform lat-lon output grid
!> @param[in] lats vector of latitudes associated with uniform lat-lon output grid
!> @param[out] interpOut interpolated scalar output
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

!> @brief Performs interpolation of a vector field at a single point.
!>
!> * Locates (lat, lon) in Delaunay triangulation
!> * Computes interpolated value using cubic Hermite polynomial on triangle containing (lat,lon) 
!>
!> @param[in] lon Longitude of interpolation output 
!> @param[in] lat Latitude of interpolation output
!> @param[in] self SSRFPACKInterface ready for interpolation (a SetSource* subroutine has already been called)
!> @param[in] aMesh @ref PolyMesh2d spherical mesh
!> @param[in] delTri Delaunay triangulation
!> @param[in] vectorField source data
!> @return interpolated vector value
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

!> @brief Performs parallel interpolation of a vector field from an LPM particle set to a uniform latitude-longitude grid.
!> 
!> Distributes work of interpolation across MPI processes by assigning each process a subset of the output points (divided by longitude).
!> Then broadcasts each processes work to the other processes.
!> 
!> @param[inout] self SSRFPACKInterface ready for interpolation (a SetSource* subroutine has already been called)
!> @param[in] aMesh @ref PolyMesh2d source mesh
!> @param[in] delTri Delaunay triangulation
!> @param[in] vectorField source data
!> @param[in] lons vector of longitudes associated with uniform lat-lon output grid
!> @param[in] lats vector of latitudes associated with uniform lat-lon output grid
!> @param[out] interpX interpolated x component of vector field
!> @param[out] interpY interpolated y component of vector field
!> @param[out] interpZ interpolated z component of vector field
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

!> @brief Interpolates the Lagrangian parameter from a set of LPM particles to a desired location on the sphere.
!> 
!> @param[in] lon longitude of interpolation output
!> @param[in] lat latitude of interpolation output
!> @param[inout] self SSRFPACKInterface ready for interpolation (a SetSourceLagrangianParameter subroutine has already been called)
!> @param[in] aMesh source @ref PolyMesh2d
!> @param[in] delTri Delaunay triangulation
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

!> @brief Builds a Delaunay triangulation of a set of LPM @ref Particles on the surface of the unit sphere using stripack.f.
!> 
!> @param[inout] self STRIPACK data structures for Delaunay triangulation
!> @param[in] aMesh spherical @ref PolyMesh2d
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

!> @brief Initializes a logger for the SSRFPACKInterface module
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
