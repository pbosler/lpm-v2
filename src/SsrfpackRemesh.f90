module SSRFPACKRemeshModule
!> @file SsrfpackRemesh.f90
!> Data structure and methods for remapping spherical LPM data using SSRFPACK.
!> @author Peter Bosler, Sandia National Laboratories, Center for Computing Research
!>
!> @defgroup SSRFPACKRemesh SSRFPACKRemesh
!> Data structure and methods for remapping spherical LPM data using SSRFPACK. @n
!>
!> Uses the interfaces provided by @ref ssrfpackInterface module.
!> For references describing the SSRFPACK package and its counterpart STRIPACK, see the @ref ssrfpackInterface module's detailed description.
!>
!> @{
use NumberKindsModule
use STDIntVectorModule
use OutputWriterModule
use LoggerModule
use ParticlesModule
use FieldModule
use EdgesModule
use FacesModule
use PlaneGeomModule
use SphereGeomModule
use PolyMesh2dModule
use SphereBVEModule
use SphereTransportModule
use RefinementModule
use MPISetupModule
use SSRFPACKInterfaceModule

implicit none

include 'mpif.h'

private
public BVERemeshSource, TransportRemesh, New, Delete
public DirectRemeshBVE, DirectRemeshTransport
public LagrangianRemeshBVEWithVorticityFunction
public LagrangianRemeshBVEToReferenceMesh
public LagrangianRemeshTransportWithFunctions
public FTLECalc

!----------------
! types and module variables
!----------------
!
type BVERemeshSource
	type(DelaunayTriangulation) :: delTri  	!< Delaunay triangulation from STRIPACK
	type(SSRFPACKInterface) :: relVortSource !< Source data for relative vorticity interpolation
	type(SSRFPACKInterface) :: absVortSource !< Source data for absolute vorticity interpolation
	type(SSRFPACKInterface) :: lagParamSource !< Source data for Lagrangian parameter interpolation
	type(SSRFPACKInterface), dimension(:), allocatable :: tracerSource !< Source data for passive tracers

	contains
		final :: deleteBVE
end type

type TransportRemesh
	type(DelaunayTriangulation) :: delTri
	type(SSRFPACKInterface) :: densitySource
	type(SSRFPACKInterface) :: lagParamSource
	type(SSRFPACKInterface), dimension(:), allocatable :: tracerSource

	contains
		final :: deleteTransport
end type

integer(kint), save :: remeshCounter = 0

!
!----------------
! interfaces
!----------------
!
interface New
	module procedure newBVE
	module procedure newTransport
end interface

interface Delete
	module procedure deleteBVE
	module procedure deleteTransport
end interface

!
!----------------
! Logging
!----------------
!
logical(klog), save :: logInit = .FALSE.
type(Logger) :: log
character(len=28), save :: logKey = 'SSRFRemesh'
character(len=MAX_STRING_LENGTH) :: logString
integer(kint), parameter :: logLevel = DEBUG_LOGGING_LEVEL

contains
!
!----------------
! public methods
!----------------
!

!> @brief Allocates memory and initializes a remapping utility for a @ref SphereBVE BVE mesh.
!>
!> @param[out] self target Remeshing data structure
!> @param[in] oldSphere source @ref SphereBVE
subroutine newBVE( self, oldSphere )
	type(BVERemeshSource), intent(out) :: self
	type(BVEMesh), intent(inout) :: oldSphere
	!
	integer(kint) :: i

	if ( .NOT. logInit ) call InitLogger(log, procRank)

	call LogMessage(log, DEBUG_LOGGING_LEVEL, trim(logKey)//" newBVERemesh : ", "entering.")

	call New(self%delTri, oldsphere%mesh)

	call New(self%lagParamSource, oldSphere%mesh, 3)
	call SetSourceLagrangianParameter(self%lagParamSource, oldSphere%mesh, self%delTri)

	call New(self%relVortSource, oldsphere%mesh, 1)
	call SetScalarSourceData(self%relVortSource, oldSphere%mesh, self%delTri, oldSphere%relVort)

	call New(self%absVortSource, oldsphere%mesh, 1)
	call SetScalarSourceData(self%absVortSource, oldSphere%mesh, self%delTri, oldSphere%absVort)

	if ( allocated(oldSphere%tracers)) then
!		call LogMessage(log, DEBUG_LOGGING_LEVEL, trim(logKey)//" newBVERemesh : found nTracers = ", size(oldSphere%tracers))
!		call StartSection( log, "newBVERemesh sees this sphere")
!		call LogStats(oldSphere,log)
!		call EndSection(log)

		allocate(self%tracerSource(size(oldSphere%tracers)))

		do i = 1, size(oldSphere%tracers)
			call New(self%tracerSource(i), oldsphere%mesh, oldSphere%tracers(i)%nDim)
			if ( oldSphere%tracers(i)%nDim == 1 ) then
				call SetScalarSourceData(self%tracerSource(i), oldSphere%mesh, self%delTri, oldSphere%tracers(i))
			else
				call SetVectorSourceData(self%tracerSource(i), oldSphere%mesh, self%delTri, oldSphere%tracers(i))
			endif
		enddo
	endif
	call LogMessage(log, DEBUG_LOGGING_LEVEL, trim(logKey)//" newBVERemesh : ", "returning.")
end subroutine

subroutine newTransport(self, oldSphere)
	type(TransportRemesh), intent(out) :: self
	type(TransportMesh), intent(inout) :: oldSphere
	!
	integer(kint) :: i

	if ( .NOT. logInit) call InitLogger(log, procRank)

	call New(self%delTri, oldSphere%mesh)
	call New(self%lagParamSource, oldSphere%mesh, 3)
	call New(self%densitySource, oldSphere%mesh, 1)

	call SetScalarSourceData(self%densitySource, oldSphere%mesh, self%delTri, oldSphere%density)

	if ( allocated(oldSphere%tracers) ) then
		allocate(self%tracerSource(size(oldSphere%tracers)))

		do i = 1, size(oldSphere%tracers)
			call New(self%tracerSource(i), oldSphere%mesh, oldSphere%tracers(i)%nDim)
			if ( oldSphere%tracers(i)%nDim == 1 ) then
				call SetScalarSourceData(self%tracerSource(i), oldSphere%mesh, self%delTri, oldSphere%tracers(i))
			else
				call SetVectorSourceData(self%tracerSource(i), oldSphere%mesh, self%delTri, oldSphere%tracers(i))
			endif
		enddo
	endif

	call SetSourceLagrangianParameter(self%lagParamSource, oldSphere%mesh, self%delTri)
end subroutine

subroutine deleteTransport(self)
	type(TransportRemesh), intent(inout) :: self
	integer(kint) :: i

	call Delete(self%delTri)
	call Delete(self%lagParamSource)
	call Delete(self%densitySource)
	if ( allocated(self%tracerSource) ) then
		do i = 1, size(self%tracerSource)
			call Delete(self%tracerSource(i))
		enddo
		deallocate(self%tracerSource)
	endif
end subroutine


!> @brief Deletes and frees memory associated with a BVE remeshing object
!> @param[inout] self Target remeshing object
subroutine deleteBVE( self )
	type(BVERemeshSource), intent(inout) :: self
	!
	integer(kint) :: i

	call Delete(self%delTri)
	call Delete(self%lagParamSource)
	call Delete(self%relVortSource)
	call Delete(self%absVortSource)
	if ( allocated(self%tracerSource) ) then
		do i = 1, size(self%tracerSource)
			call Delete(self%tracerSource(i))
		enddo
		deallocate(self%tracerSource)
	endif
end subroutine

!> @brief Performs a remesh/remap of an LPM simulation of the barotropic vorticity equation on the sphere using
!> direct interpolation of each variable.
!>
!> @param[in] self Remeshing data structure
!> @param[in] oldSphere source @ref SphereBVE mesh
!> @param[inout] newSphere target @ref SphereBVE mesh (note that this must have been allocated prior to calling this subroutine)
!> @param[in] AMR .TRUE. if adaptive refinement will be used
!> @param[in] vortFlagFn1 FlagFunction for vorticity refinement, must have same interface as refinementmodule::FlagFunction
!> @param[in] tol1 tolerance for vortFlagFn1
!> @param[in] desc1 description of AMR criterion used for vortFlagFn1
!> @param[in] flagFn2 Flag function for refinement of field2
!> @param[in] tol2 tolerance for flagFn2
!> @param[in] desc2 description of AMR criterion used for flagFn2
!> @param[inout] field2 second field to use for adaptive refinement
subroutine DirectRemeshBVE(self, oldSphere, newSphere, AMR, vortFlagFn1, tol1, desc1, flagFn2, tol2, desc2, field2 )
	type(BVERemeshSource), intent(in) :: self
	type(BVEMesh), intent(in) :: oldSphere
	type(BVEMesh), intent(inout) :: newSphere
	logical(klog), intent(in) :: AMR
	procedure(FlagFunction), optional :: vortFlagFn1
	real(kreal), intent(in), optional :: tol1
	character(len=*), intent(in), optional :: desc1
	procedure(FlagFunction), optional :: flagFn2
	real(kreal), intent(in), optional :: tol2
	character(len=*), intent(in), optional :: desc2
	type(Field), intent(inout), optional :: field2
	!
	integer(kint) :: i, j, k
	real(kreal) :: lon, lat
	real(kreal), dimension(3) :: x0, vecT
	type(RefineSetup) :: refine
	integer(kint) :: refineVariableCount
	integer(kint) :: nParticlesBefore, nParticlesAfter

	remeshCounter = remeshCounter + 1
	!
	!	Remesh to new uniform mesh
	!
	newSphere%relVort%N = newSphere%mesh%particles%N
	newSphere%absVort%N = newSphere%mesh%particles%N
	do i = 1, newSphere%mesh%particles%N
		lon = Longitude( newSphere%mesh%particles%x(i), newSphere%mesh%particles%y(i), newSphere%mesh%particles%z(i))
		lat = Latitude( newSphere%mesh%particles%x(i), newSphere%mesh%particles%y(i), newSphere%mesh%particles%z(i))
		newSphere%relVort%scalar(i)	= InterpolateScalar( lon, lat, self%relVortSource, oldSphere%mesh, &
														 self%delTri, oldSphere%relVort)
		newSphere%absVort%scalar(i)	= InterpolateScalar( lon, lat, self%absVortSource, oldSphere%mesh, &
														 self%delTri, oldSphere%absVort)
		x0 = InterpolateLagParam( lon, lat, self%lagParamSource, oldSphere%mesh, self%delTri)
		newSphere%mesh%particles%x0(i) = x0(1)
		newSphere%mesh%particles%y0(i) = x0(2)
		newSphere%mesh%particles%z0(i) = x0(3)
	enddo

	if (allocated(newSphere%tracers)) then
		do i = 1, size(newSphere%tracers)
			newSphere%tracers(i)%N = newSphere%mesh%particles%N
		enddo
		do j = 1, newSphere%mesh%particles%N
			lon = Longitude( newSphere%mesh%particles%x(j), newSphere%mesh%particles%y(j), newSphere%mesh%particles%z(j))
			lat = Latitude( newSphere%mesh%particles%x(j), newSphere%mesh%particles%y(j), newSphere%mesh%particles%z(j))
			do i = 1, size(newSphere%tracers)
				 if ( newSphere%tracers(i)%nDim == 1 ) then
				 	newSphere%tracers(i)%scalar(j) = InterpolateScalar(lon, lat, self%tracerSource(i), oldSphere%mesh, &
				 													   self%delTri, oldSphere%tracers(i))
				 else
				 	vecT = InterpolateVector(lon, lat, self%tracerSource(i), oldSphere%mesh, self%delTri, oldSphere%tracers(i))
				 	newSphere%tracers(i)%xComp(j) = vecT(1)
				 	newSphere%tracers(i)%yComp(j) = vecT(2)
				 	newSphere%tracers(i)%zComp(j) = vecT(3)
				 endif
			enddo
		enddo
	endif


	!
	!  Perform adaptive refinement, interpolate Field data to new Particles
	!
	if ( AMR ) then
		refineVariableCount = 0
		if ( present(vortFlagFn1) .AND. (present(tol1) .AND. present(desc1)) ) then
			refineVariableCount = 1
			if ( present(flagFn2) .AND. ( present(tol2) .AND. present(desc2) ) ) then
				refineVariableCount = 2
			endif
		endif

		if (refineVariableCount > 0 ) then
			call New(refine, newSphere%mesh%faces%N_Max)

			do i = 1, newSphere%mesh%amrLimit
				nParticlesBefore = newSphere%mesh%particles%N
				if ( refineVariableCount == 1 ) then
					call IterateMeshRefinementOneVariable( refine, newSphere%mesh, newSphere%relVort, &
											vortFlagFn1, tol1, desc1, nParticlesBefore, nParticlesAfter)
				elseif ( refineVariableCount == 2 ) then
					if ( present(field2) ) then
						call IterateMeshRefinementTwoVariables(refine, newSphere%mesh, newSphere%relVort, vortFlagFn1, &
								tol1, desc1, field2, flagFn2, tol2, desc2, nParticlesBefore, nParticlesAfter)
					else
						call IterateMeshRefinementTwoVariables(refine, newSphere%mesh, newSphere%relVort, vortFlagFn1, &
								tol1, desc1, newSphere%relVort, flagFn2, tol2, desc2, nParticlesBefore, nParticlesAfter)
					endif
				endif

				do j = nParticlesBefore + 1, nParticlesAfter
					lon = Longitude( newSphere%mesh%particles%x(j), newSphere%mesh%particles%y(j), &
									 newSphere%mesh%particles%z(j))
					lat = Latitude( newSphere%mesh%particles%x(j), newSphere%mesh%particles%y(j), &
									newSphere%mesh%particles%z(j))

					x0 = InterpolateLagParam( lon, lat, self%lagParamSource, oldSphere%mesh, self%delTri)
					newSphere%mesh%particles%x0(j) = x0(1)
					newSphere%mesh%particles%y0(j) = x0(2)
					newSphere%mesh%particles%z0(j) = x0(3)

					newSphere%relVort%scalar(j) = InterpolateScalar( lon, lat, self%relVortSource, oldSphere%mesh, &
																	 self%delTri, oldSphere%relVort)
					newSphere%absVort%scalar(j) = InterpolateScalar( lon, lat, self%absVortSource, oldSphere%mesh, &
																	 self%delTri, oldSphere%absVort)
					if (allocated(newSphere%tracers)) then
						do k = 1, size(newSphere%tracers)
							if ( newSphere%tracers(k)%nDim == 1 ) then
								newSphere%tracers(k)%scalar(j) = InterpolateScalar( lon, lat, self%tracerSource(k), &
													oldSphere%mesh, self%delTri, oldSphere%tracers(k) )
							else
								vecT = InterpolateVector(lon, lat, self%tracerSource(k), oldSphere%mesh, self%delTri, &
													oldSphere%tracers(k))
								newSphere%tracers(k)%xComp(j) = vecT(1)
								newSphere%tracers(k)%yComp(j) = vecT(2)
								newSphere%tracers(k)%zComp(j) = vecT(3)
							endif
						enddo
					endif
				enddo
			enddo
			newSphere%relVort%N = newSphere%mesh%particles%N
			newSphere%absVort%N = newSphere%mesh%particles%N
			if ( allocated(newSphere%tracers) ) then
				do k = 1, size(newSphere%tracers)
					newSphere%tracers(k)%N = newSphere%mesh%particles%N
				enddo
			endif
			call LoadBalance(newSphere%mpiParticles, newSphere%mesh%particles%N, numProcs)
			call Delete(refine)
		endif
	endif

	!
	!	set velocity and stream functions
	!
	call SetVelocityOnMesh( newSphere )
	call SetStreamFunctionsOnMesh( newSphere )


end subroutine

!> @brief Performs a remesh/remap of an LPM BVE simulation using indirect interpolation for the variables in a @ref SphereBVE mesh.
!> Remaps to reference time t = 0.
!>
!> @param[in] self Remeshing data structure
!> @param[in] oldSphere source @ref SphereBVE mesh
!> @param[inout] newSphere target @ref SphereBVE mesh (note that this must have been allocated prior to calling this subroutine)
!> @param[in] AMR true if adaptive refinement will be used
!> @param[in] relVortFn Vorticity distribution function, must have same interface as numberkindsmodule::scalarFnOf3DSpace
!> @param[in] flagFn1 FlagFunction for vorticity refinement, must have same interface as refinementmodule::FlagFunction
!> @param[in] tol1 tolerance for flagFn1
!> @param[in] desc1 description of AMR criterion used for flagFn1
!> @param[in] flagFn2 FlagFunction for vorticity refinement, must have same interface as refinementmodule::FlagFunction
!> @param[in] tol2 tolerance for flagFn2
!> @param[in] desc2 description of AMR criterion used for flagFn2
!> @param[in] RefineFLowMapYN True if refinement of the flow map will be used
!> @param[in] flowMapVarTol tolerance value for Lagrangian coordinate variation per face
!> @param[in] nLagTracers number of Lagrangian passive tracers (currently 0, 1, or 2 are only values allowed)
!> @param[in] tracerFn1 tracer distribution function, must have same interface as numberkindsmodule::scalarFnOf3DSpace
!> @param[in] tracerFn2 tracer distribution function, must have same interface as numberkindsmodule::scalarFnOf3DSpace
subroutine LagrangianRemeshBVEWithVorticityFunction( self, oldSphere, newSphere, AMR, relVortFn, flagFn1, tol1, desc1, &
								flagFn2, tol2, desc2, RefineFlowMapYN, flowMapVarTol, nLagTracers, tracerFn1, tracerFn2 )
	type(BVERemeshSource), intent(in) :: self
	type(BVEMesh), intent(in) :: oldSphere
	type(BVEMesh), intent(inout) :: newSphere
	logical(klog), intent(in) :: AMR
	procedure(scalarFnOf3DSpace) :: relVortFn
	procedure(FlagFunction) :: flagFn1
	real(kreal), intent(in) :: tol1
	character(len=*), intent(in) :: desc1
	procedure(FlagFunction), optional :: flagFn2
	real(kreal), intent(in), optional :: tol2
	character(len=*), intent(in), optional :: desc2
	logical(klog), intent(in), optional :: RefineFlowMapYN
	real(kreal), intent(in), optional :: flowMapVarTol
	integer(kint), intent(in), optional :: nLagTracers
	procedure(scalarFnOf3DSpace), optional :: tracerFn1
	procedure(scalarFnOf3DSpace), optional :: tracerFn2
	!
	integer(kint) :: nTracers
	logical(klog) :: refineVorticityTwice, doFlowMapRefinement
	integer(kint) :: i, j, k
	real(kreal) :: lon, lat
	real(kreal), dimension(3) :: x0
	type(RefineSetup) :: refine
	integer(kint) :: nParticlesBefore, nParticlesAfter
	real(kreal), dimension(3) :: vecT

	remeshCounter = remeshCounter + 1

	nTracers = 0
	if ( present(tracerFn1) .AND. ( nLagTracers >= 1 .AND. newSphere%tracers(1)%nDim == 1) ) then
		nTracers = 1
		if ( present(tracerFn2) .AND. ( nLagTracers == 2 .AND. newSphere%tracers(2)%nDim == 1) ) nTracers = 2
	endif

	newSphere%relVort%N = newSphere%mesh%particles%N
	newSphere%absVort%N = newSphere%mesh%particles%N
	if ( allocated(newSphere%tracers) ) then
		do i = 1, size(newSphere%tracers)
			newSphere%tracers(i)%N = newSphere%mesh%particles%N
		enddo
	endif

	do i = 1, newSphere%mesh%particles%N
		lon = Longitude( newSphere%mesh%particles%x(i), newSphere%mesh%particles%y(i), newSphere%mesh%particles%z(i))
		lat = Latitude( newSphere%mesh%particles%x(i), newSphere%mesh%particles%y(i), newSphere%mesh%particles%z(i))
		x0 = InterpolateLagParam( lon, lat, self%lagParamSource, oldSphere%mesh, self%delTri)
		newSphere%mesh%particles%x0(i) = x0(1)
		newSphere%mesh%particles%y0(i) = x0(2)
		newSphere%mesh%particles%z0(i) = x0(3)

		newSphere%absVort%scalar(i) = relVortFn( x0(1), x0(2), x0(3)) + &
				2.0_kreal * newSphere%rotationRate * x0(3) / newSphere%radius

		newSphere%relVort%scalar(i) = newSphere%absVort%scalar(i) - 2.0_kreal * &
			newSphere%rotationRate * newSphere%mesh%particles%z(i) / newSphere%radius

		if ( allocated(newSphere%tracers) ) then
			if ( nTracers == 1 ) then
				newSphere%tracers(1)%scalar(i) = tracerFn1( x0(1), x0(2), x0(3) )
			elseif ( nTracers == 2 ) then
				newSphere%tracers(1)%scalar(i) = tracerFn1( x0(1), x0(2), x0(3) )
				newSphere%tracers(2)%scalar(i) = tracerFn2( x0(1), x0(2), x0(3) )
			endif
			do j = nTracers+1, size(newSphere%tracers)
				if ( newSphere%tracers(j)%nDim == 1 ) then
					newSphere%tracers(j)%scalar(i) = InterpolateScalar(lon, lat, self%tracerSource(j), oldSphere%mesh, &
																	   self%delTri, oldSphere%tracers(j) )
				else
					vecT = InterpolateVector(lon, lat, self%tracerSource(j), oldSphere%mesh, &
											 self%delTri, oldSphere%tracers(j) )
					newSphere%tracers(j)%xComp(i) = vecT(1)
					newSphere%tracers(j)%yComp(i) = vecT(2)
					newSphere%tracers(j)%zComp(i) = vecT(3)
				endif
			enddo
		endif
	enddo

	!
	!	AMR
	!
	if ( AMR ) then
		refineVorticityTwice = ( present(flagFn2) .AND. ( present(tol2) .AND. present(desc2) ) )
		doFlowMapRefinement = ( RefineFlowMapYN .AND. present(flowMapVarTol))

		call New(refine, newSphere%mesh%faces%N_Max)

		do i = 1, newSphere%mesh%amrLimit
			nParticlesBefore = newSphere%mesh%particles%N

			if ( refineVorticityTwice ) then
				if ( doFlowMapRefinement ) then
					call IterateMeshRefinementTwoVariablesAndFlowMap(refine, newSphere%mesh, newSphere%relVort, &
							flagFn1, tol1, desc1, newSphere%relVort, flagFn2, tol2, desc2, &
							flowMapVarTol, nParticlesBefore, nParticlesAfter)
				else
					call IterateMeshRefinementTwoVariables(refine, newSphere%mesh, newSphere%relVort, &
						flagFn1, tol1, desc1, newSphere%relVort, flagFn2, tol2, desc2, nParticlesBefore, nParticlesAfter)
				endif
			else
				if ( doFlowMapRefinement ) then
					call IterateMeshRefinementOneVariableAndFlowMap(refine, newSphere%mesh, newSphere%relVort, flagFn1, &
										tol1, desc1, flowMapVarTol, nParticlesBefore, nParticlesAfter)
				else
					call IterateMeshRefinementOneVariable(refine, newSphere%mesh, newSphere%relVort, flagFn1, tol1, desc1, &
							nParticlesBefore, nParticlesAfter)
				endif
			endif

			do j = nParticlesBefore + 1, nParticlesAfter
				lon = Longitude( newSphere%mesh%particles%x(j), newSphere%mesh%particles%y(j), &
									 newSphere%mesh%particles%z(j))
				lat = Latitude( newSphere%mesh%particles%x(j), newSphere%mesh%particles%y(j), &
									newSphere%mesh%particles%z(j))
				x0 = InterpolateLagParam(lon, lat, self%lagParamSource, oldSphere%mesh, self%delTri )
				newSphere%mesh%particles%x0(j) = x0(1)
				newSphere%mesh%particles%y0(j) = x0(2)
				newSphere%mesh%particles%z0(j) = x0(3)

				newSphere%absVort%scalar(j) = relVortFn( x0(1), x0(2), x0(3) ) + &
					2.0_kreal * newSphere%rotationRate * x0(3) / newSphere%radius

				newSphere%relVort%scalar(j) = newSphere%absVort%scalar(j) - 2.0_kreal * &
					newSphere%rotationRate * newSphere%mesh%particles%z(j) / newSphere%radius

				if ( allocated(newSphere%tracers) ) then
					if ( nTracers == 1 ) then
						newSphere%tracers(1)%scalar(j) = tracerFn1( x0(1), x0(2), x0(3) )
					elseif ( nTracers == 2 ) then
						newSphere%tracers(1)%scalar(j) = tracerFn1( x0(1), x0(2), x0(3) )
						newSphere%tracers(2)%scalar(j) = tracerFn2( x0(1), x0(2), x0(3) )
					endif
					do k = nTracers+1, size(newSphere%tracers)
						if ( newSphere%tracers(k)%nDIM == 1 ) then
							newSphere%tracers(k)%scalar(j) = InterpolateScalar(lon, lat, self%tracerSource(k), &
																oldSphere%mesh, self%delTri, oldSphere%tracers(k) )
						else
							vecT = InterpolateVector( lon, lat, self%tracerSource(k), oldSphere%mesh, self%delTri, &
													   oldSphere%tracers(k) )
							newSphere%tracers(k)%xComp(j) = vecT(1)
							newSphere%tracers(k)%yComp(j) = vecT(2)
							newSphere%tracers(k)%zComp(j) = vecT(3)
						endif
					enddo
				endif
			enddo
		enddo
		newSphere%relVort%N = newSphere%mesh%particles%N
		newSphere%absVort%N = newSphere%mesh%particles%N
		if ( allocated(newSphere%tracers) ) then
			do k = 1, size(newSphere%tracers)
				newSphere%tracers(k)%N = newSphere%mesh%particles%N
			enddo
		endif

		call LoadBalance(newSphere%mpiParticles, newSphere%mesh%particles%N, numProcs)
		call Delete(refine)
	endif

	!
	!	set velocity and stream functions
	!
	call SetVelocityonMesh( newSphere )
	call SetStreamFunctionsOnMesh( newSphere )

	write(logString,'(A,I0.2,A)') 'debugRemeshOutput_',remeshCounter, '.vtk'
	call OutputToVTK(newSphere, logString)

	call LogMessage( log, DEBUG_LOGGING_LEVEL, trim(logkey)//" ", " Lagrangian Remesh Complete.")
end subroutine

! All AMR done on tracer 1
subroutine DirectRemeshTransport(self, oldSphere, newSphere, AMR, velFn, t, divFn, field1, flagFn1, tol1, desc1, &
	field2, flagFn2, tol2, desc2 )
	type(TransportRemesh), intent(in) :: self
	type(TransportMesh), intent(in) :: oldSphere
	type(TransportMesh), intent(inout) :: newSphere
	logical(klog), intent(in) :: AMR
	procedure(vectorFnOf3DSpaceAndTime) :: velFn
	real(kreal), intent(in) :: t
	procedure(scalarFnOf3DSpaceAndTime), optional :: divFn
	type(Field), intent(inout), optional :: field1
	procedure(FlagFunction), optional :: flagFn1
	real(kreal), intent(in), optional :: tol1
	character(len=*), intent(in), optional :: desc1
	type(Field), intent(inout), optional :: field2
	procedure(FlagFunction), optional :: flagFn2
	real(kreal), intent(in), optional :: tol2
	character(len=*), intent(in), optional :: desc2
	!
	integer(kint) :: i, j, k
	real(kreal) :: lon, lat
	real(kreal), dimension(3) :: x0, vec
	type(RefineSetup) :: refine
	integer(kint) :: nParticlesBefore, nParticlesAfter
	integer(kint) :: amrVarCount, mpiErrCode

	remeshCounter = remeshCounter + 1

	!
	!	interpolate to new uniform mesh
	!
	newSphere%density%N = newSphere%mesh%particles%N
	do i = 1, size(newSphere%tracers)
		newSphere%tracers(i)%N = newSphere%mesh%particles%N
	enddo
	do i = newSphere%mpiParticles%indexStart(procRank), newSphere%mpiParticles%indexEnd(procRank)
		lon = Longitude( newSphere%mesh%particles%x(i), newSphere%mesh%particles%y(i), newSphere%mesh%particles%z(i))
		lat = Latitude( newSphere%mesh%particles%x(i), newSphere%mesh%particles%y(i), newSphere%mesh%particles%z(i))

		x0 = InterpolateLagParam( lon, lat, self%lagParamSource, oldSphere%mesh, self%delTri)
		newSphere%mesh%particles%x0(i) = x0(1)
		newSphere%mesh%particles%y0(i) = x0(2)
		newSphere%mesh%particles%z0(i) = x0(3)

		newSphere%density%scalar(i) = InterpolateScalar( lon, lat, self%densitySource, oldSphere%mesh, &
											self%delTri, oldSphere%density)

		do j = 1, size(newSphere%tracers)
			if ( newSphere%tracers(j)%nDim == 1 ) then
				newSphere%tracers(j)%scalar(i) = InterpolateScalar(lon, lat, self%tracerSource(j), oldSphere%mesh, &
					self%delTri, oldSphere%tracers(j) )
			else
				vec = InterpolateVector(lon, lat, self%tracerSource(j), oldSphere%mesh, self%delTri,&
					 oldSphere%tracers(j))
				newSphere%tracers(j)%xComp(i) = vec(1)
				newSphere%tracers(j)%yComp(i) = vec(2)
				newSphere%tracers(j)%zComp(i) = vec(3)
			endif
		enddo
	enddo

	do i = 0, numProcs - 1
		! broadcast density
		call MPI_BCAST(newSphere%density%scalar(newSphere%mpiParticles%indexStart(i):newSphere%mpiParticles%indexEnd(i)), &
			newSphere%mpiParticles%messageLength(i), MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)
		! broadcast tracers
		do j = 1, size(newSphere%tracers)
			if ( newSphere%tracers(j)%nDim == 1) then
				call MPI_BCAST(newSphere%tracers(j)%scalar(newSphere%mpiParticles%indexStart(i):newSphere%mpiParticles%indexEnd(i)), &
					newSphere%mpiParticles%messageLength(i), MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)
			else
				call MPI_BCAST(newSphere%tracers(j)%xComp(newSphere%mpiParticles%indexStart(i):newSphere%mpiParticles%indexEnd(i)), &
					newSphere%mpiParticles%messageLength(i), MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)
				call MPI_BCAST(newSphere%tracers(j)%yComp(newSphere%mpiParticles%indexStart(i):newSphere%mpiParticles%indexEnd(i)), &
					newSphere%mpiParticles%messageLength(i), MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)
				call MPI_BCAST(newSphere%tracers(j)%zComp(newSphere%mpiParticles%indexStart(i):newSphere%mpiParticles%indexEnd(i)), &
					newSphere%mpiParticles%messageLength(i), MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)
			endif
		enddo
	enddo

	!
	!	AMR
	!
	if ( AMR ) then
		amrVarCount = 0
		if ( (present(field1) .AND. present(flagFn1)) .AND. (present(tol1 ) .AND. present(desc1) ) ) then
			amrVarCount = 1
			if ( (present(field2) .AND. present(flagFn2)) .AND. (present(tol2) .AND. present(desc2) ) ) then
				amrVarCount = 2
			endif
		endif

		if ( amrVarCount > 0 ) then
			call New(refine, newSphere%mesh%faces%N_Max)
			do i = 1, newSphere%mesh%amrLimit
				nParticlesBefore = newSphere%mesh%particles%N

				if ( amrVarCount == 1 ) then
					call IterateMeshRefinementOneVariable( refine, newSphere%mesh, field1, flagFn1, tol1, desc1, &
						nParticlesBefore, nParticlesAfter )
				elseif ( amrVarCount == 2) then
					call IterateMeshRefinementTwoVariables( refine, newSphere%mesh, field1, flagFn1, tol1, desc1, &
						field2, flagFn2, tol2, desc2, nParticlesBefore, nParticlesAfter )
				endif
			enddo

			do j = nParticlesBefore + 1, nParticlesAfter
				lon = Longitude( newSphere%mesh%particles%x(j), newSphere%mesh%particles%y(j), &
									 newSphere%mesh%particles%z(j))
				lat = Latitude( newSphere%mesh%particles%x(j), newSphere%mesh%particles%y(j), &
									newSphere%mesh%particles%z(j))
				x0 = InterpolateLagParam( lon, lat, self%lagParamSource, oldSphere%mesh, self%delTri)
				newSphere%mesh%particles%x0(j) = x0(1)
				newSphere%mesh%particles%y0(j) = x0(2)
				newSphere%mesh%particles%z0(j) = x0(3)

				newSphere%density%scalar(j) = InterpolateScalar( lon, lat, self%densitySource, oldSphere%mesh, &
													self%delTri, oldSphere%density)

				do k = 1, size(newSphere%tracers)
					if ( newSphere%tracers(k)%nDim == 1 ) then
						newSphere%tracers(k)%scalar(j) = InterpolateScalar( lon, lat, self%tracerSource(k), &
								 oldSphere%mesh, self%delTri, oldSphere%tracers(k) )
					else
						vec = InterpolateVector(lon, lat, self%tracerSource(k), oldSphere%mesh, self%delTri, &
							oldSphere%tracers(k) )
						newSphere%tracers(k)%xComp(j) = vec(1)
						newSphere%tracers(k)%yComp(j) = vec(2)
						newSphere%tracers(k)%zComp(j) = vec(3)
					endif
				enddo
			enddo

			newSphere%density%N = newSphere%mesh%particles%N
			do k = 1, size(newSphere%tracers)
				newSphere%tracers(k)%N = newSphere%mesh%particles%N
			enddo
			call LoadBalance( newSphere%mpiParticles, newSphere%mesh%particles%N, numProcs)

			call Delete(refine)
		endif
	endif!AMR

	call SetVelocityOnMesh( newSphere, velFn, t )
	if ( present(divFn) ) then
		call SetDivergenceOnMesh( newSphere, divFn, t )
	endif
end subroutine

subroutine LagrangianRemeshTransportWithFunctions( self, oldSphere, newSphere, AMR, velFn, t, divFn, &
	tracerFn1, flagFn1, tol1, desc1, tracerFn2, flagFn2, tol2, desc2, RefineFlowMapYN, flowMapVarTol )
	type(TransportRemesh), intent(inout) :: self
	type(TransportMesh), intent(inout) :: oldSphere
	type(TransportMesh), intent(inout) :: newSphere
	logical(klog), intent(in) :: AMR
	procedure(vectorFnOf3DSpaceAndTime) :: velFn
	real(kreal), intent(in) :: t
	procedure(scalarFnOf3DSpaceAndTime), optional :: divFn
	procedure(scalarFnOf3DSpace), optional :: tracerFn1
	procedure(FlagFunction), optional :: flagFn1
	real(kreal), intent(in), optional :: tol1
	character(len=*), intent(in), optional :: desc1
	procedure(scalarFnOf3DSpace), optional :: tracerFn2
	procedure(FlagFunction), optional :: flagFn2
	real(kreal), intent(in), optional :: tol2
	character(len=*), intent(in), optional :: desc2
	logical(klog), intent(in), optional :: RefineFlowMapYN
	real(kreal), intent(in), optional :: flowMapVarTol
	!
	integer(kint) :: nLagTracers
	logical(klog) :: doFlowMapRefinement
	logical(klog) :: twoRefinements
	integer(kint) :: i, j, k
	real(kreal) :: lon, lat
	real(kreal), dimension(3) :: x0
	type(RefineSetup) :: refine
	integer(kint) :: nParticlesBefore, nParticlesAfter
	real(kreal), dimension(3) :: vec
	type(DelaunayTriangulation) :: lagDelTri
	logical(klog) :: useLagCoords
	integer(kint) :: mpiErrCode

	remeshCounter = remeshCounter + 1

	if ( present(tracerFn1) .and. newSphere%tracers(1)%nDim == 1 ) then
		nLagTracers = 1
		if ( present(tracerFn2) .and. newSphere%tracers(2)%nDim == 1 ) then
			nLagTracers = 2
		endif
	else
		nLagTracers = 0
	endif

	newSphere%density%N = newSphere%mesh%particles%N
	if ( allocated(newSphere%tracers) ) then
		do i = 1, size(newSphere%tracers)
			newSphere%tracers(i)%N = newSphere%mesh%particles%N
		enddo
	endif


	do i = newSphere%mpiParticles%indexStart(procRank), newSphere%mpiParticles%indexEnd(procRank)
		lon = Longitude(newSphere%mesh%particles%x(i), newSphere%mesh%particles%y(i), newSphere%mesh%particles%z(i))
		lat = Latitude( newSphere%mesh%particles%x(i), newSphere%mesh%particles%y(i), newSphere%mesh%particles%z(i))
		x0 = InterpolateLagParam( lon, lat, self%lagParamSource, oldSphere%mesh, self%delTri)
		newSphere%mesh%particles%x0(i) = x0(1)
		newSphere%mesh%particles%y0(i) = x0(2)
		newSphere%mesh%particles%z0(i) = x0(3)

		!
		!	direct interpolation for density
		!
		newSphere%density%scalar(i) = InterpolateScalar( lon, lat, self%densitySource, &
			 oldSphere%mesh, self%delTri, oldSphere%density)

		if ( allocated(newSphere%tracers) ) then
			!
			!	indirect interpolation for tracers 1 and 2
			!
			if ( nLagTracers == 1 ) then
				newSphere%tracers(1)%scalar(i) = tracerFn1( x0(1), x0(2), x0(3) )
			elseif ( nLagTracers == 2 ) then
				newSphere%tracers(1)%scalar(i) = tracerFn1( x0(1), x0(2), x0(3) )
				newSphere%tracers(2)%scalar(i) = tracerFn2( x0(1), x0(2), x0(3) )
			endif
			do j = nLagTracers + 1, size(newSphere%tracers)
				if ( newSphere%tracers(j)%nDim == 1 ) then
					newSphere%tracers(j)%scalar(i) = InterpolateScalar(lon, lat, self%tracerSource(j), oldSphere%mesh, &
																	   self%delTri, oldSphere%tracers(j) )
				else
					vec = InterpolateVector(lon, lat, self%tracerSource(j), oldSphere%mesh, &
											 self%delTri, oldSphere%tracers(j) )
					newSphere%tracers(j)%xComp(i) = vec(1)
					newSphere%tracers(j)%yComp(i) = vec(2)
					newSphere%tracers(j)%zComp(i) = vec(3)
				endif
			enddo
		endif
	enddo

	do i = 0, numProcs - 1
		! broadcast density
		call MPI_BCAST(newSphere%density%scalar(newSphere%mpiParticles%indexStart(i):newSphere%mpiParticles%indexEnd(i)), &
			newSphere%mpiParticles%messageLength(i), MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)
		! broadcast tracers
		do j = 1, size(newSphere%tracers)
			if ( newSphere%tracers(j)%nDim == 1) then
				call MPI_BCAST(newSphere%tracers(j)%scalar(newSphere%mpiParticles%indexStart(i):newSphere%mpiParticles%indexEnd(i)), &
					newSphere%mpiParticles%messageLength(i), MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)
			else
				call MPI_BCAST(newSphere%tracers(j)%xComp(newSphere%mpiParticles%indexStart(i):newSphere%mpiParticles%indexEnd(i)), &
					newSphere%mpiParticles%messageLength(i), MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)
				call MPI_BCAST(newSphere%tracers(j)%yComp(newSphere%mpiParticles%indexStart(i):newSphere%mpiParticles%indexEnd(i)), &
					newSphere%mpiParticles%messageLength(i), MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)
				call MPI_BCAST(newSphere%tracers(j)%zComp(newSphere%mpiParticles%indexStart(i):newSphere%mpiParticles%indexEnd(i)), &
					newSphere%mpiParticles%messageLength(i), MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)
			endif
		enddo
	enddo

	!
	!	AMR
	!
	if ( AMR ) then
		doFlowMapRefinement = ( RefineFlowMapYN .AND. present(flowMapVarTol) )
		twoRefinements = ( present(flagFn1) .AND. present(flagFn2) )

		call New( refine, newSphere%mesh%faces%N_Max)

		do i = 1, newSphere%mesh%amrLimit
			nParticlesBefore = newSphere%mesh%particles%N

			if ( twoRefinements ) then
				if ( doFlowMapRefinement ) then
					call IterateMeshRefinementTwoVariablesAndFlowMap( refine, newSphere%mesh, newSphere%tracers(1), &
						flagFn1, tol1, desc1, newSphere%tracers(1), flagFn2, tol2, desc2, &
						flowMapVarTol, nParticlesBefore, nParticlesAfter )
				else
					call IterateMeshRefinementTwoVariables( refine, newSphere%mesh, newSphere%tracers(1), &
						flagFn1, tol1, desc1, newSphere%tracers(1), flagFn2, tol2, desc2, &
						nParticlesBefore, nParticlesAfter )
				endif
			else
				if ( doFlowMapRefinement ) then
					call IterateMeshRefinementOneVariableAndFlowMap( refine, newSphere%mesh, newSphere%tracers(1), &
						flagFn1, tol1, desc1, flowMapVarTol, nParticlesBefore, nParticlesAfter )
				else
					call IterateMeshRefinementOneVariable( refine, newSphere%mesh, newSphere%tracers(1), &
						flagFn1, tol1, desc1, nParticlesBefore, nParticlesAfter )
				endif
			endif
		enddo

		do j = nParticlesBefore + 1, nParticlesAfter
			lon = Longitude( newSphere%mesh%particles%x(j), newSphere%mesh%particles%y(j), &
								 newSphere%mesh%particles%z(j))
			lat = Latitude( newSphere%mesh%particles%x(j), newSphere%mesh%particles%y(j), &
								newSphere%mesh%particles%z(j))
			x0 = InterpolateLagParam(lon, lat, self%lagParamSource, oldSphere%mesh, self%delTri )
			newSphere%mesh%particles%x0(j) = x0(1)
			newSphere%mesh%particles%y0(j) = x0(2)
			newSphere%mesh%particles%z0(j) = x0(3)

			newSphere%density%scalar(j) = InterpolateScalar( lon, lat, self%densitySource, oldSphere%mesh, &
												self%delTri, oldSphere%density)

			if ( allocated(newSphere%tracers) ) then
				if ( nLagTracers == 1 ) then
					newSphere%tracers(1)%scalar(j) = tracerFn1( x0(1), x0(2), x0(3) )
				elseif ( nLagTracers == 2 ) then
					newSphere%tracers(1)%scalar(j) = tracerFn1( x0(1), x0(2), x0(3) )
					newSphere%tracers(2)%scalar(j) = tracerFn2( x0(1), x0(2), x0(3) )
				endif
				do k = nLagTracers + 1, size(newSphere%tracers)
					if ( newSphere%tracers(k)%nDim == 1 ) then
						newSphere%tracers(k)%scalar(j) = InterpolateScalar( lon, lat, self%tracerSource(k), &
															oldSphere%mesh, self%delTri, oldSphere%tracers(k) )
					else
						vec = InterpolateVector( lon, lat, self%tracerSource(k), oldSphere%mesh, self%delTri, &
												oldSphere%tracers(k) )
						newSphere%tracers(k)%xComp(j) = vec(1)
						newSphere%tracers(k)%yComp(j) = vec(2)
						newSphere%tracers(k)%zComp(j) = vec(3)
					endif
				enddo
			endif
		enddo

		newSphere%density%N = newSphere%mesh%particles%N
		if ( allocated(newSphere%tracers) ) then
			do k = 1, size(newSphere%tracers)
				newSphere%tracers(k)%N = newSphere%mesh%particles%N
			enddo
		endif

		call LoadBalance(newSphere%mpiParticles, newSphere%mesh%particles%N, numProcs )
		call Delete(refine)
	endif!AMR

	call SetVelocityOnMesh( newSphere, velFn, t )
	if ( present(divFn) ) then
		call SetDivergenceOnMesh( newSphere, divFn, t )
	else
		call SetDivergenceOnMesh(newSphere)
	endif

end subroutine

!> @brief Performs a remesh/remap of an LPM @ref SphereBVE simulation using indirect interpolation for vorticity variables and up to two scalar tracer variables in a BVE mesh.  Additional tracers are directly interpolated.
!> Remaps to reference time t = t_{rm}, where t_{rm} is time of definition for a reference mesh (numerical solution).
!> Note that the reference mesh's Lagrangian parameter has been reset, so that x = x0, y = y0, z = z0 at t = t_{rm}
!>
!> @param[inout] newSphere target @ref SphereBVE mesh (note that this must have been allocated prior to calling this subroutine)
!> @param[in] oldSphere source @ref SphereBVE mesh
!> @param[in] refSphere @ref SphereBVE reference mesh
!> @param[in] refRemesh remesh data structure associated with refSphere
!> @param[in] AMR .TRUE. if adaptive refinement will be used
!> @param[in] flagFn1 FlagFunction for vorticity refinement, must have same interface as refinementmodule::FlagFunction
!> @param[in] tol1 tolerance for flagFn1
!> @param[in] desc1 description of AMR criterion used for flagFn1
!> @param[in] flagFn2 FlagFunction for vorticity refinement, must have same interface as refinementmodule::FlagFunction
!> @param[in] tol2 tolerance for flagFn2
!> @param[in] desc2 description of AMR criterion used for flagFn2
!> @param[in] RefineFLowMapYN True if refinement of the flow map will be used
!> @param[in] flowMapVarTol tolerance value for Lagrangian coordinate variation per face
subroutine LagrangianRemeshBVEToReferenceMesh( newSphere, oldSphere, refSphere, refRemesh, AMR, flagFn1, tol1, desc1, &
	flagFn2, tol2, desc2, RefineFlowMapYN, flowMapVarTol )
	type(BVEMesh), intent(inout) :: newSphere
	type(BVEMesh), intent(inout) :: oldSphere
	type(BVEMesh), intent(in) :: refSphere
	type(BVERemeshSource), intent(in) :: refRemesh
	logical(klog), intent(in) :: AMR
	procedure(FlagFunction) :: flagFn1
	real(kreal), intent(in) :: tol1
	character(len=*), intent(in) :: desc1
	procedure(FlagFunction), optional :: flagFn2
	real(kreal), intent(in), optional :: tol2
	character(len=*), intent(in), optional :: desc2
	logical(klog), intent(in), optional :: RefineFlowMapYN
	real(kreal), intent(in), optional :: flowMapVarTol
	!
	integer(kint) :: i, j, k, nNew
	logical(klog) :: refineVorticityTwice, doFlowMapRefinement
	type(RefineSetup) :: refine
	type(BVERemeshSource) :: remesh
	integer(kint) :: nParticlesBefore, nParticlesAfter
	real(kreal), dimension(3) :: vecT, x0
	real(kreal) :: lon, lat, lon0, lat0

	nNew = newSphere%mesh%particles%N
	newSphere%relVort%N = nNew
	newSphere%absVort%N = nNew
	if ( allocated(newSphere%tracers)) then
		do i = 1, size(newSphere%tracers)
			newSphere%tracers(i)%N = nNew
		enddo
	endif

	call LogMessage(log, DEBUG_LOGGING_LEVEL, &
		trim(logKey)//" LagRemeshToRef : remeshing to base uniform mesh, nParticles = ", nNew)
	call LogMessage(log, DEBUG_LOGGING_LEVEL, &
		trim(logKey)//" LagRemeshToRef : found new nTracers = ", size(newSphere%tracers))
	call LogMessage(log, DEBUG_LOGGING_LEVEL, &
		trim(logKey)//" LagRemeshToRef : found old nTracers = ", size(oldSphere%tracers))

	call newBVE(remesh, oldSphere)
	call LogMessage(log, DEBUG_LOGGING_LEVEL, trim(logKey)//" LagRemeshToRef : ", "remesh object ready.")

	do i = 1, nNew
		!
		!	interpolate Lagrangian parameter from old mesh to new mesh
		!		NOTE: This Lagrangian parameter was reset (x = x0) at the time when reference sphere was defined.
		!
		lon = Longitude( newSphere%mesh%particles%x(i), newSphere%mesh%particles%y(i), newSphere%mesh%particles%z(i))
		lat = Latitude( newSphere%mesh%particles%x(i), newSphere%mesh%particles%y(i), newSphere%mesh%particles%z(i))
		x0 = InterpolateLagParam( lon, lat, remesh%lagParamSource, oldSphere%mesh, remesh%delTri)

		newSphere%mesh%particles%x0(i) = x0(1)
		newSphere%mesh%particles%y0(i) = x0(2)
		newSphere%mesh%particles%z0(i) = x0(3)

		lon0 = Longitude(x0)
		lat0 = Latitude(x0)

		!
		!	interpolate vorticity from reference mesh
		!
		newSphere%absVort%scalar(i) = InterpolateScalar( lon0, lat0, refRemesh%absVortSource, refSphere%mesh, &
			refRemesh%delTri, refSphere%absVort)
		newSphere%relVort%scalar(i) = newSphere%absVort%scalar(i) - 2.0_kreal * newSphere%rotationRate * &
			newSphere%mesh%particles%z(i) / newSphere%radius

		!
		!	interpolate tracers from reference mesh
		!
		if ( allocated(newSphere%tracers) ) then
			do j = 1, size(newSphere%tracers)
				if ( newSphere%tracers(j)%nDim == 1 ) then
					newSphere%tracers(j)%scalar(i) = InterpolateScalar(lon0, lat0, refRemesh%tracerSource(j), &
						refSphere%mesh, refRemesh%delTri, refSphere%tracers(j) )
				else
					vecT = InterpolateVector(lon0, lat0, refRemesh%tracerSource(j), refSphere%mesh, refRemesh%delTri, &
						refSphere%tracers(j) )
					newSphere%tracers(j)%xComp(i) = vecT(1)
					newSphere%tracers(j)%yComp(i) = vecT(2)
					newSphere%tracers(j)%zComp(i) = vecT(3)
				endif
			enddo
		endif
	enddo

	if ( AMR ) then
		refineVorticityTwice = ( present(flagFn2) .AND. ( present(tol2) .AND. present(desc2)))
		doFlowMapRefinement = ( RefineFlowMapYN .AND. present(flowMapVarTol) )

		call New(refine, newSphere%mesh%faces%N_Max)

		do i = 1, newSphere%mesh%amrLimit
			nParticlesBefore = newSphere%mesh%particles%N

			if ( refineVorticityTwice ) then
				if ( doFlowMapRefinement ) then
					call IterateMeshRefinementTwoVariablesAndFlowMap( refine, newSphere%mesh, newSphere%relVort, &
						flagFn1, tol1, desc1, newSphere%relVort, flagFn2, tol2, desc2, &
						flowMapVarTol, nParticlesBefore, nParticlesAfter )
				else
					call IterateMeshRefinementTwoVariables(refine, newSphere%mesh, newSphere%relVort, &
						flagFn1, tol1, desc1, newSphere%relVort, flagFn2, tol2, desc2, nParticlesBefore, nParticlesAfter)
				endif
			else
				if ( doFlowMapRefinement ) then
					call IterateMeshRefinementOneVariableAndFlowMap( refine, newSphere%mesh, newSphere%relVort, &
						flagFn1, tol1, desc1, flowMapVarTol, nParticlesBefore, nParticlesAfter)
				else
					call IterateMeshRefinementOneVariable( refine, newSphere%mesh, newSphere%relVort, flagFn1, &
						tol1, desc1, nParticlesBefore, nParticlesAfter)
				endif
			endif

			do j = nParticlesBefore + 1, nParticlesAfter
				lon = Longitude( newSphere%mesh%particles%x(j), newSphere%mesh%particles%y(j), newSphere%mesh%particles%z(j))
				lat = Latitude( newSphere%mesh%particles%x(j), newSphere%mesh%particles%y(j), newSPhere%mesh%particles%z(j))
				x0 = InterpolateLagParam( lon, lat, remesh%lagParamSource, oldSphere%mesh, remesh%delTri)
				newSphere%mesh%particles%x0(j) = x0(1)
				newSphere%mesh%particles%y0(j) = x0(2)
				newSphere%mesh%particles%z0(j) = x0(3)

				lon0 = Longitude(x0)
				lat0 = Latitude(x0)

				newSphere%absVort%scalar(j) = InterpolateScalar( lon0, lat0, refRemesh%absVortSource, refSphere%mesh, &
					refRemesh%delTri, refSphere%absVort )
				newSphere%relVort%scalar(j) = newSphere%absVort%scalar(j) - 2.0_kreal * newSphere%rotationRate * &
					newSphere%mesh%particles%z(j) / newSphere%radius

				if ( allocated(newSphere%tracers) ) then
					do k = 1, size(newSphere%tracers)
						if ( newSphere%tracers(k)%nDim == 1 ) then
							newSphere%tracers(k)%scalar(j) = InterpolateScalar( lon0, lat0, refRemesh%tracerSource(k), &
								refSphere%mesh, refRemesh%delTri, refSphere%tracers(k))
						else
							vecT = InterpolateVector( lon0, lat0, refRemesh%tracerSource(k), refSphere%mesh, &
								refRemesh%delTri, refSphere%tracers(k) )
							newSphere%tracers(k)%xComp(j) = vecT(1)
							newSphere%tracers(k)%yComp(j) = vecT(2)
							newSphere%tracers(k)%zComp(j) = vecT(3)
						endif
					enddo
				endif
			enddo
		enddo

		nNew = newSphere%mesh%particles%N
		newSphere%absVort%N = nNew
		newSphere%relVort%N = nNew
		if ( allocated(newSphere%tracers) ) then
			do i = 1, size(newSphere%tracers)
				newSphere%tracers(i)%N = nNew
			enddo
		endif

		call LoadBalance( newSphere%mpiParticles, nNew, numProcs)
		call Delete(refine)
	endif

	!
	!	set velocity and stream functions
	!
	call StartSection(log, "NEW SPHERE STATUS")
	call LogStats(newSphere,log)
	call EndSection(log)

	call LogMessage(log, DEBUG_LOGGING_LEVEL, &
		trim(logKey)//" LagRemeshToRef : ", " setting velocity on new mesh.")

	call SetVelocityonMesh( newSphere )
	call LogMessage(log, DEBUG_LOGGING_LEVEL, &
		trim(logKey)//" LagRemeshToRef : ", " velocity done.")
		call LogMessage(log, DEBUG_LOGGING_LEVEL, &
		trim(logKey)//" LagRemeshToRef : ", " setting stream functions on new mesh.")
	call SetStreamFunctionsOnMesh( newSphere )
		call LogMessage(log, DEBUG_LOGGING_LEVEL, &
		trim(logKey)//" LagRemeshToRef : ", " stream functions done.")


	call LogMessage(log, DEBUG_LOGGING_LEVEL, &
		trim(logKey)//" LagRemeshToRef : ", " deleting local remesh tool.")
	call Delete(remesh)
end subroutine

subroutine FTLECalc( amesh, faceIndex,FTLE_,FTLE_Error_ )
    type(PolyMesh2d), intent(in) :: amesh
		integer(kint), intent(in) :: faceIndex
		real(kreal),intent(out) :: FTLE_,FTLE_Error_
    !
    integer(kint) :: pIndex, i
		real(kreal),dimension(1:3) :: x0i,xi,x0i_t,xi_t,x0face,xface,x0face_t,xface_t
		real(kreal), dimension(1:2,1:2) :: FlowMapGrad
		real(kreal), dimension(1:4) :: x0,xf,y0,yf,x0filter
		real(kreal) :: DyDy0,DxDy0,DxDx0,DyDx0
		real(kreal) :: CG_11,CG_22,CG_12,Eigmin,Eigmax,x0_center,y0_center,xf_center,yf_center
		integer(kint) :: thmin1,thmin2,thmax1,thmax2,loc1,loc2,loc3,loc4,ly

		real(kreal),dimension(1:3)::Zvec,Zvec0
  	real(kreal),dimension(1:3,1:3)::Eye,S,R,Eye0,S0,R0
    real(kreal) :: pi
    real(kreal) :: Dx,Dy,Dx0,Dy0,Dxcross,Dycross,a,b

		real(kreal), dimension(1:4) :: x0_t,x_t,y0_t,y_t
		real(kreal),dimension(1:2) :: Zvec_2,xaxis_2,yaxis_2
		real(kreal),dimension(1:2,1:2)::Eye_2,S_2,R_2

		pi=atan(1.d0)*4.d0
		pIndex = amesh%faces%centerParticle(faceIndex) ! Indices of the active particles
		x0face=	LagCoord (amesh%particles, pIndex)
		xface=	PhysCoord(amesh%particles, pIndex)
		! Renormalize advected locations to remove any error
		xface=xface/dsqrt(xface(1)**2+xface(2)**2+xface(3)**2)
		! Rotate the coordinate system such that the normal aligns with
		! [0 0 1].
		! The tangent space will then be x-y plane coordinates.
		Zvec=0.d0;Zvec(3)=1.d0; Zvec=Zvec/dsqrt(Zvec(1)**2+Zvec(2)**2+Zvec(3)**2)
		Eye=0.d0; Eye(1,1)=1.d0;Eye(2,2)=1.d0;Eye(3,3)=1.d0

		! The tensor R0 is the rotation matrix in reference/ t=0 configuration
		call reflection (S,Eye,Zvec+x0face);call reflection (R0,S,Zvec)
		! The tensor R is the rotation matrix in current configuration
		call reflection (S,Eye,Zvec+xface); call reflection (R,S,Zvec)

		!#$! x0face_t=x0face;xface_t=xface

		! Rotate the face centers. Don't really need these so commented with !#$!
		! a) reference configuration
		!#$! x0face(1)=R0(1,1)*x0face_t(1)+R0(1,2)*x0face_t(2)+R0(1,3)*x0face_t(3)
		!#$! x0face(2)=R0(2,1)*x0face_t(1)+R0(2,2)*x0face_t(2)+R0(2,3)*x0face_t(3)
		!#$! x0face(3)=R0(3,1)*x0face_t(1)+R0(3,2)*x0face_t(2)+R0(3,3)*x0face_t(3)
		! b) current configuration
		!#$! xface (1)=R (1,1)*xface_t (1)+R (1,2)*xface_t (2)+R (1,3)*xface_t (3)
		!#$! xface (2)=R (2,1)*xface_t (1)+R (2,2)*xface_t (2)+R (2,3)*xface_t (3)
		!#$! xface (3)=R (3,1)*xface_t (1)+R (3,2)*xface_t (2)+R (3,3)*xface_t (3)

		!#$!x0_center=x0face(1);y0_center=x0face(3);
		!#$!xf_center=xface(1);	yf_center=xface(3);

		do i = 1, amesh%faceKind ! Coordinates of the passive particles/ vertices
        x0i = LagCoord(amesh%particles, amesh%faces%vertices(i,faceIndex))
				xi  = PhysCoord(amesh%particles, amesh%faces%vertices(i,faceIndex))

				! Normalize to remove errors
				x0i=x0i/dsqrt(x0i(1)**2+x0i(2)**2+x0i(3)**2)
				xi=xi/dsqrt(xi(1)**2+xi(2)**2+xi(3)**2)

				x0i_t=x0i;xi_t=xi
				! Rotate the face vertices
				! a) reference configuration
			  x0i(1)=R0(1,1)*x0i_t(1)+R0(1,2)*x0i_t(2)+R0(1,3)*x0i_t(3)
				x0i(2)=R0(2,1)*x0i_t(1)+R0(2,2)*x0i_t(2)+R0(2,3)*x0i_t(3)
				x0i(3)=R0(3,1)*x0i_t(1)+R0(3,2)*x0i_t(2)+R0(3,3)*x0i_t(3)
				! a) current configuration
				xi (1)=R (1,1)*xi_t (1)+R (1,2)*xi_t (2)+R (1,3)*xi_t (3)
				xi (2)=R (2,1)*xi_t (1)+R (2,2)*xi_t (2)+R (2,3)*xi_t (3)
				xi (3)=R (3,1)*xi_t (1)+R (3,2)*xi_t (2)+R (3,3)*xi_t (3)

				x0(i)=x0i(1);			y0(i)=x0i(2)
				xf(i)=xi(1);			yf (i)=xi(2)
	  enddo
		! ========================================================================
		! Locating points on the quadrilateral inscribed on the sphere
		! 1, 2, 3 and 4 in clockwise direction
		! thmin1=minloc(x0,dim=1);
		! x0filter=10.d0 ! Largest value of xf will be pi
		! do i=1,4
		! if (i.ne.thmin1) x0filter(i)=x0(i)
		! end do
		! thmin2=minloc(x0filter,dim=1);
		! do i=1,4
		! 	if( i.ne.thmin1.and.i.ne.thmin2) thmax1=i
		! end do
		! do i=1,4
		! 	if( i.ne.thmin1.and.i.ne.thmin2.and.i.ne.thmax1) thmax2=i
		! end do
		! loc1=thmin1;loc2=thmin2;
		! if (y0(thmin1)>y0(thmin2)) then
		! 	loc1=thmin2;loc2=thmin1;
		! endif
		! loc3=thmax1;loc4=thmax2;
		! if (y0(thmax2)>y0(thmax1)) then
		! 	loc3=thmax2;loc4=thmax1;
		! endif
		loc1=1;loc2=2;loc3=3;loc4=4;

		! Linear map to a properly oriented rectangle
		! OriginShift
		x0=x0-x0(loc1);y0=y0-y0(loc1);!x0_center=x0_center-x0(loc1);y0_center=y0_center-y0(loc1);
		xf=xf-xf(loc1);yf=yf-yf(loc1);!xf_center=xf_center-xf(loc1);yf_center=yf_center-yf(loc1) ;

		Eye_2=0.d0; Eye_2(1,1)=1.d0;Eye_2(2,2)=1.d0;S_2=Eye_2;R_2=Eye_2;
		Zvec_2=0.d0;Zvec_2(1)=x0(loc4)-x0(loc1); Zvec_2(2)=y0(loc4)-y0(loc1);
		Zvec_2=Zvec_2/dsqrt(Zvec_2(1)**2+Zvec_2(2)**2)
		xaxis_2=0.d0;xaxis_2(1)=1.d0;xaxis_2=xaxis_2/dsqrt(xaxis_2(1)**2+xaxis_2(2)**2)
	  call reflection_2D (S_2,Eye_2,xaxis_2+Zvec_2)
		call reflection_2D (R_2,S_2,xaxis_2)
		x0_t(loc1)=R_2(1,1)*x0(loc1)+R_2(1,2)*y0(loc1)
		y0_t(loc1)=R_2(2,1)*x0(loc1)+R_2(2,2)*y0(loc1)
		x0_t(loc2)=R_2(1,1)*x0(loc2)+R_2(1,2)*y0(loc2)
		y0_t(loc2)=R_2(2,1)*x0(loc2)+R_2(2,2)*y0(loc2)
		x0_t(loc3)=R_2(1,1)*x0(loc3)+R_2(1,2)*y0(loc3)
		y0_t(loc3)=R_2(2,1)*x0(loc3)+R_2(2,2)*y0(loc3)
		x0_t(loc4)=R_2(1,1)*x0(loc4)+R_2(1,2)*y0(loc4)
		y0_t(loc4)=R_2(2,1)*x0(loc4)+R_2(2,2)*y0(loc4)

		x_t(loc1)=R_2(1,1)*xf(loc1)+R_2(1,2)*yf(loc1)
		y_t(loc1)=R_2(2,1)*xf(loc1)+R_2(2,2)*yf(loc1)
		x_t(loc2)=R_2(1,1)*xf(loc2)+R_2(1,2)*yf(loc2)
		y_t(loc2)=R_2(2,1)*xf(loc2)+R_2(2,2)*yf(loc2)
		x_t(loc3)=R_2(1,1)*xf(loc3)+R_2(1,2)*yf(loc3)
		y_t(loc3)=R_2(2,1)*xf(loc3)+R_2(2,2)*yf(loc3)
		x_t(loc4)=R_2(1,1)*xf(loc4)+R_2(1,2)*yf(loc4)
		y_t(loc4)=R_2(2,1)*xf(loc4)+R_2(2,2)*yf(loc4)

		x0=x0_t;y0=y0_t;xf=x_t;yf=y_t
	! ========================================================================
		Dx0=x0(loc4)-x0(loc1);Dx=xf(loc4)-xf(loc1);Dycross=yf(loc4)-yf(loc1);

		ly=loc2;if (abs(y0(loc3))>abs(y0(loc2))) ly=loc2;
		Dy0=y0(ly);
		a=(x0(ly)-x0(loc1))/(x0(loc4)-x0(loc1));
		b=1.d0-a
		Dy=yf(ly)-(yf(loc4)*a+yf(loc1)*b);
		Dxcross=xf(ly)-(xf(loc4)*a+xf(loc1)*b);

		! Dy0=0.5d0*(y0(loc2)+y0(loc3))
		! x0star=0.5d0*(x0(loc2)+x0(loc3))
		! a=(x0star-x0(loc1))/(x0(loc4)-x0(loc1));
		! b=1.d0-a
		! Dy=0.5d0*(yf(loc2)+yf(loc3))-(yf(loc1)*b+yf(loc4)*a);
		! Dxcross=0.5d0*(yf(loc2)+yf(loc3))-(xf(loc1)*b+xf(loc4)*a);

		! 	x014=x0_center;y014=y0(loc1)+(y0(loc4)-y0(loc1))/(x0(loc4)-x0(loc1))*(x014-x0(loc1));
		! 	a_14=dsqrt(((y014-y0(loc1))**2+(x014-x0(loc1))**2)/((y0(loc4)-y0(loc1))**2+(x0(loc4)-x0(loc1))**2))
		! 	x14=xf(1)+a_14*(xf(loc4)-xf(loc1));y14=yf(loc1)+a_14*(yf(loc4)-yf(loc1))
		!
		! 	x023=x0_center;y023=y0(loc2)+(y0(loc3)-y0(loc2))/(x0(loc3)-x0(loc2))*(x023-x0(loc2));
		! 	a_23=dsqrt(((y023-y0(loc2))**2+(x023-x0(loc2))**2)/((y0(loc3)-y0(loc2))**2+(x0(loc3)-x0(loc2))**2))
		! 	x23=xf(2)+a_23*(xf(loc3)-xf(loc2));y23=yf(loc2)+a_23*(yf(loc3)-yf(loc2))
		!
		! 	y012=y0_center;x012=x0(loc1)+(x0(loc2)-x0(loc1))/(y0(loc2)-y0(loc1))*(y012-y0(loc1));
		! 	a_12=dsqrt(((y012-y0(loc1))**2+(x012-x0(loc1))**2)/((y0(loc2)-y0(loc1))**2+(x0(loc2)-x0(loc1))**2))
		! 	x12=xf(loc1)+a_12*(xf(loc2)-xf(loc1));y12=yf(loc1)+a_12*(yf(loc2)-yf(loc1))
		!
		! 	y043=y0_center;x043=x0(loc4)+(x0(loc3)-x0(loc4))/(y0(loc3)-y0(loc4))*(y043-y0(loc4));
		! 	a_43=dsqrt(((y043-y0(loc4))**2+(x043-x0(loc4))**2)/((y0(loc3)-y0(loc4))**2+(x0(loc3)-x0(loc4))**2))
		! 	x43=xf(loc4)+a_43*(xf(loc3)-xf(loc4));y43=yf(loc4)+a_43*(yf(loc3)-yf(loc4))
		!
		! ! ========================================================================
		!
		! 	Dx0=x043-x012;Dx=x43-x12;Dycross=y43-y12;
		! 	Dy0=y023-y014;Dy=y23-y14;Dxcross=y23-y14;

		FlowMapGrad(1,1)=Dx/Dx0
		FlowMapGrad(2,1)=Dycross/Dx0
		FlowMapGrad(1,2)=Dxcross/Dy0
		FlowMapGrad(2,2)=Dy/Dy0

		! =========== Derivative of yf/ xf current w.r.t. reference ==========
    CG_11=FlowMapGrad(1,1)*FlowMapGrad(1,1)+FlowMapGrad(2,1)*FlowMapGrad(2,1)
    CG_22=FlowMapGrad(1,2)*FlowMapGrad(1,2)+FlowMapGrad(2,2)*FlowMapGrad(2,2)
    CG_12=FlowMapGrad(1,1)*FlowMapGrad(1,2)+FlowMapGrad(2,1)*FlowMapGrad(2,2)

		Eigmin=((CG_11+CG_22)-dsqrt(4.d0*CG_12**2+(CG_11-CG_22)**2))/2.d0
		Eigmax=((CG_11+CG_22)+dsqrt(4.d0*CG_12**2+(CG_11-CG_22)**2))/2.d0

    FTLE_=dlog(Eigmax)/amesh%t/2.d0! Fix it later /(2.d0*dtLastUpdate)
		FTLE_Error_=Eigmin*Eigmax-1.d0
		! if (faceIndex.eq.370.or.faceindex.eq.371.or.faceindex.eq.383)then
		print*,faceIndex,'FTLE_',Eigmin,Eigmax,Eigmin*Eigmax-1.d0,FTLE_
		!  endif
end subroutine FTLECalc

subroutine reflection (R,A,n)

	real(kreal), dimension(1:3,1:3),intent(in) :: A
	real(kreal), dimension(1:3,1:3),intent(out) :: R
	real(kreal), dimension(1:3),intent(in) :: n

	real(kreal), dimension(1:3) :: nA
	real(kreal), dimension(1:3,1:3) :: nnA
	real(kreal) :: nn

nn=n(1)**2+n(2)**2+n(3)**2
nA(1)=n(1)*A(1,1)+n(2)*A(2,1)+n(3)*A(3,1)
nA(2)=n(1)*A(1,2)+n(2)*A(2,2)+n(3)*A(3,2)
nA(3)=n(1)*A(1,3)+n(2)*A(2,3)+n(3)*A(3,3)

nnA(1,1)=n(1)*nA(1);nnA(1,2)=n(1)*nA(2);nnA(1,3)=n(1)*nA(3);
nnA(2,1)=n(2)*nA(1);nnA(2,2)=n(2)*nA(2);nnA(2,3)=n(2)*nA(3);
nnA(3,1)=n(3)*nA(1);nnA(3,2)=n(3)*nA(2);nnA(3,3)=n(3)*nA(3);
R=A-2.d0*nnA/nn
end subroutine reflection

subroutine reflection_2D (R,A,n)

	real(kreal), dimension(1:2,1:2),intent(in) :: A
	real(kreal), dimension(1:2,1:2),intent(out) :: R
	real(kreal), dimension(1:2),intent(in) :: n

	real(kreal), dimension(1:2) :: nA
	real(kreal), dimension(1:2,1:2) :: nnA
	real(kreal) :: nn

nn=n(1)**2+n(2)**2
nA(1)=n(1)*A(1,1)+n(2)*A(2,1)
nA(2)=n(1)*A(1,2)+n(2)*A(2,2)

nnA(1,1)=n(1)*nA(1);nnA(1,2)=n(1)*nA(2)
nnA(2,1)=n(2)*nA(1);nnA(2,2)=n(2)*nA(2)
R=A-2.d0*nnA/nn
end subroutine reflection_2D

! ! ----------------------------------------------------------------------------
! ! Numerical diagonalization of 3x3 matrcies
! ! Copyright (C) 2006  Joachim Kopp
! ! ----------------------------------------------------------------------------
! ! ----------------------------------------------------------------------------
! SUBROUTINE DSYEVJ3(A, W,faceIndex)
!   ! ----------------------------------------------------------------------------
!   ! Calculates the eigenvalues and normalized eigenvectors of a symmetric 3x3
!   ! matrix A using the Jacobi algorithm.
!   ! The upper triangular part of A is destroyed during the calculation,
!   ! the diagonal elements are read but not destroyed, and the lower
!   ! triangular elements are not referenced at all.
!   ! ----------------------------------------------------------------------------
!   ! Parameters:
!   !   A: The symmetric input matrix
!   !   Q: Storage buffer for eigenvectors
!   !   W: Storage buffer for eigenvalues
!   ! ----------------------------------------------------------------------------
!   !     .. Arguments ..
!   real(kreal),dimension(1:3,1:3):: A
!   real(kreal),dimension(1:3,1:3):: Q
! 	integer :: faceIndex
!   real(kreal),dimension(1:3):: W
!
!   integer(kint) ::          N=3
!
!   real(kreal):: SD, SO,S, C, T,G, H, Z, xf,THRESH
!   integer(kint)::          I, X, Y, R
!
!   !     Initialize Q to the identitity matrix
!   !     --- This loop can be omitted if only the eigenvalues are desired ---
!   do X = 1, N
!     Q(X,X) = 1.0D0
!     do Y = 1, X-1
!       Q(X, Y) = 0.0D0
!       Q(Y, X) = 0.0D0
!     enddo
!   enddo
!
!   !     Initialize W to diag(A)
!   do X = 1, N
!     W(X) = A(X, X)
!   enddo
!
!   !     Calculate SQR(tr(A))
!   SD = 0.0D0
!   DO X = 1, N
!     SD = SD + ABS(W(X))
!   enddo
!   SD = SD**2
!
!   !     Main iteration loop
!   DO I = 1, 50
!     !       Test for convergence
!     SO = 0.0D0
!     DO  X = 1, N
!       DO Y = X+1, N
!         SO = SO + ABS(A(X, Y))
!       enddo
!     enddo
!     IF (SO .EQ. 0.0D0) THEN
!       RETURN
!     END IF
!
!     IF (I .LT. 4) THEN
!       THRESH = 0.2D0 * SO / N**2
!     ELSE
!       THRESH = 0.0D0
!     END IF
!
!     !       Do sweep
!     DO X = 1, N
!       DO Y = X+1, N
!         G = 100.0D0 * ( ABS(A(X, Y)) )
!         IF ( I .GT. 4 .AND. ABS(W(X)) + G .EQ. ABS(W(X)) &
!         .AND. ABS(W(Y)) + G .EQ. ABS(W(Y)) ) THEN
!         A(X, Y) = 0.0D0
!       ELSE IF (ABS(A(X, Y)) .GT. THRESH) THEN
!         !             Calculate Jacobi transformation
!         H = W(Y) - W(X)
!         IF ( ABS(H) + G .EQ. ABS(H) ) THEN
!           T = A(X, Y) / H
!         ELSE
!           xf = 0.5D0 * H / A(X, Y)
!           IF (xf .LT. 0.0D0) THEN
!             T = -1.0D0* (SQRT(1.0D0 + xf**2) - xf)**(-1)
!           ELSE
!             T = 1.0D0* (SQRT(1.0D0 + xf**2) + xf)**(-1)
!           END IF
!         END IF
!
!         C = 1.0D0* SQRT( 1.0D0 + T**2 )**(-1)
!         S = T * C
!         Z = T * A(X, Y)
!
!         !             Apply Jacobi transformation
!         A(X, Y) = 0.0D0
!         W(X)    = W(X) - Z
!         W(Y)    = W(Y) + Z
!         DO  R = 1, X-1
!           T       = A(R, X)
!           A(R, X) = C * T - S * A(R, Y)
!           A(R, Y) = S * T + C * A(R, Y)
!         enddo
!         DO  R = X+1, Y-1
!           T       = A(X, R)
!           A(X, R) = C * T - S * A(R, Y)
!           A(R, Y) = S * T + C * A(R, Y)
!         enddo
!         DO  R = Y+1, N
!           T       = A(X, R)
!           A(X, R) = C * T - S * A(Y, R)
!           A(Y, R) = S * T + C * A(Y, R)
!         enddo
!
!         !             Update eigenvectors
!         !             --- This loop can be omitted if only the eigenvalues are desired ---
!        ! DO R = 1, N
!        !  T       = Q(R, X)
!        !   Q(R, X) = C * T - S * Q(R, Y)
!        !   Q(R, Y) = S * T + C * Q(R, Y)
!        ! enddo
!       END IF
!     enddo
!   enddo
! enddo
!
! !PRINT *, faceIndex,"DSYEVJ3: No convergence."
!
! END SUBROUTINE DSYEVJ3

!
!----------------
! private methods
!----------------
!

!> @brief Initializes a logger for the SsrfpackRemesh module
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
