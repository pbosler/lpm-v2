module BIVARRemeshModule
!> @file BIVARRemesh.f90
!> Data structure and methods for remapping planar LPM data using the BIVAR package.
!> @author Peter Bosler, Sandia National Laboratories Center for Computing Research
!>
!>
!> @defgroup BIVARRemesh BIVARRemesh
!> Data structure and methods for remapping planar LPM data using the BIVAR package.
!> Uses the interfaces provided by @ref BIVARInterface.
!>
!> For references describing the BIVAR package, see the @ref BIVARInterface module's detailed description.
!>
!> @{
use NumberKindsModule
use OutputWriterModule
use UtilitiesModule
use LoggerModule
use ParticlesModule
use EdgesModule
use FacesModule
use PolyMesh2dModule
use FieldModule
use MPISetupModule
use PlaneGeomModule
use BIVARInterfaceModule
use PlanarIncompressibleModule
use PlanarSWEModule
use BetaPlaneMeshModule
use RefinementModule

implicit none

include 'mpif.h'

private
!
!----------------
! Module types and public declarations
!----------------
!
public DirectRemeshBetaPlane
public LagrangianRemeshBetaPlaneWithVorticityFunction
public LagrangianRemeshPlanarIncompressibleWithVorticityFunction
public LagrangianRemeshPlanarIncompressibleToReferenceMesh
public DirectRemeshPlanarSWE!, LagrangianRemeshPlanarSWE

!
!----------------
! Module interfaces
!----------------
!

!
!----------------
! Logging
!----------------
!
logical(klog), save :: logInit = .FALSE.
type(Logger) :: log
character(len=28), save :: logKey = 'BIVARRemesh'
integer(kint), parameter :: logLevel = DEBUG_LOGGING_LEVEL
character(len=MAX_STRING_LENGTH) :: logString

contains
!
!----------------
! public methods
!----------------
!

!> @brief Performs a remesh/remap of an LPM simulation using direct interpolation of all variables in a @ref BetaPlane mesh.
!>
!> @param[in] oldBetaPlane @ref BetaPlane source
!> @param[inout] newBetaPlane @ref BetaPlane target
!> @param[in] AMR .TRUE. if adaptive refinement will be used
!> @param[in] vortFlagFn1 Flag function for vorticity, must have same interface as refinementmodule::FlagFunction
!> @param[in] tol1 tolerance for first flag function
!> @param[in] desc1 description of first type of refinement
!> @param[in] flagFn2 Flag function for field2, must have same interface as refinementmodule::FlagFunction
!> @param[in] tol2 tolerance for second flag function
!> @param[in] desc2 description of second type of refinement
!> @param[in] field2 @ref Field to use for second type of refinement
subroutine DirectRemeshBetaPlane( oldBetaPlane, newBetaPlane, AMR, vortFlagFn1, tol1, desc1, flagFn2, tol2, desc2, field2)
	type(BetaPlaneMesh), intent(in) :: oldBetaPlane
	type(BetaPlaneMesh), intent(inout) :: newBetaPlane
	logical(klog), intent(in) :: AMR
	procedure(FlagFunction), optional :: vortFlagFn1
	real(kreal), intent(in), optional :: tol1
	character(len=*), intent(in), optional :: desc1
	procedure(FlagFunction), optional :: flagFn2
	real(kreal), intent(in), optional :: tol2
	character(len=*), intent(in), optional :: desc2
	type(Field), intent(inout), optional :: field2
	!
	integer(kint) :: i, j, k, nn
	real(kreal), dimension(2) :: x0, vecT
	type(RefineSetup) :: refine
	integer(kint) :: refineVariableCount
	integer(kint) :: nParticlesBefore, nParticlesAfter
	type(BIVARInterface) :: bivar

	if ( .NOT. logInit ) call InitLogger(log, procRank)

	call New(bivar, oldBetaPlane%mesh%particles )

	!
	!	Remesh to a uniform mesh
	!
	nn = newBetaPlane%mesh%particles%N
	newBetaPlane%relVort%N = nn
	newBetaPlane%absVort%N = nn
	if ( allocated(newBetaPlane%tracers) ) then
		do i = 1, size(newBetaPlane%tracers)
			newBetaPlane%tracers(i)%N = nn
		enddo
	endif

	call InterpolateScalar( newBetaPlane%relVort%scalar(1:nn), newBetaPlane%mesh%particles%x(1:nn), &
		newBetaPlane%mesh%particles%y(1:nn), bivar, oldBetaPlane%mesh%particles, oldBetaPlane%relVort )
	!call SetBIVARMD(bivar, 3)
	call InterpolateScalar( newBetaPlane%absVort%scalar(1:nn), newBetaPlane%mesh%particles%x(1:nn), &
		newBetaPlane%mesh%particles%y(1:nn), bivar, oldBetaPlane%mesh%particles, oldBetaPlane%absVort )
	call InterpolateLagParam( newBetaPlane%mesh%particles%x0(1:nn), newBetaPlane%mesh%particles%y0(1:nn), &
			newBetaPlane%mesh%particles%x(1:nn), newBetaPlane%mesh%particles%y(1:nn), bivar, oldBetaPlane%mesh%particles)

	if ( allocated(newBetaPlane%tracers) ) then
		do i = 1, size(newBetaPlane%tracers)
			if ( newBetaPlane%tracers(i)%nDim == 1 ) then
				call InterpolateScalar( newBetaPlane%tracers(i)%scalar(1:nn), newBetaPlane%mesh%particles%x(1:nn), &
					newBetaPlane%mesh%particles%y(1:nn), bivar, oldBetaPlane%mesh%particles, oldBetaPlane%tracers(i) )
			else
				call InterpolateVector( newBetaPlane%tracers(i)%xComp(1:nn), newBetaPlane%tracers(i)%yComp(1:nn), &
					newBetaPlane%mesh%particles%x(1:nn), newBetaPlane%mesh%particles%y(1:nn), bivar, &
					oldBetaPlane%mesh%particles, oldBetaPlane%tracers(i) )
			endif
		enddo
	endif

	!
	!	adaptive refinement
	!
	if ( AMR ) then
		refineVariableCount = 0
		if ( present( vortFlagFn1 ) .AND. ( present(tol1) .AND. present(desc1) ) ) then
			refineVariableCount = 1
			if ( ( present(flagFn2) .AND. present(field2)) .AND. ( present(tol2) .AND. present(desc2) ) ) then
				refineVariableCount = 2
			endif
		endif

		if ( refineVariableCount > 0 ) then

			call New(refine, newBetaPlane%mesh%faces%N_Max)

			do i = 1, newBetaPlane%mesh%amrLimit

				if ( refineVariableCount == 1 ) then
					call IterateMeshRefinementOneVariable( refine, newBetaPlane%mesh, newBetaPlane%relVort, &
						vortFlagFn1, tol1, desc1, nParticlesBefore, nParticlesAfter)
				elseif ( refineVariableCount == 2 ) then
					call IterateMeshRefinementTwoVariables( refine, newBetaPlane%mesh, newBetaPlane%relVort, &
						vortFlagFn1, tol1, desc1, field2, flagFn2, tol2, desc2, nParticlesBefore, nParticlesAfter)
				endif

				if ( nParticlesAfter > nParticlesBefore ) then
					!call SetBIVARMD( bivar, 2 )

					call InterpolateScalar( newBetaPlane%relVort%scalar(nParticlesBefore + 1 : nParticlesAfter), &
							newBetaPlane%mesh%particles%x(nParticlesBefore + 1 : nParticlesAfter), &
							newBetaPlane%mesh%particles%y(nParticlesBefore + 1 : nParticlesAfter), &
							bivar, oldBetaPlane%mesh%particles, oldBetaPlane%relVort)

					!call SetBIVARMD( bivar, 3 )

					call InterpolateScalar( newBetaPlane%absVort%scalar(nParticlesBefore + 1 : nParticlesAfter), &
							newBetaPlane%mesh%particles%x(nParticlesBefore + 1 : nParticlesAfter), &
							newBetaPlane%mesh%particles%y(nParticlesBefore + 1 : nParticlesAfter), &
							bivar, oldBetaPlane%mesh%particles, oldBetaPlane%absVort)
					call InterpolateLagParam( newBetaPlane%mesh%particles%x0(nParticlesBefore + 1 : nParticlesAfter), &
							newBetaPlane%mesh%particles%y0(nParticlesBefore + 1 : nParticlesAfter), &
							newBetaPlane%mesh%particles%x(nParticlesBefore + 1 : nParticlesAfter), &
							newBetaPlane%mesh%particles%y(nParticlesBefore + 1 : nParticlesAfter), &
							bivar, oldBetaPlane%mesh%particles )

					if ( allocated(newBetaPlane%tracers) ) then
						do j = 1, size(newBetaPlane%tracers)
							if ( newBetaPlane%tracers(j)%nDim == 1 ) then
								call InterpolateScalar( newBetaPlane%tracers(j)%scalar(nParticlesBefore + 1 : nParticlesAfter), &
									newBetaPlane%mesh%particles%x(nParticlesBefore + 1 : nParticlesAfter), &
									newBetaPlane%mesh%particles%y(nParticlesBefore + 1 : nParticlesAfter), &
									bivar, oldBetaPlane%mesh%particles, oldBetaPlane%tracers(j))
							else
								call InterpolateVector( newBetaPlane%tracers(j)%xComp(nParticlesBefore + 1 : nParticlesAfter), &
									newBetaPlane%tracers(j)%yComp(nParticlesBefore + 1 : nParticlesAfter), &
									newBetaPlane%mesh%particles%x(nParticlesBefore + 1 : nParticlesAfter), &
									newBetaPlane%mesh%particles%y(nParticlesBefore + 1 : nParticlesAfter), &
									bivar, oldBetaPlane%mesh%particles, oldBetaPlane%tracers(j))
							endif
						enddo
					endif
				endif
			enddo
			newBetaPlane%relVort%N = newBetaPlane%mesh%particles%N
			newBetaPlane%absVort%N = newBetaPlane%mesh%particles%N
			if ( allocated(newBetaPlane%tracers) ) then
				do i = 1, size(newBetaPlane%tracers)
					newBetaPlane%tracers(i)%N = newBetaPlane%mesh%particles%N
				enddo
			endif
			call LoadBalance(newBetaPlane%mpiParticles, newBetaPlane%mesh%particles%N, numPRocs)
			call Delete(refine)
		endif
	endif

	call SetVelocityOnMesh( newBetaPlane )
	call SetStreamFunctionsOnMesh( newBetaPlane )
	call Delete(bivar)
end subroutine

subroutine DirectRemeshPlanarSWE(oldPlane, newPlane, useAMR, vortFlagFn1, tol1, desc1, flagFn2, tol2, desc2, field2)
	type(SWEMesh), intent(in) :: oldPlane
	type(SWEMesh), intent(inout) :: newPlane
	logical(klog), intent(in) :: useAMR
	procedure(FlagFunction), optional :: vortFlagFn1
	real(kreal), optional, intent(in) :: tol1
	character(len=*), intent(in), optional :: desc1
	procedure(FlagFunction), optional :: flagFn2
	real(kreal), intent(in), optional :: tol2
	character(len=*), intent(in), optional :: desc2
	type(Field), intent(inout), optional :: field2
	!
	integer(kint) :: i, j, k, nn
	real(kreal), dimension(2) :: x0, vecT
	type(RefineSetup) :: refine
	integer(kint) :: nParticlesBefore, nParticlesAfter
	type(BIVARInterface) :: bivar

	if ( .NOT. logInit) call InitLogger(log, procRank)

	call LogMessage(log, DEBUG_LOGGING_LEVEL, trim(logKey)//" DirectRemeshPlanarSWE : ", "entering.")

	!
	!	Remesh to uniform new mesh
	!
	call New(bivar, oldPlane%mesh%particles)
	nn = newPlane%mesh%particles%N
	newPlane%relVort%N = nn
	newPlane%potVort%N = nn
	newPlane%divergence%N = nn
	newPlane%velocity%N = nn
	newPlane%h%N = nn
!	if ( allocated(newPlane%tracers)) then
!		do i = 1, size(newPlane%tracers)
!			newPlane%tracers(i)%N = nn
!		enddo
!	endif

	call SetBIVARMD(bivar, 1)

	call InterpolateLagParam(newPlane%mesh%particles%x0(1:nn), newPlane%mesh%particles%y0(1:nn), &
		newPlane%mesh%particles%x(1:nn), newPlane%mesh%particles%y(1:nn), bivar, oldPlane%mesh%particles)
	call LogMessage(log, DEBUG_LOGGING_LEVEL, trim(logKey)//" DirectRemeshPlanarSWE : ", "LagParam interpolation done.")

	call SetBIVARMD(bivar, 3)

	call InterpolateScalar(newPlane%relVort%scalar(1:nn), newPlane%mesh%particles%x(1:nn), &
		newPlane%mesh%particles%y(1:nn), bivar, oldPlane%mesh%particles, oldPlane%relVort)
	call LogMessage(log, DEBUG_LOGGING_LEVEL, trim(logKey)//" DirectRemeshPlanarSWE : ", "relvort interpolation done.")

	call InterpolateScalar(newPlane%potVort%scalar(1:nn), newPlane%mesh%particles%x(1:nn), &
		newPlane%mesh%particles%y(1:nn), bivar, oldPlane%mesh%particles, oldPlane%potVort)
	call LogMessage(log, DEBUG_LOGGING_LEVEL, trim(logKey)//" DirectRemeshPlanarSWE : ", "potvort interpolation done.")

	call InterpolateScalar(newPlane%divergence%scalar(1:nn), newPlane%mesh%particles%x(1:nn), &
		newPlane%mesh%particles%y(1:nn), bivar, oldPlane%mesh%particles, oldPlane%divergence)
	call LogMessage(log, DEBUG_LOGGING_LEVEL, trim(logKey)//" DirectRemeshPlanarSWE : ", "divergence interpolation done.")

	call InterpolateScalar(newPlane%h%scalar(1:nn), newPlane%mesh%particles%x(1:nn), &
		newPlane%mesh%particles%y(1:nn), bivar, oldPlane%mesh%particles, oldPlane%h)
	call LogMessage(log, DEBUG_LOGGING_LEVEL, trim(logKey)//" DirectRemeshPlanarSWE : ", "h interpolation done.")

!	if (allocated(newPlane%tracers)) then
!		do i = 1, size(newPlane%tracers)
!			if (newPlane%tracers(i)%nDim == 1) then
!				call InterpolateScalar(newPlane%tracers(i)%scalar(1:nn), newPlane%mesh%particles%x(1:nn), &
!					newPlane%mesh%particles%y(1:nn), bivar, oldPlane%mesh%particles, oldPlane%tracers(i))
!			else
!				call InterpolateVector(newPlane%tracers(i)%xComp(1:nn), newPlane%tracers(i)%yComp(1:nn), &
!					newPlane%mesh%particles%x(1:nn), newPlane%mesh%particles%y(1:nn), bivar, oldPlane%mesh%particles, &
!					oldPlane%tracers(i))
!			endif
!		enddo
!	endif

	if ( useAMR ) then
		call LogMessage(log, TRACE_LOGGING_LEVEL, trim(logkey)//" DirectRemeshPlanarSWE WARNING : ", "AMR not implemented yet.")
		call SetBIVARMD(bivar, 2)
	endif

	call SetVelocityOnMesh(newPlane)

	call Delete(bivar)
end subroutine

!> @brief Performs a remesh/remap of an LPM simulation using indirect interpolation of all variables in a @ref BetaPlane mesh.
!> Remaps to reference time t = 0.
!>
!> @param[in] oldBetaPlane @ref BetaPlane source
!> @param[inout] newBetaPlane @ref BetaPlane target
!> @param[in] AMR .TRUE. if adaptive refinement will be used
!> @param[in] relVortFn Vorticity distribution function, must have same interface as numberkindsmodule::scalarFnOf2DSpace
!> @param[in] flagFn1 Flag function for vorticity, must have same interface as refinementmodule::FlagFunction
!> @param[in] tol1 tolerance for first flag function
!> @param[in] desc1 description of first type of refinement
!> @param[in] flagFn2 Flag function for field2, must have same interface as refinementmodule::FlagFunction
!> @param[in] tol2 tolerance for second flag function
!> @param[in] desc2 description of second type of refinement
!> @param[in] RefineFLowMapYN True if refinement of the flow map will be used
!> @param[in] flowMapVarTol tolerance value for Lagrangian coordinate variation per face
!> @param[in] nLagTracers number of Lagrangian passive tracers (currently 0, 1, or 2 are only values allowed)
!> @param[in] tracerFn1 tracer distribution function, must have same interface as numberkindsmodule::scalarFnOf2DSpace
!> @param[in] tracerFn2 tracer distribution function, must have same interface as numberkindsmodule::scalarFnOf2DSpace
subroutine LagrangianRemeshBetaPlaneWithVorticityFunction( oldBetaPlane, newBetaPlane, AMR, relVortFn, &
			flagFn1, tol1, desc1, flagFn2, tol2, desc2, RefineFLowMapYN, flowMapVarTol, nLagTracers, tracerFn1, tracerFn2)
	type(BetaPlaneMesh), intent(in) :: oldBetaPlane
	type(BetaPlaneMesh), intent(inout) :: newBetaPlane
	logical(klog), intent(in) :: AMR
	procedure(scalarFnOf2DSpace) :: relVortFn
	procedure(FlagFunction) :: flagFn1
	real(kreal), intent(in) :: tol1
	character(len=*), intent(in) :: desc1
	procedure(FlagFunction), optional :: flagFn2
	real(kreal), intent(in), optional :: tol2
	character(len=*), intent(in), optional :: desc2
	logical(klog), intent(in), optional :: RefineFLowMapYN
	real(kreal), intent(in), optional :: flowMapVarTol
	integer(kint), intent(in), optional :: nLagTracers
	procedure(scalarFnOf2DSpace), optional :: tracerFn1
	procedure(scalarFnOf2DSpace), optional :: tracerFn2
	!
	integer(kint) :: nTracers
	logical(klog) :: refineVorticityTwice, doFlowMapRefinement
	integer(kint) :: i, j, k, nn
	real(kreal), dimension(2) :: x0, vecT
	type(RefineSetup) :: refine
	integer(kint) :: nParticlesBefore, nParticlesAfter
	real(kreal) :: zeta0
	type(BIVARInterface) :: bivar

	if ( .NOT. logInit ) call InitLogger(log, procRank)

	call New(bivar, oldBetaPlane%mesh%particles)

	! Determine how many (0, 1, or 2) Lagrangian scalar tracers exist
	nTracers = 0
	if ( allocated(oldBetaPlane%tracers) ) then
		if ( present(tracerFn1) .AND. ( nLagTracers >=1 .AND. newBetaPlane%tracers(1)%nDim == 1 ) ) then
			nTracers = 1
			if ( present(tracerFn2) .AND. (nLagTracers == 2 .AND. newBetaPlane%tracers(2)%nDim == 1 ) ) nTracers = 2
		endif
	endif

	!
	!	Remesh to a uniform mesh
	!
	nn = newBetaPlane%mesh%particles%N
	newBetaPlane%relVort%N = nn
	newBetaPlane%absVort%N = nn
	if ( allocated(newBetaPlane%tracers) ) then
		do i = 1, size(newBetaPlane%tracers)
			newBetaPlane%tracers(i)%N = nn
		enddo
	endif

	call InterpolateLagParam( newBetaPlane%mesh%particles%x0(1:nn), newBetaPlane%mesh%particles%y0(1:nn), &
			newBetaPlane%mesh%particles%x(1:nn), newBetaPlane%mesh%particles%y(1:nn), bivar, &
			oldBetaPlane%mesh%particles )
	do i = 1, nn
		zeta0 = relVortFn( newBetaPlane%mesh%particles%x0(i), newBetaPlane%mesh%particles%y0(i) )

		newBetaPlane%absVort%scalar(i) = zeta0 + newBetaPlane%f0 + newBetaPlane%beta * newBetaPlane%mesh%particles%y0(i)
		newBetaPlane%relVort%scalar(i) = zeta0 + newBetaPlane%beta * ( newBetaPlane%mesh%particles%y0(i) - &
			newBetaPlane%mesh%particles%y(i) )

		if ( nTracers == 1 ) then
			newBetaPlane%tracers(1)%scalar(i) = tracerFn1( newBetaPlane%mesh%particles%x0(i), &
														   newBetaPlane%mesh%particles%y0(i) )
		elseif ( nTracers == 2 ) then
			newBetaPlane%tracers(1)%scalar(i) = tracerFn1( newBetaPlane%mesh%particles%x0(i), &
														   newBetaPlane%mesh%particles%y0(i) )
			newBetaPlane%tracers(2)%scalar(i) = tracerFn2( newBetaPlane%mesh%particles%x0(i), &
														   newBetaPlane%mesh%particles%y0(i) )
		endif
	enddo

	if ( allocated( newBetaPlane%tracers ) ) then
		!call SetBIVARMD( bivar, 3)
		do i = nTracers + 1, size(newBetaPlane%tracers)
			if ( newBetaPlane%tracers(i)%nDim == 1 ) then
				call InterpolateScalar( newBetaPlane%tracers(i)%scalar(1:nn), newBetaPlane%mesh%particles%x(1:nn), &
						newBetaPlane%mesh%particles%y(1:nn), bivar, oldBetaPlane%mesh%particles, oldBetaPlane%tracers(i))
			else
				call InterpolateVector( newBetaPlane%tracers(i)%xComp(1:nn), newBetaPlane%tracers(i)%yComp(1:nn), &
						newBetaPlane%mesh%particles%x(1:nn), newBetaPlane%mesh%particles%y(1:nn), bivar, &
						oldBetaPlane%mesh%particles, oldBetaPlane%tracers(i))
			endif
		enddo
	endif
	call LogMessage(log, DEBUG_LOGGING_LEVEL, trim(logKey)//" lagRemesh : ", "uniform mesh ready.")

	!
	!	AMR
	!
	if ( AMR ) then
		refineVorticityTwice = ( present(flagFn2) .AND. (present(tol2) .AND. present(desc2)))
		doFlowMapRefinement = ( RefineFLowMapYN .AND. present(flowMapVarTol))

		call New(refine, newBetaPlane%mesh%faces%N_Max)

		do i = 1, newBetaPlane%mesh%amrLimit

			!call SetBIVARMD(bivar, 2)

			if ( refineVorticityTwice) then
				if ( doFlowMapRefinement ) then
					call IterateMeshRefinementTwoVariablesAndFlowMap( refine, newBetaPlane%mesh, newBetaPlane%relVort, &
							flagFn1, tol1, desc1, newBetaPlane%relVort, flagFn2, tol2, desc2, flowMapVarTol, &
							nParticlesBefore, nParticlesAfter)
				else
					call IterateMeshRefinementTwoVariables( refine, newBetaPlane%mesh, newBetaPlane%relVort, &
							flagFn1, tol1, desc1, newBetaPlane%relVort, flagFn2, tol2, desc2, &
							nParticlesBefore, nParticlesAfter)
				endif
			else
				if ( doFlowMapRefinement ) then
					call IterateMeshRefinementOneVariableAndFlowMap( refine, newBetaPlane%mesh, newBetaPlane%relVort, &
						flagFn1, tol1, desc1, flowMapVarTol, nParticlesBefore, nParticlesAfter)
				else
					call IterateMeshRefinementOneVariable( refine, newBetaPlane%mesh, newBetaPlane%relVort, flagFn1, &
						tol1, desc1, nParticlesBefore, nParticlesAfter)
				endif
			endif

			if ( nParticlesAfter > nParticlesBefore ) then
				call InterpolateLagParam( newBetaPlane%mesh%particles%x0(nParticlesBefore + 1 : nParticlesAfter), &
							newBetaPlane%mesh%particles%y0(nParticlesBefore + 1 : nParticlesAfter), &
							newBetaPlane%mesh%particles%x(nParticlesBefore + 1 : nParticlesAfter), &
							newBetaPlane%mesh%particles%y(nParticlesBefore + 1 : nParticlesAfter), &
							bivar, oldBetaPlane%mesh%particles )
				!call SetBIVARMD(bivar, 3)

				do j = nParticlesBefore + 1, nParticlesAfter
					zeta0 = relVortFn( newBetaPlane%mesh%particles%x0(j), newBetaPlane%mesh%particles%y0(j) )

					newBetaPlane%absVort%scalar(j) = zeta0 + newBetaPlane%f0 + &
						newBetaPlane%beta * newBetaPlane%mesh%particles%y0(j)
					newBetaPlane%relVort%scalar(j) = zeta0 + newBetaPlane%beta * &
						(newBetaPlane%mesh%particles%y0(j) - newBetaPlane%mesh%particles%y(j) )

					if ( nTracers == 1 ) then
						newBetaPlane%tracers(1)%scalar(j) = tracerFn1( newBetaPlane%mesh%particles%x0(j), &
																	   newBetaPlane%mesh%particles%y0(j) )
					elseif ( nTracers == 2 ) then
						newBetaPlane%tracers(1)%scalar(j) = tracerFn1( newBetaPlane%mesh%particles%x0(j), &
																	   newBetaPlane%mesh%particles%y0(j) )
						newBetaPlane%tracers(2)%scalar(j) = tracerFn2( newBetaPlane%mesh%particles%x0(j), &
																	   newBetaPlane%mesh%particles%y0(j) )
					endif
				enddo

				if ( allocated( newBetaPlane%tracers) ) then
					do j = nTracers + 1, size(newBetaPlane%tracers)
						if ( newBetaPlane%tracers(j)%nDim == 1 ) then
							call InterpolateScalar( newBetaPlane%tracers(j)%scalar(nParticlesBefore + 1 : nParticlesAfter), &
									newBetaPlane%mesh%particles%x(nParticlesBefore + 1 : nParticlesAfter), &
									newBetaPlane%mesh%particles%y(nParticlesBefore + 1 : nParticlesAfter), &
									bivar, oldBetaPlane%mesh%particles, oldBetaPlane%tracers(j) )
						else
							call InterpolateVector( newBetaPlane%tracers(j)%xComp(nParticlesBefore + 1 : nParticlesAfter), &
									newBetaPlane%tracers(j)%yComp(nParticlesBefore + 1 : nParticlesAfter), &
									newBetaPlane%mesh%particles%x(nParticlesBefore + 1 : nParticlesAfter), &
									newBetaPlane%mesh%particles%y(nParticlesBefore + 1 : nParticlesAfter), &
									bivar, oldBetaPlane%mesh%particles, oldBetaPlane%tracers(j) )
						endif
					enddo
				endif
			endif
		enddo
		newBetaPlane%relVort%N = newBetaPlane%mesh%particles%N
		newBetaPlane%absVort%N = newBetaPlane%mesh%particles%N
		if ( allocated(newBetaPlane%tracers) ) then
			do i = 1, size(newBetaPlane%tracers)
				newBetaPlane%tracers(i)%N = newBetaPlane%mesh%particles%N
			enddo
		endif

		call LoadBalance(newBetaPlane%mpiParticles, newBetaPlane%mesh%particles%N, numProcs)
		call Delete(refine)
	endif

	call SetVelocityOnMesh(newBetaPlane)
	call SetStreamFunctionsOnMesh(newBetaPlane)
	call LogMessage(log, DEBUG_LOGGING_LEVEL, trim(logKey)//" lagRemesh : ", "velocity and stream functions set on new mesh.")

	call Delete(bivar)
end subroutine


!> @brief Performs a remesh/remap of an LPM simulation using indirect interpolation of all variables in a @ref PlanarIncompressible mesh.
!> Remaps to reference time t = 0.
!>
!> @param[inout] newPlane @ref PlanarIncompressible target
!> @param[in] oldPlane @ref PlanarIncompressible source
!> @param[in] vortFn Vorticity distribution function, must have same interface as numberkindsmodule::scalarFnOf2DSpace
!> @param[in] flagFn1 Flag function for vorticity, must have same interface as refinementmodule::FlagFunction
!> @param[in] tol1 tolerance for first flag function
!> @param[in] desc1 description of first type of refinement
!> @param[in] flagFn2 Flag function for field2, must have same interface as refinementmodule::FlagFunction
!> @param[in] tol2 tolerance for second flag function
!> @param[in] desc2 description of second type of refinement
!> @param[in] RefineFLowMapYN True if refinement of the flow map will be used
!> @param[in] flowMapVarTol tolerance value for Lagrangian coordinate variation per face
!> @param[in] nLagTracers number of Lagrangian passive tracers (currently 0, 1, or 2 are only values allowed)
!> @param[in] tracerFn1 tracer distribution function, must have same interface as numberkindsmodule::scalarFnOf2DSpace
!> @param[in] tracerFn2 tracer distribution function, must have same interface as numberkindsmodule::scalarFnOf2DSpace
subroutine LagrangianRemeshPlanarIncompressibleWithVorticityFunction( newPlane, oldPlane, vortFn, flagFn1, tol1, desc1, &
				flagFn2, tol2, desc2, RefineFLowMapYN, flowMapVarTol, nLagTracers, tracerFn1, tracerFn2 )
	type(PlaneMeshIncompressible), intent(inout) :: newPlane
	type(PlaneMeshIncompressible), intent(in) :: oldPlane
	procedure(scalarFnOf2DSpace) :: vortFn
	procedure(FlagFunction) :: flagFn1
	real(kreal), intent(in) :: tol1
	character(len=*), intent(in) :: desc1
	procedure(FlagFunction), optional :: flagFn2
	real(kreal), intent(in), optional :: tol2
	character(len=*), intent(in), optional :: desc2
	logical(klog), intent(in), optional :: RefineFLowMapYN
	real(kreal), intent(in), optional :: flowMapVarTol
	integer(kint), intent(in), optional :: nLagTracers
	procedure(scalarFnOf2DSpace), optional :: tracerFn1
	procedure(scalarFnOf2DSpace), optional :: tracerFn2
	!
	integer(kint) :: nTracers
	logical(klog) :: refineVorticityTwice, doFlowMapRefinement
	integer(kint) :: i, j, nn
	real(kreal), dimension(2) :: vecT
	type(RefineSetup) :: refine
	integer(kint) :: nParticlesBefore, nParticlesAfter
	type(BIVARInterface) :: bivar

	if ( .NOT. logInit ) call InitLogger(log, procRank)

	call New(bivar, oldPlane%mesh%particles)

	! Determine how many (0, 1, or 2) Lagrangian scalar tracers exist
	nTracers = 0
	if ( allocated( oldPlane%tracers ) .AND. allocated(newPlane%tracers) ) then
		if ( present( tracerFn1) .AND. ( nLagTracers >= 1 .AND. newPlane%tracers(1)%nDIm == 1 ) ) then
			nTracers = 1
			if ( present(tracerFn2) .AND. ( nLagTracers == 2 .AND. newPlane%tracers(2)%nDim == 2 ) ) nTracers = 2
		endif
	endif

	!
	!	Remesh to a uniform mesh
	!
	nn = newPlane%mesh%particles%N
	newPlane%vorticity%N = nn
	if ( allocated(newPlane%tracers) ) then
		do i = 1, size(newPlane%tracers)
			newPlane%tracers(i)%n = nn
		enddo
	endif

	call InterpolateLagParam( newPlane%mesh%particles%x0(1:nn), newPlane%mesh%particles%y0(1:nn), &
					newPlane%mesh%particles%x(1:nn), newPlane%mesh%particles%y(1:nn), bivar, oldPlane%mesh%particles)
	do i = 1, nn
		newPlane%vorticity%scalar(i) = vortFn( newPlane%mesh%particles%x0(i), newPlane%mesh%particles%y0(i) )

		if ( nTracers == 1 ) then
			newPlane%tracers(1)%scalar(i) = tracerFn1( newPlane%mesh%particles%x0(i), newPlane%mesh%particles%y0(i) )
		elseif (nTracers == 2) then
			newPlane%tracers(1)%scalar(i) = tracerFn1( newPlane%mesh%particles%x0(i), newPlane%mesh%particles%y0(i) )
			newPlane%tracers(2)%scalar(i) = tracerFn2( newPlane%mesh%particles%x0(i), newPlane%mesh%particles%y0(i) )
		endif
	enddo

	if ( allocated(newPlane%tracers) ) then
		call SetBIVARMD(bivar, 3)
		do i = nTracers + 1, size(newPlane%tracers)
			if ( newPlane%tracers(i)%nDim == 1 ) then
				call InterpolateScalar( newPlane%tracers(i)%scalar(1:nn), newPlane%mesh%particles%x(1:nn), &
					newPlane%mesh%particles%y(1:nn), bivar, oldPlane%mesh%particles, oldPlane%tracers(i))
			else
				call InterpolateVector( newPlane%tracers(i)%xComp(1:nn), newPlane%tracers(i)%yComp(1:nn), &
					newPlane%mesh%particles%x(1:nn), newPlane%mesh%particles%y(1:nn), bivar, &
					oldPlane%mesh%particles, oldPlane%tracers(i) )
			endif
		enddo
	endif

	call LogMessage(log, DEBUG_LOGGING_LEVEL, trim(logKey)//" lagRemesh : ", "uniform mesh ready.")

	!
	! AMR
	!
	if ( newPlane%useAMR ) then
		refineVorticityTwice = ( present(flagFn2) .AND. (present(tol2) .AND. present(desc2)))
		doFlowMapRefinement = ( RefineFLowMapYN .AND. present(flowMapVarTol))

		call New(refine, newPlane%mesh%faces%N_Max)

		do i = 1, newPlane%mesh%amrLimit

			call SetBIVARMD( bivar, 1)

			if ( refineVorticityTwice ) then
				if ( doFlowMapRefinement ) then
					call IterateMeshRefinementTwoVariablesAndFlowMap( refine, newPlane%mesh, newPlane%vorticity, &
							flagFn1, tol1, desc1, newPlane%vorticity, flagFn2, tol2, desc2, flowMapVarTol, &
							nParticlesBefore, nParticlesAfter)
				else
					call IterateMeshRefinementTwoVariables( refine, newPlane%mesh, newPlane%vorticity, &
							flagFn1, tol1, desc1, newPlane%vorticity, flagFn2, tol2, desc2, &
							nParticlesBefore, nParticlesAfter)
				endif
			else
				if ( doFlowMapRefinement ) then
					call IterateMeshRefinementOneVariableAndFlowMap( refine, newPlane%mesh, newPlane%vorticity, &
						flagFn1, tol1, desc1, flowMapVarTol, nParticlesBefore, nParticlesAfter)
				else
					call IterateMeshRefinementOneVariable( refine, newPlane%mesh, newPlane%vorticity, flagFn1, &
						tol1, desc1, nParticlesBefore, nParticlesAfter)
				endif
			endif


			if ( nParticlesAfter > nParticlesBefore ) then
				write(logString,'(A, I4, A, I6, A)') "amr loop ", i, ": ", nParticlesAfter - nParticlesBefore, " particles added."
				call LogMessage(log, DEBUG_LOGGING_LEVEL, trim(logkey), trim(logString) )

				call InterpolateLagParam( newPlane%mesh%particles%x0(nParticlesBefore + 1 : nParticlesAfter), &
						newPlane%mesh%particles%y0(nParticlesBefore + 1 : nParticlesAfter), &
						newPlane%mesh%particles%x(nParticlesBefore + 1 : nParticlesAfter), &
						newPlane%mesh%particles%y(nParticlesBefore + 1 : nParticlesAfter), &
						bivar, oldPlane%mesh%particles)
				call SetBIVARMD(bivar, 3)

				do j = nParticlesBefore + 1, nParticlesAfter
					newPlane%vorticity%scalar(j) = vortFn( newPlane%mesh%particles%x0(j), newPlane%mesh%particles%y0(j))

					if ( nTracers == 1) then
						newPlane%tracers(1)%scalar(i) = tracerFn1( newPlane%mesh%particles%x0(j), &
																   newPlane%mesh%particles%y0(j))
					elseif ( nTracers == 2 ) then
						newPlane%tracers(1)%scalar(i) = tracerFn1( newPlane%mesh%particles%x0(j), &
																   newPlane%mesh%particles%y0(j))
						newPlane%tracers(2)%scalar(i) = tracerFn2( newPlane%mesh%particles%x0(j), &
																   newPlane%mesh%particles%y0(j))
					endif
				enddo


				if ( allocated(newPlane%tracers) ) then
					do j = ntracers + 1, size(newPlane%tracers)
						if ( newPlane%tracers(j)%nDim == 1) then
							call InterpolateScalar( newPlane%tracers(j)%scalar(nParticlesBefore + 1 : nParticlesAfter), &
									newPlane%mesh%particles%x(nParticlesBefore + 1 : nParticlesAfter), &
									newPlane%mesh%particles%y(nParticlesBefore + 1 : nParticlesAfter), &
									bivar, oldPlane%mesh%particles, oldPlane%tracers(j) )
						else
							call InterpolateVector( newPlane%tracers(j)%xComp(nParticlesBefore + 1 : nParticlesAfter), &
										newPlane%tracers(j)%yComp(nParticlesBefore + 1 : nParticlesAfter), &
										newPlane%mesh%particles%x(nParticlesBefore + 1 : nParticlesAfter), &
										newPlane%mesh%particles%y(nParticlesBefore + 1 : nParticlesAfter), &
										bivar, oldPlane%mesh%particles, oldPlane%tracers(j) )
						endif
					enddo
				endif
			else
				exit
			endif
		enddo

		nn = newPlane%mesh%particles%N
		newPlane%vorticity%N = nn
		if ( allocated(newPlane%tracers) ) then
			do i = 1, size(newPlane%tracers)
				newPlane%tracers(i)%N = nn
			enddo
		endif

		call LoadBalance(newPlane%mpiParticles, nn, numProcs)

		call Delete(refine)
	endif

	call SetVelocityOnMesh(newPlane)
	call SetStreamFunctionOnMesh(newPlane)

	call Delete(bivar)
end subroutine


!> @brief Performs a remesh/remap of an LPM simulation using indirect interpolation of all variables in a @ref PlanarIncompressible mesh.
!> Remaps to reference time t = t_{rm}, where t_{rm} is time of definition for a reference mesh (numerical solution).
!> Note that the reference mesh's Lagrangian parameter has been reset, so that x = x0, y = y0 at t = t_{rm}
!>
!> @param[inout] newPlane @ref PlanarIncompressible target
!> @param[in] oldPlane @ref PlanarIncompressible source
!> @param[in] refPlane @ref PlanarIncompressible reference
!> @param[in] flagFn1 Flag function for vorticity, must have same interface as refinementmodule::FlagFunction
!> @param[in] tol1 tolerance for first flag function
!> @param[in] desc1 description of first type of refinement
!> @param[in] flagFn2 Flag function for field2, must have same interface as refinementmodule::FlagFunction
!> @param[in] tol2 tolerance for second flag function
!> @param[in] desc2 description of second type of refinement
!> @param[in] RefineFLowMapYN True if refinement of the flow map will be used
!> @param[in] flowMapVarTol tolerance value for Lagrangian coordinate variation per face
subroutine LagrangianRemeshPlanarIncompressibleToReferenceMesh( newPlane, oldPlane, refPlane, flagFn1, tol1, desc1, &
				flagFn2, tol2, desc2, RefineFLowMapYN, flowMapVarTol )
	type(PlaneMeshIncompressible), intent(inout) :: newPlane
	type(PlaneMeshIncompressible), intent(in) :: oldPlane
	type(PlaneMeshIncompressible), intent(in) :: refPlane
	procedure(FlagFunction) :: flagFn1
	real(kreal), intent(in) :: tol1
	character(len=*), intent(in) :: desc1
	procedure(FlagFunction), optional :: flagFn2
	real(kreal), intent(in), optional :: tol2
	character(len=*), intent(in), optional :: desc2
	logical(klog), intent(in), optional :: RefineFLowMapYN
	real(kreal), intent(in), optional :: flowMapVarTol
	!
	integer(kint) :: i, j, nn
	type(BIVARInterface) :: bivar
	logical(klog) :: refineVorticityTwice
	logical(klog) :: doFlowMapRefinement
	type(RefineSetup) :: refine
	integer(kint) :: nParticlesBefore, nParticlesAfter

	!if ( .NOT. logInit ) call InitLogger(log, procRank)

	call New(bivar, oldPlane%mesh%particles)

	!
	! remesh to uniform new mesh
	!
	nn = newPlane%mesh%particles%N
	newPlane%vorticity%N = nn
	if ( allocated(newPlane%tracers) ) then
		do i = 1, size(newPlane%tracers)
			newPlane%tracers(i)%N = nn
		enddo
	endif

	! interpolate Lagrangian parameter from old mesh to new mesh
	call InterpolateLagParam( newPlane%mesh%particles%x0(1:nn), newPlane%mesh%particles%y0(1:nn), &
			newPlane%mesh%particles%x(1:nn), newPlane%mesh%particles%y(1:nn), bivar, oldPlane%mesh%particles)

	! set vorticity from reference mesh
	call InterpolateScalar( newPlane%vorticity%scalar(1:nn), newPlane%mesh%particles%x0(1:nn), &
			newPlane%mesh%particles%y0(1:nn), bivar, refPlane%mesh%particles, refPlane%vorticity)

	call SetBIVARMD( bivar, 3 )

	! set tracers from reference mesh
	if ( allocated(newPlane%tracers) .AND. allocated(oldPlane%tracers) ) then
		do i = 1, size(oldPlane%tracers)
			if ( oldPlane%tracers(i)%nDim == 1 ) then
				call InterpolateScalar( newPlane%tracers(i)%scalar(1:nn), newPlane%mesh%particles%x0(1:nn), &
						newPlane%mesh%particles%y0(1:nn), bivar, refPlane%mesh%particles, refPlane%tracers(i) )
			else
				call InterpolateVector( newPlane%tracers(i)%xComp(1:nn), newPlane%tracers(i)%yComp(1:nn), &
						newPlane%mesh%particles%x0(1:nn), newPlane%mesh%particles%y0(1:nn), bivar, &
						refPlane%mesh%particles, refPlane%tracers(i) )
			endif
		enddo
	endif

!	call LogMessage(log, DEBUG_LOGGING_LEVEL, trim(logKey)//" lagRemeshToRef ", "uniform mesh completed.")

	!
	!	AMR
	!
	if ( newPlane%useAMR ) then
		refineVorticityTwice = ( present(flagFn2) .AND. (present(tol2) .AND. present(desc2)))
		doFlowMapRefinement = ( RefineFLowMapYN .AND. present(flowMapVarTol))

		call New( refine, newPlane%mesh%faces%N_Max)

		do i = 1, newPlane%mesh%amrLimit
			call SetBIVARMD( bivar, 1)

!			call LogMessage(log, DEBUG_LOGGING_LEVEL, trim(logKey)//" lagRemeshToRef start amrLoop = ", i)

			if ( refineVorticityTwice ) then
				if ( doFlowMapRefinement ) then
					call IterateMeshRefinementTwoVariablesAndFlowMap( refine, newPlane%mesh, newPlane%vorticity, &
							flagFn1, tol1, desc1, newPlane%vorticity, flagFn2, tol2, desc2, flowMapVarTol, &
							nParticlesBefore, nParticlesAfter)
				else
					call IterateMeshRefinementTwoVariables( refine, newPlane%mesh, newPlane%vorticity, &
							flagFn1, tol1, desc1, newPlane%vorticity, flagFn2, tol2, desc2, &
							nParticlesBefore, nParticlesAfter)
				endif
			else
				if ( doFlowMapRefinement ) then
					call IterateMeshRefinementOneVariableAndFlowMap( refine, newPlane%mesh, newPlane%vorticity, &
						flagFn1, tol1, desc1, flowMapVarTol, nParticlesBefore, nParticlesAfter)
				else
					call IterateMeshRefinementOneVariable( refine, newPlane%mesh, newPlane%vorticity, flagFn1, &
						tol1, desc1, nParticlesBefore, nParticlesAfter)
				endif
			endif

!			call LogMessage(log, DEBUG_LOGGING_LEVEL, trim(logKey)//" lagRemeshToRef refinement complete = ", i)

			if ( nParticlesAfter > nParticlesBefore ) then

				call InterpolateLagParam( newPlane%mesh%particles%x0(nParticlesBefore + 1 : nParticlesAfter), &
						newPlane%mesh%particles%y0(nParticlesBefore + 1 : nParticlesAfter), &
						newPlane%mesh%particles%x(nParticlesBefore + 1 : nParticlesAfter), &
						newPlane%mesh%particles%y(nParticlesBefore + 1 : nParticlesAfter), &
						bivar, oldPlane%mesh%particles)

				call LogMessage(log, DEBUG_LOGGING_LEVEL, trim(logKey)//" lagRemeshToRef lagParam set = ", i)

				call InterpolateScalar( newPlane%vorticity%scalar(nParticlesBefore + 1 : nParticlesAfter), &
						newPlane%mesh%particles%x0(nParticlesBefore + 1 : nParticlesAfter), &
						newPlane%mesh%particles%y0(nParticlesBefore + 1 : nParticlesAfter), &
						bivar, refPlane%mesh%particles, refPlane%vorticity)
				call LogMessage(log, DEBUG_LOGGING_LEVEL, trim(logKey)//" lagRemeshToRef vorticity set = ", i)

				call SetBIVARMD(bivar, 3)

				if ( allocated(newPlane%tracers) .AND. allocated(oldPlane%tracers) ) then
					do j = 1, size(oldPlane%tracers)
						if ( oldPlane%tracers(j)%nDim == 1 ) then
							call InterpolateScalar( newPlane%tracers(j)%scalar(nParticlesBefore + 1 : nParticlesAfter), &
									newPlane%mesh%particles%x0(nParticlesBefore + 1 : nParticlesAfter), &
									newPlane%mesh%particles%y0(nParticlesBefore + 1 : nParticlesAfter), &
									bivar, refPlane%mesh%particles, refPlane%tracers(j))
						else
							call InterpolateVector( newPlane%tracers(j)%xComp(nParticlesBefore + 1 : nParticlesAfter), &
									newPlane%tracers(j)%yComp(nParticlesBefore + 1 : nParticlesAfter), &
									newPlane%mesh%particles%x0(nParticlesBefore + 1 : nParticlesAfter), &
									newPlane%mesh%particles%y0(nParticlesBefore + 1 : nParticlesAfter), &
									bivar, refPlane%mesh%particles, refPlane%tracers(j) )
						endif
!						call LogMessage(log, DEBUG_LOGGING_LEVEL, trim(logKey)//" lagRemeshToRef tracer set = ", i)
					enddo
				endif
			endif
		enddo

!		call LogMessage(log, DEBUG_LOGGING_LEVEL, trim(logKey)//" lagRemeshToRef : ", "interpolation complete.")

		nn = newPlane%mesh%particles%N
		newPlane%vorticity%N = nn
		if ( allocated(newPlane%tracers) ) then
			do i = 1, size(newPlane%tracers)
				newPlane%tracers(i)%N = nn
			enddo
		endif

!		call LogMessage(log, DEBUG_LOGGING_LEVEL, trim(logKey)//" lagRemeshToRef : ", " calling LoadBalance.")
		call LoadBalance(newPlane%mpiParticles, nn, numProcs)

		call Delete(refine)
	endif

	call SetVelocityOnMesh(newPlane)
	call SetStreamFunctionOnMesh(newPlane)

	call Delete(bivar)
end subroutine

!
!----------------
! private methods
!----------------
!

!> @brief Initializes a logger for the BIVARRemesh module
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
