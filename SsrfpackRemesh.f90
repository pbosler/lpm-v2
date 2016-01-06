module SSRFPACKRemeshModule
@file SsrfpackRemesh.f90
@author Peter Bosler, Sandia National Laboratories, Center for Computing Research
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
use RefinementModule
use MPISetupModule
use SSRFPACKInterfaceModule

implicit none

private
public BVERemeshSource, New, Delete
public DirectRemeshBVE
public LagrangianRemeshBVEWithVorticityFunction
public LagrangianRemeshBVEToReferenceMesh


!----------------
! types and module variables
!----------------
!
type BVERemeshSource
	type(DelaunayTriangulation) :: delTri
	type(SSRFPACKInterface) :: relVortSource
	type(SSRFPACKInterface) :: absVortSource
	type(SSRFPACKInterface) :: lagParamSource
	type(SSRFPACKInterface), dimension(:), allocatable :: tracerSource
	
	contains
		final :: deleteBVE
end type

integer(kint), save :: remeshCounter = 0

!
!----------------
! interfaces
!----------------
!
interface New
	module procedure newBVE
end interface

interface Delete
	module procedure deleteBVE
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
						call IterateMeshRefinementTwoVariables( refine, newSphere%mesh, newSphere%relVort, vortFlagFn1, &
								tol1, desc1, field2, flagFn2, tol2, desc2, nParticlesBefore, nParticlesAfter)
					else
						call IterateMeshRefinementTwoVariables( refine, newSphere%mesh, newSphere%relVort, vortFlagFn1, &
								tol1, desc1, newSphere%relVort, flagFn2, tol2, desc2, nParticlesBefore, nParticlesAfter)
					endif
				endif

				do j = nParticlesBefore + 1, nParticlesAfter
					lon = Longitude( newSphere%mesh%particles%x(j), newSphere%mesh%particles%y(j), &
									 newSphere%mesh%particles%z(j))
					lat = Latitude( newSphere%mesh%particles%x(j), newSphere%mesh%particles%y(j), &
									newSphere%mesh%particles%z(j))
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
	integer(kint) :: refinementType, nTracers
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

!
!----------------
! private methods
!----------------
!
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