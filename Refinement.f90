module RefinementModule
!> @file Refinement.f90
!> Data structure and methods for adaptively refining @ref PolyMesh2d objects.
!> @author Peter Bosler, Sandia National Laboratories Center for Computing Research
!> 
!>
!> @defgroup Refinement Refinement
!> Data structure and methods for adaptively refining @ref PolyMesh2d objects.
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

implicit none

private
public RefineSetup, New, Delete, FlagFunction
public ScalarIntegralRefinement, ScalarVariationRefinement
public IterateMeshRefinementForFlowMap
public IterateMeshRefinementOneVariable
public IterateMeshRefinementOneVariableAndFlowMap
public IterateMeshRefinementTwoVariables
public IterateMeshRefinementTwoVariablesAndFlowMap

!----------------
! types and module variables
!----------------
!
type RefineSetup
	logical(klog), pointer :: refineFlag(:) => null()  !< Face i will be refined if `refineSetup%%refineFlag(i) = .TRUE.`
	contains
		final :: deletePrivate
end type
!
!----------------
! interfaces
!----------------
!
interface New
	module procedure newPrivate
end interface

interface Delete
	module procedure deletePrivate
end interface

!> @brief Generic interface to a function that flags @ref Faces for refinement.
interface
	function FlagFunction( mesh, dataField, tolerance, faceIndex )
		use NumberKindsModule
		use PolyMesh2dModule
		use FieldModule
		logical(klog) :: FlagFunction
		type(PolyMesh2d), intent(in) :: mesh
		type(Field), intent(in) :: dataField
		real(kreal), intent(in) :: tolerance
		integer(kint), intent(in) :: faceIndex
	end function 
end interface

!
!----------------
! Logging
!----------------
!
logical(klog), save :: logInit = .FALSE.
type(Logger) :: log
character(len=28), save :: logKey = 'Refine'
character(len=MAX_STRING_LENGTH) :: logString
integer(kint), parameter :: logLevel = DEBUG_LOGGING_LEVEL

contains

!
!----------------
! public methods
!----------------
!
!> @brief Allocates memory for a RefineSetup object.
!> @param[out] self Target RefineSetup object
!> @param[in] nMaxFaces number of maximum faces allowed in memory
subroutine newPrivate(self, nMaxFaces )
	type(RefineSetup), intent(out) :: self
	integer(kint), intent(in) :: nMaxFaces
	
	if (.NOT. logInit ) call InitLogger(log, procRank)
	
	allocate(self%refineFlag(nMaxFaces))
	self%refineFlag = .FALSE.	
end subroutine

!> @brief Deletes and frees memory associated with a RefineSetup object
!> @param[inout] self Target RefineSetup object
subroutine deletePrivate(self)
	type(RefineSetup), intent(inout) :: self
	if ( associated(self%refineFlag)) deallocate(self%refineFlag)
end subroutine

!> @brief Calls the appropriate divide face subroutine on every flagged face.
!> Note:  New @ref Field values will need to be set after calling this routine for every new particle added here.
!> 
!> @param[in] self Target RefineSetup object
!> @param[inout] mesh @ref PolyMesh2d to be refined
subroutine DivideFlaggedFaces( self, mesh )
	type(RefineSetup), intent(in) :: self
	type(PolyMesh2d), intent(inout) :: mesh
	!
	integer(kint) :: i
	integer(kint) :: nFacesIn
	integer(kint) :: spaceLeft
	logical(klog) :: limitReached
	integer(kint) :: flagCount, refineCount
	
	nFacesIn = mesh%faces%N
	flagCount = count(self%refineFlag)
	spaceLeft = mesh%faces%N_Max - mesh%faces%N
	if ( spaceLeft / 4 < flagCount ) then
		call LogMessage(log, WARNING_LOGGING_LEVEL, trim(logKey)//" DivideFlaggedFaces WARNING : ", &
			"not enough memory to continue AMR.")
		return	
	endif
	
	call LogMessage(log, DEBUG_LOGGING_LEVEL, trim(logKey)//" entering :", " DivideFlaggedFaces")
	
	refineCount = 0
	limitReached = .FALSE.
	do i = 1, nFacesIn
		if ( self%refineFlag(i) ) then
			if ( mesh%initNest + mesh%amrLimit > CountParents(mesh%faces, i) ) then
				if ( mesh%faceKind == TRI_PANEL) then
					call DivideTriFace( mesh%faces, i, mesh%particles, mesh%edges)
				elseif ( mesh%faceKind == QUAD_PANEL ) then
					call DivideQuadFace( mesh%faces, i, mesh%particles, mesh%edges)
				endif
				refineCount = refineCount + 1
			else
				limitReached = .TRUE.
			endif
			self%refineFlag(i) = .FALSE.
		endif
	enddo
	if ( limitReached ) then
		call LogMessage(log, TRACE_LOGGING_LEVEL, trim(logKey)//" DivideFlaggedFaces REMARK : ", &
			"AMR limit reached. Some flagged faces were not divided.")
	endif
	write(logString,'(2(A,I6),A)') "refined ", refineCount, " of ", flagCount, " flagged faces."
	call LogMessage(log, TRACE_LOGGING_LEVEL, trim(logKey)//" DivideFlaggedFaces ", logString)
end subroutine

!> @brief Performs one iteration of mesh refinement using one criterion for one field variable. 
!> 
!> One iteration of mesh refinement is performed in the following steps:
!> * Loops over all active faces in a @ref PolyMesh2d, records faces flagged by the flag function
!> * Divides flagged faces
!> 
!> Following a call to this routine, all @ref Field objects that accompany the @ref PolyMesh2d that was refined need to 
!> be updated with new values at all new particles.
!> 
!> @param[inout] self target RefineSetup object
!> @param[inout] aMesh target @PolyMesh2d mesh to be refined
!> @param[in] aField @Field whose values will be used to flag faces
!> @param[in] flagFn Flag function, must have same interface as refinementmodule::FlagFunction
!> @param[in] tol tolerance value for use with Flag Function
!> @param[in] desc String describing the type of refinement
!> @param[out] nParticlesStart the number of particles in the mesh upon entry to this subroutine
!> @param[out] nParticlesEnd the number of particles in the mesh upon exit from this subroutine
subroutine IterateMeshRefinementOneVariable(self, aMesh, aField, flagFn, tol, desc, nParticlesStart, nParticlesEnd )
	type(RefineSetup), intent(inout) :: self
	type(PolyMesh2d), intent(inout) :: aMesh
	type(Field), intent(in) :: aField
	procedure(FlagFunction) :: flagFn
	real(kreal), intent(in) :: tol
	character(len=*), intent(in) :: desc
	integer(kint), intent(out) :: nParticlesStart
	integer(kint), intent(out) :: nParticlesEnd
	!
	integer(kint) :: counter
	integer(kint) :: i
	
	counter = 0
	nParticlesStart = aMesh%particles%N
	do i = 1, aMesh%faces%N
		if ( .NOT. aMesh%faces%hasChildren(i) ) then
			if ( flagFn(aMesh, aField, tol, i) ) then
				self%refineFlag(i) = .TRUE.
				counter = counter + 1			
			endif
		endif
	enddo
	
	write(logString, '(2A,I8,A)') trim(desc), ": flagged ", counter, " faces."
	call LogMessage(log, TRACE_LOGGING_LEVEL, trim(logKey)//" ", logString)
	
	call DivideFlaggedFaces( self, aMesh)
	nParticlesEnd = aMesh%particles%N
end subroutine

!> @brief Performs one iteration of mesh refinement using two criteria with two field variables. 
!> 
!> One iteration of mesh refinement is performed in the following steps:
!> * Loops over all active faces in a @ref PolyMesh2d, records faces flagged by the flag functions
!> * Divides flagged faces
!> 
!> Following a call to this routine, all @ref Field objects that accompany the @ref PolyMesh2d that was refined need to 
!> be updated with new values at all new particles.
!> 
!> @param[inout] self target RefineSetup object
!> @param[inout] aMesh target @PolyMesh2d mesh to be refined
!> @param[in] field1  First @Field whose values will be used to flag faces
!> @param[in] flagFn1 Flag function for use with field1, must have same interface as refinementmodule::FlagFunction
!> @param[in] tol1 tolerance value for use with first Flag Function
!> @param[in] desc1 String describing the first type of refinement
!> @param[in] field2 Second @Field whose values will be used to flag faces
!> @param[in] flagFn2 Flag function for use with field2, must have same interface as refinementmodule::FlagFunction
!> @param[in] tol2 tolerance value for use with second Flag Function
!> @param[in] desc2 String describing the second type of refinement
!> @param[out] nParticlesStart the number of particles in the mesh upon entry to this subroutine
!> @param[out] nParticlesEnd the number of particles in the mesh upon exit from this subroutine
subroutine IterateMeshRefinementTwoVariables( self, aMesh, field1, flagFn1, tol1, desc1, field2, flagFn2, tol2, desc2, &
	nParticlesStart, nParticlesEnd )
	type(RefineSetup), intent(inout) :: self
	type(PolyMesh2d), intent(inout) :: aMesh
	type(Field), intent(in) :: field1
	procedure(FlagFunction) :: flagFn1
	real(kreal), intent(in) :: tol1
	character(len=*), intent(in) :: desc1
	type(Field), intent(in) :: field2
	procedure(FlagFunction) :: flagFn2
	real(kreal), intent(in) :: tol2
	character(len=*), intent(in) :: desc2
	integer(kint), intent(out) :: nParticlesStart
	integer(kint), intent(out) :: nParticlesEnd
	!
	integer(kint), dimension(2) :: counters
	integer(kint) :: i
	
	counters = 0
	nParticlesStart = aMesh%particles%N
	do i = 1, aMesh%faces%N
		if ( .NOT. aMesh%faces%hasChildren(i) ) then
			if ( flagFn1(aMesh, field1, tol1, i) ) then
				self%refineFlag(i) = .TRUE.
				counters(1) = counters(1) + 1			
			endif
			if ( flagFn2(aMesh, field2, tol2, i ) ) then
				self%refineFlag(i) = .TRUE.
				counters(2) = counters(2) + 1
			endif
		endif
	enddo
	
	write(logString,'(2A,I8,A)') trim(desc1), ": flagged ", counters(1), " faces."
	call LogMessage(log, TRACE_LOGGING_LEVEL, trim(logKey)//" ", logString)
	
	write(logString,'(2A,I8,A)') trim(desc2), ": flagged ", counters(2), " faces."
	call LogMessage(log, TRACE_LOGGING_LEVEL, trim(logKey)//" ", logString)
	
	call DivideFlaggedFaces( self, aMesh)
	nParticlesEnd = aMesh%particles%N
end subroutine

!> @brief Performs one iteration of mesh refinement using two criteria with two field variables and flow map refinement. 
!> 
!> One iteration of mesh refinement is performed in the following steps:
!> * Loops over all active faces in a @ref PolyMesh2d, records faces flagged by the flag functions
!> * Divides flagged faces
!> 
!> Following a call to this routine, all @ref Field objects that accompany the @ref PolyMesh2d that was refined need to 
!> be updated with new values at all new particles.
!> 
!> @param[inout] self target RefineSetup object
!> @param[inout] aMesh target @PolyMesh2d mesh to be refined
!> @param[in] field1  First @Field whose values will be used to flag faces
!> @param[in] flagFn1 Flag function for use with field1, must have same interface as refinementmodule::FlagFunction
!> @param[in] tol1 tolerance value for use with first Flag Function
!> @param[in] desc1 String describing the first type of refinement
!> @param[in] field2 Second @Field whose values will be used to flag faces
!> @param[in] flagFn2 Flag function for use with field2, must have same interface as refinementmodule::FlagFunction
!> @param[in] tol2 tolerance value for use with second Flag Function
!> @param[in] desc2 String describing the second type of refinement
!> @param[in] flowMapVarTol tolerance for maximum Lagrangian coordinate variation per face
!> @param[out] nParticlesStart the number of particles in the mesh upon entry to this subroutine
!> @param[out] nParticlesEnd the number of particles in the mesh upon exit from this subroutine
subroutine IterateMeshRefinementTwoVariablesAndFlowMap( self, aMesh, field1, flagFn1, tol1, desc1, &
														field2, flagFn2, tol2, desc2, flowMapVarTol, &
														nParticlesStart, nParticlesEnd )
	type(RefineSetup), intent(inout) :: self
	type(PolyMesh2d), intent(inout) :: aMesh
	type(Field), intent(in) :: field1
	procedure(FlagFunction) :: flagFn1
	real(kreal), intent(in) :: tol1
	character(len=*), intent(in) :: desc1
	type(Field), intent(in) :: field2
	procedure(FlagFunction) :: flagFn2
	real(kreal), intent(in) :: tol2
	character(len=*), intent(in) :: desc2
	real(kreal), intent(in) :: flowMapVarTol
	integer(kint), intent(out) :: nParticlesStart
	integer(kint), intent(out) :: nParticlesEnd
	!
	integer(kint), dimension(3) :: counters
	integer(kint) :: i
	
!	call StartSection(log, trim(logKey)//" meshAMR status")
	
	counters = 0
	nParticlesStart = aMesh%particles%N
	do i = 1, aMesh%faces%N
		if ( .NOT. aMesh%faces%hasChildren(i) ) then
			if ( flagFn1(aMesh, field1, tol1, i) ) then
				self%refineFlag(i) = .TRUE.
				counters(1) = counters(1) + 1
			endif
			if ( flagFn2(aMesh, field2, tol2, i) ) then
				self%refineFlag(i) = .TRUE.
				counters(2) = counters(2) + 1
			endif
			if ( FlowMapVariationRefinement( aMesh, flowMapVarTol, i ) ) then
				self%refineFlag(i) = .TRUE.
				counters(3) = counters(3) + 1
			endif
		endif
	enddo
	
	write(logString,'(2A,I8,A)') trim(desc1), ": flagged ", counters(1), " faces."
	call LogMessage(log, TRACE_LOGGING_LEVEL, trim(logKey)//" ", logString)
	
	write(logString,'(2A,I8,A)') trim(desc2), ": flagged ", counters(2), " faces."
	call LogMessage(log, TRACE_LOGGING_LEVEL, trim(logKey)//" ", logString)
	
	write(logString,'(A,I8,A)') "flow map refinement : flagged ", counters(3), " faces."
	call LogMessage(log, TRACE_LOGGING_LEVEL, trim(logKey)//" ", logString)
	
	call DivideFlaggedFaces( self, aMesh)
	nParticlesEnd = aMesh%particles%N
	
!	call EndSection(log)
end subroutine

!> @brief Performs one iteration of mesh refinement using one criterion for one field variable and flow map refinement. 
!> 
!> One iteration of mesh refinement is performed in the following steps:
!> * Loops over all active faces in a @ref PolyMesh2d, records faces flagged by the flag function
!> * Divides flagged faces
!> 
!> Following a call to this routine, all @ref Field objects that accompany the @ref PolyMesh2d that was refined need to 
!> be updated with new values at all new particles.
!> 
!> @param[inout] self target RefineSetup object
!> @param[inout] aMesh target @PolyMesh2d mesh to be refined
!> @param[in] aField @Field whose values will be used to flag faces
!> @param[in] flagFn Flag function, must have same interface as refinementmodule::FlagFunction
!> @param[in] tol tolerance value for use with Flag Function
!> @param[in] desc String describing the type of refinement
!> @param[in] flowMapVarTol tolerance for maximum Lagrangian coordinate variation per face
!> @param[out] nParticlesStart the number of particles in the mesh upon entry to this subroutine
!> @param[out] nParticlesEnd the number of particles in the mesh upon exit from this subroutine
subroutine IterateMeshRefinementOneVariableAndFlowMap( self, aMesh, aField, flagFn, tol, desc, flowMapVarTol, &
	nParticlesStart, nParticlesEnd )
	type(RefineSetup), intent(inout) :: self
	type(PolyMesh2d), intent(inout) :: aMesh
	type(Field), intent(in) :: aField
	procedure(FlagFunction) :: flagFn
	real(kreal), intent(in) :: tol
	character(len=*), intent(in) :: desc
	real(kreal), intent(in) :: flowMapVarTol
	integer(kint), intent(out) :: nParticlesStart
	integer(kint), intent(out) :: nParticlesEnd
	!
	integer(kint), dimension(2) :: counters
	integer(kint) :: i
	
	call LogMessage(log, DEBUG_LOGGING_LEVEL, trim(logKey)//" entering ", "IterateMeshRefinementOneVariableAndFlowMap...")
	counters = 0
	call StartSection(log, trim(logKey)//" mesh status")
	call LogStats(aMesh, log)
	nParticlesStart = aMesh%particles%N
	do i = 1, aMesh%faces%N
		if ( .NOT. aMesh%faces%hasChildren(i) ) then
			!call LogMessage(log, DEBUG_LOGGING_LEVEL, trim(logKey)//" crit 1 face i = ", i)
			if ( flagFn( aMesh, aField, tol, i ) ) then
				self%refineFlag(i) = .TRUE.
				counters(1) = counters(1) + 1
			endif
			!call LogMessage(log, DEBUG_LOGGING_LEVEL, trim(logKey)//" crit 2 face i = ", i)
			if ( FlowMapVariationRefinement( aMesh, flowMapVarTol, i) ) then
				self%refineFlag(i) = .TRUE.
				counters(2) = counters(2) + 1
			endif
		endif
	enddo
	
	write(logString,'(2A,I8,A)') trim(desc), ": flagged ", counters(1), " faces."
	call LogMessage(log, TRACE_LOGGING_LEVEL, trim(logKey)//" ", logString)
	
	write(logString,'(A,I8,A)') "flow map refinement : flagged ", counters(2), " faces."
	call LogMessage(log, TRACE_LOGGING_LEVEL, trim(logKey)//" ", logString)
	
	call DivideFlaggedFaces( self, aMesh)
	nParticlesEnd = aMesh%particles%N
	call EndSection(log)
end subroutine

!> @brief Performs one iteration of mesh refinement using flow map refinement. 
!> 
!> One iteration of mesh refinement is performed in the following steps:
!> * Loops over all active faces in a @ref PolyMesh2d, records faces flagged by the flag function
!> * Divides flagged faces
!> 
!> Following a call to this routine, all @ref Field objects that accompany the @ref PolyMesh2d that was refined need to 
!> be updated with new values at all new particles.
!> 
!> @param[inout] self target RefineSetup object
!> @param[inout] aMesh target @PolyMesh2d mesh to be refined
!> @param[in] flowMapVarTol tolerance for maximum Lagrangian coordinate variation per face
!> @param[out] nParticlesStart the number of particles in the mesh upon entry to this subroutine
!> @param[out] nParticlesEnd the number of particles in the mesh upon exit from this subroutine
subroutine IterateMeshRefinementForFlowMap( self, aMesh, flowMapVarTol, nParticlesStart, nParticlesEnd)
	type(RefineSetup), intent(inout) :: self
	type(PolyMesh2d), intent(inout) :: aMesh
	real(kreal), intent(in) :: flowMapVarTol
	integer(kint), intent(out) :: nParticlesStart
	integer(kint), intent(out) :: nParticlesEnd
	!
	integer(kint) :: counter
	integer(kint) :: i
	
	counter = 0
	nParticlesStart = aMesh%particles%N
	do i = 1, aMesh%faces%N
		if ( .NOT. aMesh%faces%hasChildren(i) ) then
			if ( FlowMapVariationRefinement( aMesh, flowMapVarTol, i ) ) then
				self%refineFlag(i) = .TRUE.
				counter = counter + 1
			endif
		endif
	enddo
	
	write(logString,'(A,I8,A)') "flow map refinement : flagged ", counter, " faces."
	call LogMessage(log, TRACE_LOGGING_LEVEL, trim(logKey)//" ", logString)
	
	call DivideFlaggedFaces(self, aMesh)
	nParticlesEnd = aMesh%particles%N
end subroutine

!> @brief Example Flag Function to refine panels whose scalar @ref Field integrals exceed a given tolerance.
!> 
!> Scalar integrals on face i are approximated using the midpoint rule 
!> @f$ \int_{F_i} f(x)\,dx \approx @f$ `field%%scalar(i) * field%%area(i)`.
!> If this value exceeds the tolerance, this function returns .TRUE.
!>
!> @param[in] mesh @ref PolyMesh2d 
!> @param[in] dataField scalar @ref Field 
!> @param[in] tol tolerance value
!> @param[in] faceIndex face to flag, if tolerance is exceeded
!> @return Flag value
function ScalarIntegralRefinement( mesh, dataField, tol, faceIndex )
	logical(klog) :: ScalarIntegralRefinement
	type(PolyMesh2d), intent(in) :: mesh
	type(Field), intent(in) :: dataField
	real(kreal), intent(in) :: tol
	integer(kint), intent(in) :: faceIndex
	!
	integer(kint) :: pIndex
	pIndex = mesh%faces%centerParticle(faceIndex)
	ScalarIntegralRefinement = ( abs(dataField%scalar(pIndex)) * mesh%particles%area(pIndex) > tol )
end function

!> @brief Example Flag Function to refine panels whose scalar @Field variation per face exceed a given tolerance.
!> 
!> The scalar variation scalar f(x) of face F_i is computed as @f$ \max_{F_i} f(x) - \min_{F_i} f(x) @f$,
!> where the max and min are simply the maximum and minimum values of the scalar field carried by all particles associated
!> with face F_i (center and vertices).
!>
!> @param[in] mesh @PolyMesh2d 
!> @param[in] dataField scalar @Field 
!> @param[in] tol tolerance value
!> @param[in] faceIndex face to flag, if tolerance is exceeded
!> @return Flag value
function ScalarVariationRefinement( mesh, dataField, tol, faceIndex )
	logical(klog) :: ScalarVariationRefinement
	type(PolyMesh2d), intent(in) :: mesh
	type(Field), intent(in) :: dataField
	real(kreal), intent(in) :: tol
	integer(kint), intent(in) :: faceIndex
	!
	integer(kint) :: pIndex, i
	real(kreal) :: maxScalar, minScalar
	type(STDIntVector) :: faceVerts
	
	pIndex = mesh%faces%centerParticle(faceIndex)
	call CCWVerticesAroundFace( mesh, faceVerts, faceindex)
	maxScalar = dataField%scalar(pIndex)
	minScalar = maxScalar
	do i = 1, faceVerts%N
		if ( dataField%scalar( faceVerts%int(i) ) > maxScalar ) maxScalar = dataField%scalar( faceVerts%int(i))
		if ( dataField%scalar( faceVerts%int(i) ) < minScalar ) minScalar = dataField%scalar( faceVerts%int(i))
	enddo
	ScalarVariationRefinement = ( maxScalar - minScalar > tol )
end function

!> @brief Flag Function to refine panels whose Lagrangian coordinate variation per face exceed a given tolerance.
!> 
!> The variation of the Lagrangian coordinates in face F_j is computed as @f$ \sum_{i=1}^3 \left(\max_{F_j} x0_i - \min_{F_j} x0_i\right) @f$,
!> where the max and min are simply the maximum and minimum values of the Lagrangian coordinate components carried by all particles associated
!> with face F_j (center and vertices).
!>
!> @param[in] mesh @PolyMesh2d 
!> @param[in] tol tolerance value
!> @param[in] faceIndex face to flag, if tolerance is exceeded
!> @return Flag value
function FlowMapVariationRefinement( mesh, tol, faceIndex )
	logical(klog) :: FlowMapVariationRefinement
	type(PolyMesh2d), intent(in) :: mesh
	real(kreal), intent(in) :: tol
	integer(kint), intent(in) :: faceIndex
	!
	integer(kint) :: pIndex, i
	real(kreal), dimension(3) :: maxX0, minX0, x0i
	type(STDIntVector) :: faceVerts
	
	pIndex = mesh%faces%centerParticle(faceIndex)
	!call CCWVerticesAroundFace(mesh, faceVerts, faceIndex)
	maxX0 = LagCoord(mesh%particles, pIndex)
	minX0 = LagCoord(mesh%particles, pIndex)
	do i = 1, mesh%faceKind
		!x0i = LagCoord(mesh%particles, faceVerts%int(i) )
		x0i = LagCoord(mesh%particles, mesh%faces%vertices(i,faceIndex))
		if ( x0i(1) < minX0(1) ) minX0(1) = x0i(1)
		if ( x0i(2) < minX0(2) ) minX0(2) = x0i(2)
		if ( x0i(3) < minX0(3) ) minX0(3) = x0i(3)
		if ( x0i(1) > maxX0(1) ) maxX0(1) = x0i(1)
		if ( x0i(2) > maxX0(2) ) maxX0(2) = x0i(2)
		if ( x0i(3) > maxX0(3) ) maxX0(3) = x0i(3)
	enddo
	FlowMapVariationRefinement = ( sum(maxX0 - minX0) > tol )
end function 
!
!----------------
! private methods
!----------------
!
!> @brief Initializes a logger for the Refinement module
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