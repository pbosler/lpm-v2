module RefinementModule

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
	logical(klog), pointer :: refineFlag(:) => null()
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
subroutine newPrivate(self, nMaxFaces )
	type(RefineSetup), intent(out) :: self
	integer(kint), intent(in) :: nMaxFaces
	
	if (.NOT. logInit ) call InitLogger(log, procRank)
	
	allocate(self%refineFlag(nMaxFaces))
	self%refineFlag = .FALSE.	
end subroutine

subroutine deletePrivate(self)
	type(RefineSetup), intent(inout) :: self
	if ( associated(self%refineFlag)) deallocate(self%refineFlag)
end subroutine

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