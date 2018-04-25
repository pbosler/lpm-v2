program EllipticVortexDriver

use NumberKindsModule
use OutputWriterModule
use LoggerModule
use ParticlesModule
use EdgesModule
use FacesModule
use PolyMesh2dModule
use FieldModule
use MPISetupModule
use PlaneGeomModule
use PlanarIncompressibleModule
use PlanarIncompressibleSolverModule
use RefinementModule
use BIVARRemeshModule

implicit none

include 'mpif.h'

! mesh variables
type(PlaneMeshIncompressible) :: plane 
integer(kint) :: initNest
integer(kint) :: maxNest
integer(kint) :: amrLimit
integer(kint) :: faceKind
integer(kint) :: meshSeed
real(kreal) :: meshRadius
integer(kint) :: testID

! AMR variables
type(RefineSetup) :: refinement
logical(klog) :: useAMR
logical(klog) :: doFlowMapRefine
integer(kint) :: nParticlesBefore
integer(kint) :: nParticlesAfter
real(kreal) :: circulationTol
real(kreal) :: flowMapVarTol
real(kreal) :: vortVarTol

! remeshing variables
integer(kint) :: remeshInterval
integer(kint) :: resetLagParamInterval
integer(kint) :: remeshCounter
type(PlaneMeshIncompressible) :: tempMesh
type(PlaneMeshIncompressible) :: refMesh

! timestepping variables
type(PlaneSolver) :: solver
real(kreal) :: t
real(kreal) :: dt
real(kreal) :: tfinal
integer(kint) :: timeJ
integer(kint) :: nTimesteps
real(kreal), allocatable, dimension(:) :: kineticEnergy, enstrophy

! i/o variables
character(len=MAX_STRING_LENGTH) :: outputDir
character(len=MAX_STRING_LENGTH) :: outputRoot
character(len=MAX_STRING_LENGTH) :: vtkFile
character(len=MAX_STRING_LENGTH) :: matlabFile
character(len=MAX_STRING_LENGTH) :: vtkRoot
character(len=MAX_STRING_LENGTH) :: meshString
integer(kint) :: frameOut
integer(kint) :: frameCounter

namelist /meshDefine/ faceKind, initNest, maxNest, amrLimit, meshRadius, testID, &
	 circulationTol, flowMapVarTol, vortVarTol
namelist /timestepping/ dt, tfinal, remeshInterval, resetLagParamInterval
namelist /fileIO/ outputDir, outputRoot, frameOut

! computing environment / general
type(Logger) :: exeLog
logical(klog), save :: logInit = .FALSE.
character(len=MAX_STRING_LENGTH) :: logString
integer(kint), parameter :: logLevel = DEBUG_LOGGING_LEVEL
character(len=28) :: logKey = "ellipticVortex"
integer(kint) :: mpiErrCode
real(kreal) :: timeStart, timeEnd
integer(kint) :: i

! vorticity profile
real(kreal), parameter :: lambda = 20.0_kreal
real(kreal), parameter :: R0 = 0.8_kreal
real(kreal), parameter :: q = 2.56085_kreal

!--------------------------------
!	initialize : setup computing environment
!--------------------------------
call MPI_INIT(mpiErrCode)
call MPI_COMM_SIZE(MPI_COMM_WORLD, numProcs, mpiErrCode)
call MPI_COMM_RANK(MPI_COMM_WORLD, procRank, mpiErrCode)

call InitLogger(exeLog, procRank)

timeStart = MPI_WTIME()

call ReadNamelistFile( procRank )

!
!	initialize mesh and spatial fields
!
call New( plane, initNest, maxNest, meshSeed, amrLimit, meshRadius)
if ( testID == 1 ) then
	call SetInitialVorticityOnMesh( plane, ellipticVortex1 )
else
	call SetInitialVorticityOnMesh( plane, ellipticVortex2 )
endif

useAMR = ( amrLimit > 0 )
if ( useAMR ) then

	doFlowMapRefine = .TRUE.
	
	call New(refinement, plane%mesh%faces%N_Max)
	
	call SetAbsoluteTolerances( plane, circulationTol, flowMapVarTol )
	call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logKey)//" circulation tol = ", circulationTol)
	call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logKey)//" flowMapVarTol = ", flowMapVarTol)
	
	do i = 1, amrLimit
		call IterateMeshRefinementOneVariableAndFlowMap( refinement, plane%mesh, plane%vorticity, &
			ScalarIntegralRefinement, circulationTol, "circulation refinement", flowMapVarTol, &
			nParticlesBefore, nParticlesAfter)
		if ( testID == 1 ) then
			call SetInitialVorticityOnMesh( plane, ellipticVortex1 )
		else
			call SetInitialVorticityOnMesh( plane, ellipticVortex2 )
		endif
	enddo
	
	call LoadBalance( plane%mpiParticles, plane%mesh%particles%N, numProcs)
	
	call Delete(refinement)
endif

!call AddTracers( plane, 1, [2])
!plane%tracers(1)%name = "initX"
!call StoreLagParamAsTracer( plane, 1 )
call SetStreamFunctionOnMesh( plane )
call SetVelocityOnMesh( plane )

!
!	output initial state
!
frameCounter = 0
t = 0.0_kreal
if ( procRank == 0 ) then
	call LogStats(plane, exeLog)
	
	call MeshSeedString( meshString, meshSeed)
	
	if ( useAMR ) then
		write( meshString, '(A,I1,A,I1,A)') trim(meshString)//"_AMR_", initNest, 'to', initNest+amrLimit, '_'
	else
		write( meshString, '(A,I1,A)') trim(meshString)//'_', initNest, '_'
	endif
	
	write(vtkRoot,'(4A)') trim(outputDir), '/vtkOut/', trim(outputRoot), trim(meshString)
	write(vtkFile,'(A,I0.4,A)') trim(vtkRoot), frameCounter, '.vtk'
	write(vtkRoot,'(4A)') trim(outputDir), '/vtkOut/', trim(outputRoot), trim(meshString)
	write(vtkFile,'(A,I0.4,A)') trim(vtkRoot), frameCounter, '.vtk'
	
	call OutputToVTK(plane, vtkFile )
	
	frameCounter = frameCounter + 1
	
	call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logkey)//" t = ", plane%mesh%t)
endif

!
!	initialize timestepping
!
call New(solver, plane)
nTimesteps = floor(tfinal/dt)
allocate(kineticEnergy(nTimesteps+1))
allocate(enstrophy(nTimesteps+1))
kineticEnergy(1) = TotalKE(plane)
enstrophy(1) = TotalEnstrophy(plane)
remeshCounter = 0

!--------------------------------
!	run : evolve the problem in time 
!--------------------------------

do timeJ = 0, nTimesteps - 1
	!
	! remesh
	!
	if ( mod( timeJ+1, remeshInterval) == 0 ) then
		remeshCounter = remeshCounter + 1
		!
		!	choose remeshing procedure
		!
		if ( remeshCounter < resetLagParamInterval ) then
			!
			!	remesh to t = 0
			!
			call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logKey)//" remeshing to t = 0, remeshCounter = ", remeshCounter)
			
			call New(tempMesh, initNest, maxNest, meshSeed, amrLimit, meshRadius)
!			call AddTracers( tempMesh, 1, [2])
			
			if ( testID == 1 ) then
				call LagrangianRemeshPlanarIncompressibleWithVorticityFunction( tempMesh, plane, ellipticVortex1, &
					ScalarIntegralRefinement, circulationTol, "circulation refinement", &
					ScalarVariationRefinement, vortVarTol, "vorticity variation refinement", &
					RefineFlowMapYN = doFlowMapRefine, &
					flowMapVarTol = flowMapVarTol )
			else
				call LagrangianRemeshPlanarIncompressibleWithVorticityFunction( tempMesh, plane, ellipticVortex2, &
					ScalarIntegralRefinement, circulationTol, "circulation refinement", &
					ScalarVariationRefinement, vortVarTol, "vorticity variation refinement", &
					RefineFlowMapYN = doFlowMapRefine, &
					flowMapVarTol = flowMapVarTol )
			endif

!			call StoreLagParamAsTracer( tempMesh, 1)
		
			call Copy( plane, tempMesh)
			call Delete(tempMesh)
		elseif ( remeshCounter == resetLagParamInterval ) then
			!
			!	remesh to t = 0, keep old mesh as reference defined for t = t_remesh
			!
			call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logKey)//" remeshing to t = 0, remeshCounter = ", remeshCounter)
			
			t = plane%mesh%t
			call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logKey)//" setting new ref_t = ", t)
			
			call New( refMesh, initNest, maxNest, meshSeed, amrLimit, meshRadius )
!			call AddTracers( refMesh, 1, [2])
!			refMesh%tracers(1)%name = "initX"
			
			if ( testID == 1 ) then
				call LagrangianRemeshPlanarIncompressibleWithVorticityFunction( refMesh, plane, ellipticVortex1, &
					ScalarIntegralRefinement, circulationTol, "circulation refinement", &
					ScalarVariationRefinement, vortVarTol, "vorticity variation refinement", &
					RefineFlowMapYN = doFlowMapRefine, &
					flowMapVarTol = flowMapVarTol )
			else
				call LagrangianRemeshPlanarIncompressibleWithVorticityFunction( refMesh, plane, ellipticVortex2, &
					ScalarIntegralRefinement, circulationTol, "circulation refinement", &
					ScalarVariationRefinement, vortVarTol, "vorticity variation refinement", &
					RefineFlowMapYN = doFlowMapRefine, &
					flowMapVarTol = flowMapVarTol )
			endif
			
!			call StoreLagParamAsTracer(refMesh, 1)
			
			call Copy(plane, refMesh)
			refMesh%mesh%t = t		
			
			call ResetLagrangianParameter(refMesh%mesh%particles)
			call ResetLagrangianParameter(plane%mesh%particles)
		elseif ( remeshCounter > resetLagParamInterval .AND. mod(remeshCounter, resetLagParamInterval) == 0 ) then
			!
			!	remesh to previous reference time, then keep current mesh as new reference
			!
			call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logKey)//" remeshing to t = ", refMesh%mesh%t)
			t = plane%mesh%t
			call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logKey)//" setting new ref_t = ", t)
				
			call New(tempMesh, initNest, maxNest, meshSeed, amrLimit, meshRadius )
!			call AddTracers( tempMesh, 1, [2])
!			tempMesh%tracers(1)%name = "initX"
			
			call LogMessage(exeLog, DEBUG_LOGGING_LEVEL, trim(logKey)//" calling remesh routine = ", t)
			
			call LagrangianRemeshPlanarIncompressibleToReferenceMesh( tempMesh, plane, refMesh, ScalarIntegralRefinement, &
					circulationTol, "circulation refinement", &
					ScalarVariationRefinement, vortVarTol, "vorticity variation refinement", &
					RefineFlowMapYN = doFlowMapRefine, &
					flowMapVarTol = flowMapVarTol )
			
			call LogMessage(exeLog, DEBUG_LOGGING_LEVEL, trim(logKey)//" returned from remesh = ", t)
			
			call Copy(plane, tempMesh)
			call Copy(refMesh, tempMesh)
			refMesh%mesh%t = t
			
!			call ResetLagrangianParameter(refMesh%mesh%particles)
!			call ResetLagrangianParameter(plane%mesh%particles)
			call Delete(tempMesh)
		else
			!
			!	remesh to reference time
			!
			call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logKey)//" remeshing to t = ", refMesh%mesh%t)
			call New(tempMesh, initNest, maxNest, meshSeed, amrLimit, meshRadius )
!			call AddTracers( tempMesh, 1, [2])
!			tempMesh%tracers(1)%name = "initX"
			
			call LagrangianRemeshPlanarIncompressibleToReferenceMesh( tempMesh, plane, refMesh, ScalarIntegralRefinement, &
					circulationTol, "circulation refinement", RefineFlowMapYN = doFlowMapRefine, &
					flowMapVarTol = flowMapVarTol)
					
			call Copy(plane, tempMesh)
			call Delete(tempMesh)
		endif
		
		
		call Delete(solver)
		call New(solver, plane)
	endif
	
	if ( timeJ == 2 ) call LogStats(plane, exeLog)
	
	!
	! advance time step
	!
	call Timestep(solver, plane, dt)
	plane%mesh%t = real( timeJ + 1, kreal) * dt
	
	kineticEnergy(timeJ + 2) = TotalKE(plane)
	enstrophy(timeJ + 2) = TotalEnstrophy(plane)
	
	if ( procRank == 0 .AND. mod(timeJ + 1, frameOut) == 0 ) then
		write(vtkFile,'(A,I0.4,A)') trim(vtkRoot), frameCounter, '.vtk'
		call OutputToVTK(plane, vtkFile)
		frameCounter = frameCounter + 1
		call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logKey)//" t = ", plane%mesh%t)
	endif
enddo

!
!	write t = tfinal output
!
if ( procRank == 0 ) then
	write(matlabFile, '(5A)') trim(outputDir), '/', trim(outputRoot), trim(meshString), '.m'
	open(unit=WRITE_UNIT_1, file=matlabFile, status='REPLACE', action='WRITE')
		write(WRITE_UNIT_1,'(A,F12.9,A,F12.6,A)') "t = 0:", dt,":", tfinal, ";"
		call WriteToMatlab(kineticEnergy, WRITE_UNIT_1, "kineticEnergy")
		call WriteToMatlab(enstrophy, WRITE_UNIT_1, "enstrophy")
	close(WRITE_UNIT_1)
endif

!--------------------------------
!	finalize : clean up
!--------------------------------

timeEnd = MPI_WTIME()
write(logstring,'(A,F12.2,A)') "PROGRAM COMPLETE : elapsed time ", timeEnd - timeStart, " seconds."
call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logKey)//" ", logstring)

deallocate(kineticEnergy)
deallocate(enstrophy)
call Delete(refMesh)
call Delete(solver)
call Delete(plane)

call MPI_FINALIZE(mpiErrCode)

contains

pure function ellipticVortex1( x, y )
	real(kreal) :: ellipticVortex1
	real(kreal), intent(in) :: x, y
	real(kreal) :: r, z
	
	r = sqrt( x*x + y*y / 4.0_kreal)
	if ( r <= R0 ) then
		z = r / R0
		ellipticVortex1 = lambda * (1.0_kreal - exp( - (q / z) * exp( 1.0_kreal / (z - 1.0_kreal) ) ) )
	else
		ellipticVortex1 = 0.0_kreal
	endif
end function

pure function ellipticVortex2( x, y )
	real(kreal) :: ellipticVortex2
	real(kreal), intent(in) :: x, y
	real(kreal) :: r, z
	
	r = sqrt( x*x + y*y / 4.0 )
	if ( r <= R0 ) then
		z = r / R0
		ellipticVortex2 = lambda * ( 1 - z ** 4)
	else
		ellipticVortex2 = 0.0_kreal
	endif
end function

subroutine ReadNamelistFile( rank )
	integer(kint), intent(in) :: rank
	!
	character(len=MAX_STRING_LENGTH) :: namelistFilename
	integer(kint), parameter :: initBcast_intSize = 8
	integer(kint), parameter :: initBcast_realSize = 6
	integer(kint), dimension(initBcast_intSize) :: bcastIntegers
	real(kreal), dimension(initBcast_realSize) :: bcastReals
	integer(kint) :: mpiErrCode, readStat
	
	if ( COMMAND_ARGUMENT_COUNT() /= 1 ) then
		call LogMessage(exeLog, ERROR_LOGGING_LEVEL, trim(logKey), " ERROR: expected namelist file as 1st argument.")
		stop
	endif
	
!namelist /meshDefine/ faceKind, initNest, maxNest, amrLimit, meshRadius, testID
!namelist /fileIO/ outputDir, outputRoot, frameOut
!namelist /timestepping/ dt, tfinal, remeshInterval, resetLagParamInterval	
	
	if ( rank == 0 ) then
		call GET_COMMAND_ARGUMENT(1, namelistFilename)
		
		open(unit=READ_UNIT, file=namelistFilename, status='OLD', action='READ', iostat=readStat)
			if ( readStat /= 0 ) then
				call LogMessage(exeLog, ERROR_LOGGING_LEVEL, trim(logKey), " ERROR: cannot read namelist file.")
				stop
			endif
		
			read(READ_UNIT, nml=meshDefine)
			rewind(READ_UNIT)
			read(READ_UNIT, nml=timestepping)
			rewind(READ_UNIT)
			read(READ_UNIT, nml=fileIO)
		close(READ_UNIT)
		
		if ( faceKind == 3 ) then
			meshSeed = TRI_HEX_SEED
		elseif ( faceKind == 4) then
			meshSeed = QUAD_RECT_SEED
		else
			call LogMessage(exeLog, WARNING_LOGGING_LEVEL, trim(logkey)//" ReadNamelistFile WARNING : ", &
				" invalid faceKind -- using triangles.")
			meshSeed = QUAD_RECT_SEED
		endif
		
		if ( testID < 1 .OR. testID > 2 ) then
			testID = 1
			call LogMessage(exeLog, WARNING_LOGGING_LEVEL, trim(logkey)//" ReadNamelistFile WARNING : ", &
				" invalid testID -- using vorticity profile 1 from Koumoutsakos 1997.")
		endif
		
		bcastIntegers(1) = initNest
		bcastIntegers(2) = maxNest
		bcastIntegers(3) = amrLimit
		bcastIntegers(4) = frameOut
		bcastIntegers(5) = remeshInterval
		bcastIntegers(6) = meshSeed
		bcastIntegers(7) = resetLagParamInterval
		bcastIntegers(8) = testID
		
		bcastReals(1) = dt
		bcastReals(2) = tfinal
		bcastReals(3) = circulationTol
		bcastReals(4) = flowMapVarTol
		bcastReals(5) = meshRadius
		bcastReals(6) = vortVarTol
	endif
	
	call MPI_BCAST(bcastIntegers, initBCAST_intSize, MPI_INTEGER, 0, MPI_COMM_WORLD, mpiErrCode)
	if ( mpiErrCode /= 0 ) then
		call LogMessage(exeLog, ERROR_LOGGING_LEVEL, trim(logKey)//" bcastIntegers, MPI_BCAST ERROR : ", mpiErrCode)
	endif
	
	call MPI_BCAST(bcastReals, initBCAST_realSize, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpiErrCode)
	if ( mpiErrCode /= 0 ) then
		call LogMessage(exeLog, ERROR_LOGGING_LEVEL, trim(logKey)//" bcastReals, MPI_BCAST ERROR : ", mpiErrCode)
	endif
	
	initNest = bcastIntegers(1)
	maxNest = bcastIntegers(2)
	amrLimit = bcastIntegers(3)
	frameOut = bcastIntegers(4)
	remeshInterval = bcastIntegers(5)
	meshSeed = bcastIntegers(6)
	resetLagParamInterval = bcastIntegers(7)
	testID = bcastIntegers(8)
	
	dt = bcastReals(1)
	tfinal = bcastReals(2)
	circulationTol = bcastReals(3)
	flowMapVarTol = bcastReals(4)
	meshRadius = bcastReals(5)
	vortVarTol = bcastReals(6)
end subroutine

subroutine StoreLagParamAsTracer( aPlane, tracerID )
	type(PlaneMeshIncompressible), intent(inout) :: aPlane
	integer(kint), intent(in) :: tracerID
	!
	integer(kint) :: i
	
	do i = 1, aPlane%mesh%particles%N
		aPlane%tracers(tracerID)%xComp(i) = aPlane%mesh%particles%x0(i)
		aPlane%tracers(tracerID)%yCOmp(i) = aPlane%mesh%particles%y0(i)
	enddo
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

end program