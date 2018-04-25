program BVERossbyHaurwitz54

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
use SphereBVESolverModule
use RefinementModule
use SSRFPACKRemeshModule

implicit none

include 'mpif.h'

! mesh variables
type(BVEMesh) :: sphere
type(RefineSetup) :: refinement
integer(kint) :: initNest
integer(kint) :: maxNest
integer(kint) :: faceKind
integer(kint) :: meshSeed
integer(kint) :: amrLimit
real(kreal), parameter :: radius = 1.0_kreal
logical(klog) :: continueAMR
logical(klog) :: AMR
integer(kint) :: nParticlesBefore
integer(kint) :: nParticlesAfter
real(kreal) :: flowMapVarTol
real(kreal) :: circulationTol
namelist /meshDefine/ faceKind, initNest, amrLimit

! test case variables
real(kreal), allocatable, dimension(:) :: kineticEnergy, enstrophy
real(kreal) :: rotationRate
real(kreal) :: rhWaveAmplitude
real(kreal) :: zonalWind
namelist /rhWave/ rotationRate, rhWaveAmplitude, flowMapVarTol, circulationTol, zonalWind
!real(kreal) :: legendre54
!real(kreal) :: RossbyHaurwitz54Vorticity
!external legendre54, RossbyHaurwitz54Vorticity

! remeshing
type(BVERemeshSource) :: remesh
integer(kint) :: remeshInterval
integer(kint) :: remeshCounter
type(BVEMesh) :: tempSphere
type(BVEMesh) :: refSphere
real(kreal) :: refTime

! timestepping
type(BVESolver) :: solver
real(kreal) :: dt
real(kreal) :: t
real(kreal) :: tfinal
integer(kint) :: nTimesteps
integer(kint) :: timeJ
integer(kint) :: resetLagParamInterval
namelist /timestepping/ dt, tfinal, remeshInterval, resetLagParamInterval

! i/o
character(len=MAX_STRING_LENGTH) :: outputDir
character(len=MAX_STRING_LENGTH) :: outputRoot
character(len=MAX_STRING_LENGTH) :: vtkFile
character(len=MAX_STRING_LENGTH) :: matlabFile
character(len=MAX_STRING_LENGTH) :: vtkRoot
character(len=MAX_STRING_LENGTH) :: meshString
integer(kint) :: frameOut
integer(kint) :: frameCounter
namelist /fileIO/ outputDir, outputRoot, frameOut

! computing environment / general
type(Logger) :: exeLog
character(len=MAX_STRING_LENGTH) :: logString
integer(kint), parameter :: logLevel = DEBUG_LOGGING_LEVEL
character(len=28) :: logKey = "RossbyHaurwitz54"
integer(kint) :: mpiErrCode
real(kreal) :: timeStart, timeEnd
integer(kint) :: i
real(kreal), dimension(3) :: vec3

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
t = 0.0_kreal
call New( sphere, meshSeed, initNest, maxNest, amrLimit, radius, rotationRate )
call SetInitialVorticityOnMesh( sphere, RossbyHaurwitz54Vorticity)
call LogMessage(exeLog, DEBUG_LOGGING_LEVEL, trim(logKey)//" ", "uniform mesh ready...")

AMR = (amrLimit > 0 )
if ( AMR ) then
	call LogMessage(exeLog, DEBUG_LOGGING_LEVEL, trim(logKey)//" ", "starting adaptive refinement...")
	call New(refinement, sphere%mesh%faces%N_Max)
	call SetAbsoluteTolerances( sphere, circulationTol, flowMapVarTol)
	
	do i = 1, amrLimit
		write(logString,"(A,I2)") 'amr loop ', i
		call StartSection(exeLog, trim(logString) )
		call IterateMeshRefinementOneVariableAndFlowMap( refinement, sphere%mesh, sphere%relVort, ScalarIntegralRefinement, &
			circulationTol, "circulation refinement", flowMapVarTol, nParticlesBefore, nParticlesAfter)
		call SetInitialVorticityOnMesh( sphere, RossbyHaurwitz54Vorticity)
		call EndSection(exeLog)
	enddo
	
	call Delete(refinement)
	call LoadBalance(sphere%mpiParticles, sphere%mesh%particles%N, numProcs)
endif

call SetStreamFunctionsOnMesh(sphere)
call SetVelocityOnMesh(sphere)

call AddTracers(sphere, 2, [1, 1])
sphere%tracers(1)%name = "initLat"
sphere%tracers(1)%units = "radians"
sphere%tracers(2)%name = "lat(x0)"
sphere%tracers(2)%units = "radians"
do i = 1, sphere%mesh%particles%N
	call InsertScalarToField(sphere%tracers(1), Latitude( sphere%mesh%particles%x0(i), sphere%mesh%particles%y0(i), &
		sphere%mesh%particles%z0(i)) )
	call InsertScalarToField(sphere%tracers(2), Latitude( sphere%mesh%particles%x0(i), sphere%mesh%particles%y0(i), &
		sphere%mesh%particles%z0(i)))
enddo

call LogMessage(exeLog, DEBUG_LOGGING_LEVEL, trim(logKey)//" ", &
	"initial setup complete; writing initial mesh to data file...")
frameCounter = 0
remeshCounter = 0
if ( procRank == 0 ) then
	call LogStats(sphere, exeLog)
	if (meshSeed == ICOS_TRI_SPHERE_SEED) then
		if ( initNest == maxNest) then
			write(meshString, '(A,I1,A)') '_icosTri', initNest, '_'
		else
			write(meshString, '(2(A,I1),A)') '_icosTriAMR', initNest, 'to', maxNest, '_'
		endif
	elseif (meshSeed == CUBED_SPHERE_SEED ) then
		if ( initNest == maxNest ) then
			write(meshString, '(A,I1,A)') '_cubedSphere', initNest, '_'
		else
			write(meshString, '(2(A,I1),A)') '_cubedSphereAMR', initNest, 'to', maxNest, '_'
		endif
	endif
	
	write(vtkRoot,'(4A)') trim(outputDir), '/vtkOut/', trim(outputRoot), trim(meshString)
	write(vtkFile,'(A,I0.4,A)') trim(vtkRoot), frameCounter, '.vtk'
	
	call OutputToVTK(sphere, vtkFile)
	frameCounter = frameCounter + 1
	
	call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logkey)//" t = ", t)
endif

!
!	initialize time stepping
!
call New(solver, sphere)
nTimesteps = floor(tfinal / dt)
allocate(kineticEnergy(nTimesteps+1))
allocate(enstrophy(nTimesteps+1))
kineticEnergy(1) = TotalKE(sphere)
enstrophy(1) = TotalEnstrophy(sphere)
refTime = 0.0_kreal

!--------------------------------
!	run : evolve the problem in time 
!--------------------------------
call LogMessage(exeLog, DEBUG_LOGGING_LEVEL, trim(logkey)//" ", "starting timestepping loop.")
do timeJ = 0, nTimesteps - 1
	!
	!	remesh if necessary
	!
	if ( mod(timeJ+1, remeshInterval) == 0 ) then
		remeshCounter = remeshCounter + 1
		!
		!	choose remeshing procedure
		!
		if ( remeshCounter < resetLagParamInterval ) then
			!
			!	remesh to t = 0
			!
			call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logKey)//" remeshing to t = 0, remeshCounter = ", remeshCounter)
			
			call New(remesh, sphere)
			
			call New( tempSphere, meshSeed, initNest, maxNest, amrLimit, radius, rotationRate )
			call AddTracers(tempSphere, 2, [1, 1])
			tempsphere%tracers(1)%name = "initLat"
			tempsphere%tracers(1)%units = "radians"
			tempsphere%tracers(2)%name = "lat(x0)"
			tempsphere%tracers(2)%units = "radians"
			
			call LagrangianRemeshBVEWithVorticityFunction( remesh, sphere, tempSphere, AMR, RossbyHaurwitz54Vorticity, &
				ScalarIntegralRefinement, circulationTol, "circulation refinement", &
				RefineFlowMapYN = .TRUE., flowMapVarTol = flowMapVarTol, nLagTracers = 1, tracerFn1 = InitLatTracer)
			
!			call LagrangianRemeshBVEWithVorticityFunction(remesh, sphere, tempSphere, AMR, RossbyHaurwitz54Vorticity, &
!					scalarIntegralRefinement, circulationTol, "circulation refinement", &
!					RefineFlowMapYN = .TRUE., flowMapVarTol = flowMapVarTol )
				
			call Copy(sphere, tempSphere)
			
			sphere%tracers(2)%N = sphere%mesh%particles%N		
			do i = 1, sphere%mesh%particles%N
				sphere%tracers(2)%scalar(i) = Latitude( sphere%mesh%particles%x0(i), sphere%mesh%particles%y0(i), &
					sphere%mesh%particles%z0(i) )
			enddo
			
			call Delete(tempSphere)
			call Delete(remesh)
		elseif ( remeshCounter == resetLagParamInterval ) then
			!
			!	remesh to t = 0, keep old mesh as reference for next set of remeshings
			!
			call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logKey)//" remeshing to t = 0, remeshCounter = ", remeshCounter)
			refTime = sphere%mesh%t		
			call New(remesh, sphere)
			
			call New( refSphere, meshSeed, initNest, maxNest, amrLimit, radius, rotationRate )
			
			call AddTracers(refSphere, 2, [1, 1])
			refsphere%tracers(1)%name = "initLat"
			refsphere%tracers(1)%units = "radians"
			refsphere%tracers(2)%name = "lat(x0)"
			refsphere%tracers(2)%units = "radians"
			
			call LagrangianRemeshBVEWithVorticityFunction( remesh, sphere, refSphere, AMR, RossbyHaurwitz54Vorticity, &
				ScalarIntegralRefinement, circulationTol, "circulation refinement", &
				RefineFlowMapYN = .TRUE., flowMapVarTol = flowMapVarTol, nLagTracers = 1, tracerFn1 = InitLatTracer)
			
!			call LagrangianRemeshBVEWithVorticityFunction(remesh, sphere, refSphere, AMR, RossbyHaurwitz54Vorticity, &
!					scalarIntegralRefinement, circulationTol, "circulation refinement", &
!					RefineFlowMapYN = .TRUE., flowMapVarTol = flowMapVarTol )
			
			call Copy(sphere, refSphere)
			
			call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logKey)//" setting new ref_t = ", refTime)
			refSphere%mesh%t = refTime
			call ResetLagrangianParameter( sphere%mesh%particles )
			call ResetLagrangianParameter( refSphere%mesh%particles )
			
!			call Delete(remesh)
!			call New(remesh, sphere )
			
			sphere%tracers(2)%N = sphere%mesh%particles%N
			refSphere%tracers(2)%N = refSphere%mesh%particles%N		
			do i = 1, sphere%mesh%particles%N
				sphere%tracers(2)%scalar(i) = Latitude( sphere%mesh%particles%x0(i), sphere%mesh%particles%y0(i), &
					sphere%mesh%particles%z0(i) )
			enddo
		elseif ( remeshCounter > resetLagParamInterval .AND. mod(remeshCounter, resetLagParamInterval) == 0 ) then
			!
			!	remesh to previous reference mesh, then keep current mesh as new reference
			!
			refTime = sphere%mesh%t
			call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logKey)//" remeshing to t = ", refSphere%mesh%t)
			
			call New( tempSphere, meshSeed, initNest, maxNest, amrLimit, radius, rotationRate )
			call AddTracers(tempSphere, 2, [1, 1])
			tempsphere%tracers(1)%name = "initLat"
			tempsphere%tracers(1)%units = "radians"
			tempsphere%tracers(2)%name = "lat(x0)"
			tempsphere%tracers(2)%units = "radians"
			
			call LagrangianRemeshBVEToReferenceMesh( tempSphere, sphere, refSphere, remesh, &
					AMR, ScalarIntegralRefinement, circulationTol, "circulation refinement", &
					RefineFlowMapYN = .TRUE., flowMapVarTol = flowMapVarTol )
			call LogMessage(exeLog,DEBUG_LOGGING_LEVEL, trim(logKey)//" returned from LagRemesh ", "DEBUG2")		
			
			call Copy(sphere, tempSphere)
			call Copy(refSphere, tempSphere)
			
			call LogMessage(exeLog,DEBUG_LOGGING_LEVEL, trim(logKey)//" sphere copying complete ", "DEBUG4")		
			
			call Delete(remesh)
			call LogMessage(exeLog,DEBUG_LOGGING_LEVEL, trim(logKey)//" remesh delete complete ", "DEBUG5")
			call Delete(tempSphere)
			call LogMessage(exeLog,DEBUG_LOGGING_LEVEL, trim(logKey)//" tempSphere delete complete ", "DEBUG6")
			
			call ResetLagrangianParameter( sphere%mesh%particles )
			call ResetLagrangianParameter( refSphere%mesh%particles )
			
			call LogMessage(exeLog,DEBUG_LOGGING_LEVEL, trim(logKey)//" resetLagParam complete ", "DEBUG7")
			
			call New(remesh, refSphere)
			
			call LogMessage(exeLog,DEBUG_LOGGING_LEVEL, trim(logKey)//" refRemeshSource build complete ", "DEBUG8")
			
			refSphere%mesh%t = refTime
			call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logKey)//" setting new ref_t = ", refTime)
			sphere%tracers(2)%N = sphere%mesh%particles%N
			refsphere%tracers(2)%N = refsphere%mesh%particles%N		
			do i = 1, sphere%mesh%particles%N
				sphere%tracers(2)%scalar(i) = Latitude( sphere%mesh%particles%x0(i), sphere%mesh%particles%y0(i), &
					sphere%mesh%particles%z0(i) )
			enddo
		else
			!
			!	remesh to reference mesh
			!
			call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logKey)//" remeshing to t = ", refSphere%mesh%t)
			
			call StartSection(exeLog, "REFERENCE SPHERE STATUS")
			call LogStats(refSphere,exeLog)
			call EndSection(exeLog)
			
			call New( tempSphere, meshSeed, initNest, maxNest, amrLimit, radius, rotationRate )
			call AddTracers(tempSphere, 2, [1, 1])
			tempsphere%tracers(1)%name = "initLat"
			tempsphere%tracers(1)%units = "radians"
			tempsphere%tracers(2)%name = "lat(x0)"
			tempsphere%tracers(2)%units = "radians"
						
			call LogMessage(exeLog,DEBUG_LOGGING_LEVEL, trim(logKey)//" calling LagRemesh ", "DEBUG3")
			call LagrangianRemeshBVEToReferenceMesh( tempSphere, sphere, refSphere, remesh, &
					AMR, ScalarIntegralRefinement, circulationTol, "circulation refinement", &
					RefineFlowMapYN = .TRUE., flowMapVarTol = flowMapVarTol )
			call LogMessage(exeLog,DEBUG_LOGGING_LEVEL, trim(logKey)//" returned from LagRemesh ", "DEBUG1")
			
			call Copy(sphere, tempSphere)
			call Delete(tempSphere)
			
			sphere%tracers(2)%N = sphere%mesh%particles%N		
			do i = 1, sphere%mesh%particles%N
				sphere%tracers(2)%scalar(i) = Latitude( sphere%mesh%particles%x0(i), sphere%mesh%particles%y0(i), &
					sphere%mesh%particles%z0(i) )
			enddo
		endif
		
		call LogStats(sphere, exeLog)
		
		call Delete(solver)
		call New(solver, sphere)
	endif
	

	!
	!	advance time step
	!
	call Timestep(solver, sphere, dt)
	sphere%mesh%t = real(timeJ+1, kreal) * dt
	
	enstrophy(timeJ+2) = TotalEnstrophy(sphere)
	kineticEnergy(timeJ+2) = TotalKE(sphere)
	
	if ( procRank == 0 .AND. mod(timeJ+1, frameOut) == 0 ) then
		write(vtkFile,'(A,I0.4,A)') trim(vtkRoot), frameCounter, '.vtk'
		call OutputToVTK(sphere, vtkFile)
		frameCounter = frameCounter + 1
		call LogMessage(exelog, TRACE_LOGGING_LEVEL, trim(logKey)//" t = ", sphere%mesh%t)
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
call Delete(remesh)
call Delete(refSphere)
call Delete(solver)
call Delete(sphere)

call MPI_FINALIZE(mpiErrCode)

contains

subroutine SetAbsoluteTolerances( aBVEMesh, circTol, lagVarTol )
	type(BVEMesh), intent(in) :: aBVEMesh
	real(kreal), intent(inout) :: circTol
	real(kreal), intent(inout) :: lagVarTol
	!
	circTol = circTol * MaxCirculationMagnitudePerFace(aBVEMesh)
	lagVarTol = lagVarTol * MaxLagrangianVariationPerFace(aBVEMesh%mesh)
end subroutine

pure function legendre54( z )
	real(kreal) :: legendre54
	real(kreal), intent(in) :: z
	legendre54 = z * ( z * z - 1.0_kreal ) * ( z * z - 1.0_kreal)
end function

pure function RossbyHaurwitz54Vorticity( x, y, z )
	real(kreal) :: RossbyHaurwitz54Vorticity
	real(kreal), intent(in) :: x, y, z
	!
	real(kreal) :: lon
	
	lon = Longitude( x, y, z)
	RossbyHaurwitz54Vorticity = 2.0_kreal * zonalWind * z / radius + 30.0_kreal * rhWaveAmplitude * &
		cos( 4.0_kreal * lon ) * legendre54( z / radius ) / radius
end function

pure function InitLatTracer( x0, y0, z0 )
	real(kreal) :: InitLatTracer
	real(kreal), intent(in) :: x0, y0, z0
	InitLatTracer = Latitude( x0, y0, z0 )
end function

subroutine ReadNamelistFile( rank )
	integer(kint), intent(in) :: rank
	!
	character(len=MAX_STRING_LENGTH) :: namelistFilename
	integer(kint), parameter :: initBcast_intSize = 7
	integer(kint), parameter :: initBcast_realSize = 7
	integer(kint), dimension(initBcast_intSize) :: bcastIntegers
	real(kreal), dimension(initBcast_realSize) :: bcastReals
	integer(kint) :: mpiErrCode, readStat
	
	if ( COMMAND_ARGUMENT_COUNT() /= 1 ) then
		call LogMessage(exeLog, ERROR_LOGGING_LEVEL, trim(logKey), " ERROR: expected namelist file as 1st argument.")
		stop
	endif
	
	if ( rank == 0 ) then
		call GET_COMMAND_ARGUMENT(1, namelistFilename)
		
		open(unit=READ_UNIT, file=namelistFilename, status='OLD', action='READ', iostat=readStat)
			if ( readStat /= 0 ) then
				call LogMessage(exeLog, ERROR_LOGGING_LEVEL, trim(logKey), " ERROR: cannot read namelist file.")
				stop
			endif
		
			read(READ_UNIT, nml=meshDefine)
			rewind(READ_UNIT)
			read(READ_UNIT, nml=rhWave)
			rewind(READ_UNIT)
			read(READ_UNIT, nml=timestepping)
			rewind(READ_UNIT)
			read(READ_UNIT, nml=fileIO)
		close(READ_UNIT)
		
		if ( faceKind == 3 ) then
			meshSeed = ICOS_TRI_SPHERE_SEED
		elseif ( faceKind == 4) then
			meshSeed = CUBED_SPHERE_SEED
		else
			call LogMessage(exeLog, WARNING_LOGGING_LEVEL, trim(logkey)//" ReadNamelistFile WARNING : ", &
				" invalid faceKind -- using triangles.")
			meshSeed = ICOS_TRI_SPHERE_SEED
		endif
		
		maxNest = initNest + amrLimit

		bcastIntegers(1) = meshSeed
		bcastIntegers(2) = initNest
		bcastIntegers(3) = maxNest
		bcastIntegers(4) = amrLimit
		bcastIntegers(5) = frameOut
		bcastIntegers(6) = remeshInterval
		bcastIntegers(7) = resetLagParamInterval
		
		bcastReals(1) = dt
		bcastReals(2) = tfinal
		bcastReals(3) = rotationRate
		bcastReals(4) = rhWaveAmplitude
		bcastReals(5) = circulationTol
		bcastReals(6) = flowMapVarTol
		bcastReals(7) = zonalWind
	endif
	
	call MPI_BCAST(bcastIntegers, initBCAST_intSize, MPI_INTEGER, 0, MPI_COMM_WORLD, mpiErrCode)
	if ( mpiErrCode /= 0 ) then
		call LogMessage(exeLog, ERROR_LOGGING_LEVEL, trim(logKey)//" bcastIntegers, MPI_BCAST ERROR : ", mpiErrCode)
	endif
	
	call MPI_BCAST(bcastReals, initBCAST_realSize, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpiErrCode)
	if ( mpiErrCode /= 0 ) then
		call LogMessage(exeLog, ERROR_LOGGING_LEVEL, trim(logKey)//" bcastReals, MPI_BCAST ERROR : ", mpiErrCode)
	endif
	
	meshSeed = bcastIntegers(1)
	initNest = bcastIntegers(2)
	maxNest = bcastIntegers(3)
	amrLimit = bcastIntegers(4)
	frameOut = bcastIntegers(5)
	remeshInterval = bcastIntegers(6)
	resetLagParamInterval = bcastIntegers(7)
	
	dt = bcastReals(1)
	tfinal = bcastReals(2)
	rotationRate = bcastReals(3)
	rhWaveAmplitude = bcastReals(4)
	circulationTol = bcastReals(5)
	flowMapVarTol = bcastReals(6)
	zonalWind = bcastReals(7)
end subroutine


subroutine InitLogger(log, rank)
	type(Logger), intent(inout) :: log
	integer(kint), intent(in) :: rank
	if ( rank == 0 ) then
		call New(log, logLevel)
	else
		call New(log, ERROR_LOGGING_LEVEL)
	endif
	write(logKey,'(A,I0.2,A)') trim(logKey)//"_", rank, ":"
end subroutine


end program

