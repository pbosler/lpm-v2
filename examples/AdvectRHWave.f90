program RHWaveTransport

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
use RefinementModule
use SSRFPACKInterfaceModule
use SSRFPACKRemeshModule
use SphereTransportSolverModule
use SphereTracersModule
use SphereTransportModule
use SphereTransportVelocitiesModule, velFn => RH4Velocity

implicit none

include 'mpif.h'

! mesh variables
type(TransportMesh) :: sphere
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
namelist /meshDefine/ faceKind, initNest, amrLimit

! test case variables
real(kreal), allocatable, dimension(:) :: trMass
integer(kint), parameter :: nTracers = 2
integer(kint), dimension(2), parameter :: tracerDims = [1,1]

! remeshing variables
type(TransportRemesh) :: remesh
integer(kint) :: remeshInterval
integer(kint) :: remeshCounter
type(TransportMesh) :: tempSphere

! timestepping
type(TransportSolver) :: solver
real(kreal) :: dt
real(kreal) :: t
real(kreal) :: tfinal
integer(kint) :: nTimesteps
integer(kint) :: timeJ
namelist /timestepping/ dt, tfinal, remeshInterval

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
character(len=28) :: logKey = "AdvectRHWave"
integer(kint) :: mpiErrCode
real(kreal) :: timeStart, timeEnd
integer(kint) :: i
real(kreal), dimension(3) :: vec

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
call New( sphere, meshSeed, initNest, maxNest, amrLimit, radius, .FALSE.)
call AddTracers(sphere, nTracers, tracerDims)
sphere%tracers(1)%name = "initialLatitude"
sphere%tracers(1)%units = "radians"
sphere%tracers(2)%name = "vorticity"
sphere%tracers(2)%units = "1/s"
call SetInitialDensityOnMesh(sphere)
do i = 1, sphere%mesh%particles%N
	vec = LagCoord(sphere%mesh%particles, i)
	call InsertScalarToField( sphere%tracers(1), Latitude(vec) )
enddo
call SetTracerOnMesh( sphere, 2, RH54Vorticity )
call SetVelocityOnMesh( sphere, velFn, t)
call SetDivergenceOnMesh(sphere)

call LogStats(sphere%mpiParticles, exeLog)

! TO DO : AMR, initial refinement

!
!	Output initial data
!
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
nTimesteps = floor( tfinal/dt )
allocate(trMass( nTimesteps + 1))
trMass(1) = TracerMass(sphere,1)

!--------------------------------
!	run : evolve the problem in time 
!--------------------------------
call LogMessage(exeLog, DEBUG_LOGGING_LEVEL, trim(logkey)//" ", "starting timestepping loop.")

do timeJ = 0, nTimesteps - 1

	if ( mod(timeJ+1, remeshInterval) == 0 ) then
		call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logkey)//" remesh triggered by remesh interval : ",& 
			 remeshCounter)
		call New(remesh, sphere)
		
		call New(tempSphere, meshSeed, initNest, maxNest, amrLimit, radius, .FALSE.)
		call SetDivergenceOnMesh(tempSphere)
		
		call AddTracers(tempSphere, nTracers, tracerDims)
		tempSphere%tracers(1)%name = "initialLatitude"
		tempSphere%tracers(1)%units = "radians"
		tempSphere%tracers(1)%units = "radians"
		tempSphere%tracers(2)%name = "vorticity"
		tempSphere%tracers(2)%units = "1/s"
		
		! TO DO : Adaptive remeshing
		
		call LagrangianRemeshTransportWithFunctions(remesh, sphere, tempSphere, .FALSE., velFn, t, &
			tracerFn1 = InitLatTracer, tracerFn2 = RH54Vorticity )
		
		call LogMessage(exeLog, DEBUG_LOGGING_LEVEL, "remesh loop: ", "Remesh::LagRemesh has returned")
		
		call Copy(sphere, tempSphere)	
		remeshCounter = remeshCounter + 1
		
		call Delete(tempSphere)
		call Delete(remesh)
		
		call Delete(solver)
		call New(solver, sphere)
	endif

	call Timestep(solver, sphere, t, dt, velFn)
	t = real(timeJ +1, kreal) * dt
	sphere%mesh%t = t
	
	trMass(timeJ+2) = TracerMass(sphere, 1)

	if ( procRank == 0 .AND. mod(timeJ+1, frameOut) == 0 ) then
		write(vtkFile,'(A,I0.4,A)') trim(vtkRoot), frameCounter, '.vtk'
		call OutputToVTK(sphere, vtkFile)
		frameCounter = frameCounter + 1
		call LogMessage(exelog, TRACE_LOGGING_LEVEL, trim(logKey)//" t = ", t)
	endif
enddo 

!--------------------------------
!	finalize : clean up
!--------------------------------

timeEnd = MPI_WTIME()
write(logstring,'(A,F12.2,A)') "PROGRAM COMPLETE : elapsed time ", timeEnd - timeStart, " seconds."
call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logKey)//" ", logstring)

deallocate(trMass)
call Delete(solver)

call MPI_FINALIZE(mpiErrCode)

contains 

!--------------------------------
!	subprograms used by this application
!--------------------------------

subroutine ReadNamelistFile( rank )
	integer(kint), intent(in) :: rank
	!
	character(len=MAX_STRING_LENGTH) :: namelistFilename
	integer(kint), parameter :: initBcast_intSize = 6
	integer(kint), parameter :: initBcast_realSize = 2
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

!		namelist /meshDefine/ faceKind, initNest, maxNest, amrLimit, radius
!		namelist /timestepping/ dt, tfinal		
!		namelist /fileIO/ outputDir, outputRoot, frameOut

		bcastIntegers(1) = meshSeed
		bcastIntegers(2) = initNest
		bcastIntegers(3) = maxNest
		bcastIntegers(4) = amrLimit
		bcastIntegers(5) = frameOut
		bcastIntegers(6) = remeshInterval
		
		bcastReals(1) = dt
		bcastReals(2) = tfinal
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
	
	dt = bcastReals(1)
	tfinal = bcastReals(2)
end subroutine

!> @brief Initializes a @ref Logger for this executable program.
!> 
!> Output is controlled by message priority level and MPI rank.
!> 
!> @param[in] log @ref Logger to initialize
!> @param[in] rank MPI rank

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
