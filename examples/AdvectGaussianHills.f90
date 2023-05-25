program AdvectGaussianHillsDriver

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
use SphereTransportVelocitiesModule, velFn => LauritzenEtalDeformationalVelocity

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
real(kreal), allocatable, dimension(:) :: ghMass
real(kreal) :: l2Err, l2Denom
real(kreal) :: lInfErr
real(kreal) :: qMinTrue, qMaxTrue
real(kreal) :: qMinComp, qMaxComp
real(kreal) :: qMinErr, qMaxErr
real(kreal), parameter :: qRange = 0.95_kreal
integer(kint), parameter :: nTracers = 3
integer(kint), dimension(3), parameter :: tracerDims = [1,1,1]

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
logical(klog) :: useDirectRemesh
namelist /timestepping/ dt, tfinal, remeshInterval, useDirectRemesh

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
character(len=28) :: logKey = "AdvectGaussianHills"
integer(kint) :: mpiErrCode
real(kreal) :: timeStart, timeEnd
real(kreal) :: rmTimeStart, rmTimeEnd, rmTimeTotal
real(kreal) :: stepTimeStart, stepTimeEnd, stepTimeTotal
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
sphere%tracers(1)%name = "gaussianHills"
sphere%tracers(2)%name = "initialLatitude"
sphere%tracers(3)%name = "relError"
call SetInitialDensityOnMesh(sphere)
call SetTracerOnMesh( sphere, 1, GaussianHillsTracer )
sphere%tracers(3)%N = sphere%mesh%particles%N
call SetFieldToZero(sphere%tracers(3))

call SetVelocityOnMesh( sphere, velFn, t)
call SetDivergenceOnMesh(sphere)

do i = 1, sphere%mesh%particles%N
	vec = LagCoord(sphere%mesh%particles, i)
	call InsertScalarToField( sphere%tracers(2), Latitude(vec) )
enddo

! TO DO : AMR

!
!	Output initial data
!
frameCounter = 0
remeshCounter = 1
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
allocate(ghMass( nTimesteps + 1))
ghMass(1) = TracerMass(sphere, 1)

qMinTrue = MinScalarVal(sphere%tracers(1))
qMaxTrue = MaxScalarVal(sphere%tracers(1))

!--------------------------------
!	run : evolve the problem in time
!--------------------------------

call LogMessage(exeLog, DEBUG_LOGGING_LEVEL, trim(logkey)//" ", "starting timestepping loop.")

rmTimeTotal = 0.0_kreal
stepTimeTotal = 0.0_kreal
do timeJ = 0, nTimesteps - 1

	if ( mod(timeJ+1, remeshInterval) == 0 ) then
		rmTimeStart = MPI_WTIME()
		call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logkey)//" remesh triggered by remesh interval : ",&
			 remeshCounter)
		call New(remesh, sphere)

		call New(tempSphere, meshSeed, initNest, maxNest, amrLimit, radius, .FALSE.)
		call AddTracers(tempSphere, nTracers, tracerDims)
		tempSphere%tracers(1)%name = "gaussianHills"
		tempSphere%tracers(2)%name = "initialLatitude"
		tempSphere%tracers(3)%name = "relError"

		if ( useDirectRemesh ) then
			call DirectRemeshTransport( remesh, sphere, tempSphere, .FALSE., velFn, t )
		else
			call LagrangianRemeshTransportWithFunctions(remesh, sphere, tempSphere, .FALSE., velFn, t, &
				tracerFn1 = GaussianHillsTracer, tracerFn2 = InitLatTracer )
		endif

		call Copy(sphere, tempSphere)
		remeshCounter = remeshCounter + 1

		call Delete(tempSphere)
		call Delete(remesh)

		call Delete(solver)
		call New(solver, sphere)
		rmTimeEnd = MPI_WTIME()
		rmTimeTotal = rmTimeTotal + (rmTimeEnd - rmTimeStart)
	endif

	stepTimeStart = MPI_WTIME()
	call Timestep(solver, sphere, t, dt, velFn)
	stepTimeEnd = MPI_WTIME()
	stepTimeTotal = stepTimeTotal + (stepTimeEnd - stepTimeStart)

	t = real(timeJ +1, kreal) * dt
	sphere%mesh%t = t

	ghMass(timeJ+2) = TracerMass(sphere, 1)

	if ( timeJ+1 == nTimesteps ) then
	!--------------------------------
	!	calculate error at each particle
	!--------------------------------
		do i = 1, sphere%mesh%particles%N
			sphere%tracers(3)%scalar(i) = abs(GaussianHillsTracer( sphere%mesh%particles%x(i), &
				 sphere%mesh%particles%y(i), sphere%mesh%particles%z(i) ) - sphere%tracers(1)%scalar(i) ) / qRange
		enddo
	endif

	if ( procRank == 0 .AND. mod(timeJ+1, frameOut) == 0 ) then
		write(vtkFile,'(A,I0.4,A)') trim(vtkRoot), frameCounter, '.vtk'
		call OutputToVTK(sphere, vtkFile)
		frameCounter = frameCounter + 1
		call LogMessage(exelog, TRACE_LOGGING_LEVEL, trim(logKey)//" t = ", t)
	endif
enddo



!
!	write t = tfinal output
!
if ( procRank == 0 ) then
	write(matlabFile, '(5A)') trim(outputDir), '/', trim(outputRoot), trim(meshString), '.m'

	l2Err = 0.0_kreal
	l2Denom = 0.0_kreal
	do i = 1, sphere%mesh%particles%N
		if ( sphere%mesh%particles%isActive(i) ) then
			l2Err = l2Err + sphere%tracers(3)%scalar(i) * sphere%tracers(3)%scalar(i) * sphere%mesh%particles%area(i)
			l2Denom = l2Denom + GaussianHillsTracer( sphere%mesh%particles%x(i), sphere%mesh%particles%y(i), &
				  sphere%mesh%particles%z(i) ) ** 2 * sphere%mesh%particles%area(i)
		endif
	enddo
	l2Err = sqrt(l2Err / l2Denom)
	lInfErr = MaxScalarVal(sphere%tracers(3))

	qMinComp = MinScalarVal(sphere%tracers(1))
	qMaxComp = MaxScalarVal(sphere%tracers(1))

	qMinErr = (qMinComp - qMinTrue) / qRange
	qMaxErr = (qMaxComp - qMaxTrue) / qRange

	call StartSection(exeLog, "Final Errors: "//meshString )
		call LogMessage(exeLog, TRACE_LOGGING_LEVEL, "l2Err = ", l2Err )
		call LogMessage(exeLog, TRACE_LOGGING_LEVEL, "lInfErr = ", lInfErr )
		call LogMessage(exeLog, TRACE_LOGGING_LEVEL, "qMinErr = ", qMinErr )
		call LogMessage(exeLog, TRACE_LOGGING_LEVEL, "qMaxErr = ", qMaxErr )
		call LogMessage(exeLog, TRACE_LOGGING_LEVEL, "rel. tracer mass change = ", &
											maxval(abs(ghMass - ghMass(1))) / ghMass(1) )
		call LogMessage(exeLog, TRACE_LOGGING_LEVEL, " ", " ")
	call EndSection(exeLog)

	open(unit=WRITE_UNIT_1, file=matlabFile, status='REPLACE', action='WRITE')
		write(WRITE_UNIT_1,'(A,F12.9,A,F12.6,A)') "t = 0:", dt,":", tfinal, ";"
		call WriteToMatlab(l2Err, WRITE_UNIT_1, "l2Err")
		call WriteToMatlab(lInfErr, WRITE_UNIT_1, "lInfErr")
		call WriteToMatlab(ghMass, WRITE_UNIT_1, "ghMass")
	close(WRITE_UNIT_1)
endif

!--------------------------------
!	finalize : clean up
!--------------------------------

write(logString, '(A,2(F12.6,A))') "TOTAL REMESH TIME = ", rmTimeTotal, " seconds. Time per remesh = ", &
	rmTimeTotal / real(remeshCounter, kreal), " seconds."
call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logKey)//" ", logstring)

write(logString, '(A,2(F12.6,A))') "TOTAL RK4 TIME = ", stepTimeTotal, " seconds. Time per time step = ", &
	stepTimeTotal / real(nTimesteps, kreal), " seconds."
call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logKey)//" ", logstring)

timeEnd = MPI_WTIME()
write(logstring,'(A,F12.2,A)') "PROGRAM COMPLETE : elapsed time ", timeEnd - timeStart, " seconds."
call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logKey)//" ", logstring)

deallocate(ghMass)

call MPI_FINALIZE(mpiErrCode)

contains

subroutine ReadNamelistFile( rank )
	integer(kint), intent(in) :: rank
	!
	character(len=MAX_STRING_LENGTH) :: namelistFilename
	integer(kint), parameter :: initBcast_intSize = 7
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
		if (useDirectRemesh) then
			bcastIntegers(7) = 1
		else
			bcastIntegers(7) = 0
		endif

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
	if ( bcastIntegers(7) > 0 ) then
		useDirectRemesh = .TRUE.
	else
		useDirectRemesh = .FALSE.
	endif

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
