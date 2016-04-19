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
real(kreal), allocatable, dimension(:) :: l2Err
real(kreal), allocatable, dimension(:) :: lInfErr

! TO DO : remeshing for transport
integer(kint) :: remeshInterval
integer(kint) :: remeshCounter

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
character(len=28) :: logKey = "AdvectGaussianHills"
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
call AddTracers(sphere, 2, [1,1])
sphere%tracers(1)%name = "gaussianHills"
sphere%tracers(2)%name = "initialLatitude"
call SetInitialDensityOnMesh(sphere)
call SetTracerOnMesh( sphere, 1, GaussianHillsTracer )

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
allocate(ghMass( nTimesteps + 1))
allocate(lInfErr( nTimesteps + 1))
allocate(l2Err( nTimesteps + 1))
ghMass(1) = TracerMass(sphere, 1)


!--------------------------------
!	run : evolve the problem in time 
!--------------------------------

call LogMessage(exeLog, DEBUG_LOGGING_LEVEL, trim(logkey)//" ", "starting timestepping loop.")

do timeJ = 0, nTimesteps - 1
	
	call Timestep(solver, sphere, t, dt, velFn)
	
	t = real(timeJ +1, kreal) * dt
	sphere%mesh%t = t
	
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
	open(unit=WRITE_UNIT_1, file=matlabFile, status='REPLACE', action='WRITE')
		write(WRITE_UNIT_1,'(A,F12.9,A,F12.6,A)') "t = 0:", dt,":", tfinal, ";"
		call WriteToMatlab(ghMass, WRITE_UNIT_1, "ghMass")
	close(WRITE_UNIT_1)
endif

!--------------------------------
!	finalize : clean up
!--------------------------------

timeEnd = MPI_WTIME()
write(logstring,'(A,F12.2,A)') "PROGRAM COMPLETE : elapsed time ", timeEnd - timeStart, " seconds."
call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logKey)//" ", logstring)

deallocate(ghMass)
deallocate(lInfErr)
deallocate(l2Err)


call MPI_FINALIZE(mpiErrCode)

contains 

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
