program betaPlaneGaussian

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
use BetaPlaneMeshModule
use BetaPlaneSolverModule
use RefinementModule
use BIVARRemeshModule

implicit none

include 'mpif.h'

! mesh variables
type(BetaPlaneMesh) :: betaPlane
integer(kint) :: initNest
integer(kint) :: maxNest
integer(kint) :: amrLimit
real(kreal) :: lat0
real(kreal) :: rotationRate

! AMR variables
type(RefineSetup) :: refinement
logical(klog) :: AMR
integer(kint) :: nParticlesBefore
integer(kint) :: nParticlesAfter

namelist /meshDefine/ initNest, amrLimit, maxNest, lat0, rotationRate

! test case variables
real(kreal), allocatable, dimension(:) :: kineticEnergy, enstrophy
real(kreal) :: shapeParam
real(kreal) :: vorticityMax
real(kreal), parameter :: vortInitX = 0.5_kreal
real(kreal) :: vortInitY
real(kreal) :: flowMapVarTol
real(kreal) :: circulationTol
namelist /gaussVort/ shapeParam, vorticityMax, vortInitY, flowMapVarTol, circulationTol

! remeshing
integer(kint) :: remeshInterval
type(BetaPlaneMesh) :: tempMesh

! timestepping
type(BetaPlaneSolver) :: solver
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
character(len=28) :: logKey = "bPlaneGaussVort"
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
call New( betaPlane, initNest, maxNest, amrLimit, lat0, rotationRate )
call SetInitialVorticityOnMesh( betaPlane, GaussianVortexVorticity)

AMR = (amrLimit > 0 )
if ( AMR ) then
	!call LogMessage(exeLog, DEBUG_LOGGING_LEVEL, trim(logKey)//" ", "starting adaptive refinement...")

	call StartSection(exeLog, "--- Initial refinement ---")
	
	call New(refinement, betaPlane%mesh%faces%N_Max)

	call SetAbsoluteTolerances( betaPlane, circulationTol, flowMapVarTol )
	call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logKey)//" circulation tol = ", circulationTol)
	call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logKey)//" flowMapVarTol = ", flowMapVarTol)
	
	do i = 1, amrLimit
		call IterateMeshRefinementOneVariableAndFlowMap( refinement, betaPlane%mesh, betaPlane%relVort, &
			ScalarIntegralRefinement, circulationTol, "circulation refinement", flowMapVarTol, &
			nParticlesBefore, nParticlesAfter)
		call SetInitialVorticityOnMesh( betaPlane, GaussianVortexVorticity)
	enddo
	call LoadBalance(betaPlane%mpiParticles, betaPlane%mesh%particles%N, numProcs)
	call EndSection(exeLog)
endif

call SetStreamFunctionsOnMesh(betaPlane)
call SetVelocityOnMesh(betaPlane)

!
!	output initial state
!
frameCounter = 0
if ( procRank == 0 ) then
	call LogStats(betaPlane, exeLog)
	
	if ( AMR ) then
		write(meshString,'(A,I1,A,I2,A)') "_betaPlaneAMR", initNest, "to", initNest + amrLimit, "_"
	else
		write(meshString,'(A,I1,A)') "_betaPlane", initNest, "_"
	endif
	
	write(vtkRoot,'(4A)') trim(outputDir), '/vtkOut/', trim(outputRoot), trim(meshString)
	write(vtkFile,'(A,I0.4,A)') trim(vtkRoot), frameCounter, '.vtk'
	write(vtkRoot,'(4A)') trim(outputDir), '/vtkOut/', trim(outputRoot), trim(meshString)
	write(vtkFile,'(A,I0.4,A)') trim(vtkRoot), frameCounter, '.vtk'
	
	call OutputToVTK(betaPlane, vtkFile)
	frameCounter = frameCounter + 1
	
	call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logkey)//" t = ", t)
endif

!
!	initialize timestepping
!
call New(solver, betaPlane)
nTimesteps = floor(tfinal / dt)
allocate(kineticEnergy(nTimesteps+1))
allocate(enstrophy(nTimesteps+1))
kineticEnergy(1) = TotalKE(betaPlane)
enstrophy(1) = TotalEnstrophy(betaPlane)

!--------------------------------
!	run : evolve the problem in time 
!--------------------------------

do timeJ = 0, nTimesteps - 1
	
	
	if ( mod(timeJ + 1, remeshInterval) == 0 ) then
		call LogMessage(exeLog, DEBUG_LOGGING_LEVEL, trim(logKey)//" remesh : ", "remesh built.")
		
		call New( tempMesh, initNest, maxNest, amrLimit, lat0, rotationRate )
		
		call LogMessage(exeLog, DEBUG_LOGGING_LEVEL, trim(logKey)//" remesh : ", "new mesh built.")
		
		call LagrangianRemeshBetaPlaneWithVorticityFunction( betaPlane, tempMesh, AMR, GaussianVortexVorticity, &
				ScalarIntegralRefinement, circulationTol, "circulation refinement", RefineFlowMapYN = .TRUE., &
				flowMapVarTol = flowMapVarTol )
		
		call Copy(betaPlane, tempMesh)
		
		call Delete(solver)
		call Delete(tempMesh)
		
		call New(solver, betaPlane)
	endif
	
	call Timestep(solver, betaPlane, dt)
	t = real(timeJ + 1, kreal) * dt
	betaPlane%mesh%t = t
	
	kineticEnergy(timeJ + 2) = TotalKE(betaPlane)
	enstrophy(timeJ + 2) = TotalEnstrophy(betaPlane)
	
	if ( procRank == 0 .AND. mod(timeJ+1, frameOut) == 0 ) then
		write(vtkFile,'(A,I0.4,A)') trim(vtkRoot), frameCounter, '.vtk'
		call OutputToVTK(betaPlane, vtkFile)
		frameCounter = frameCounter + 1
		call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logKey)//" t = ", t)
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
call Delete(solver)
call Delete(betaPlane)

call MPI_FINALIZE(mpiErrCode)

contains

subroutine SetAbsoluteTolerances( betaPlane, circTol, lagVarTol )
	type(BetaPlaneMesh), intent(in) :: betaPlane
	real(kreal), intent(inout) :: circTol
	real(kreal), intent(inout) :: lagVarTol
	!
	circTol = circTol * MaxCirculationMagnitudePerFace(betaPlane)
	lagVarTol = lagVarTol * MaxLagrangianVariationPerFace(betaPlane%mesh)
end subroutine

function GaussianVortexVorticity( x, y )
	real(kreal) :: GaussianVortexVorticity
	real(kreal), intent(in) :: x, y
	GaussianVortexVorticity = vorticityMax * exp( - shapeParam * shapeParam * ( (x - vortInitX)**2  + (y - vortInitY)**2 ))
end function

subroutine ReadNamelistFile( rank )
	integer(kint), intent(in) :: rank
	!
	character(len=MAX_STRING_LENGTH) :: namelistFilename
	integer(kint), parameter :: initBcast_intSize = 5
	integer(kint), parameter :: initBcast_realSize = 9
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
			read(READ_UNIT, nml=gaussVort)
			rewind(READ_UNIT)
			read(READ_UNIT, nml=timestepping)
			rewind(READ_UNIT)
			read(READ_UNIT, nml=fileIO)
		close(READ_UNIT)
		
		bcastIntegers(1) = initNest
		bcastIntegers(2) = maxNest
		bcastIntegers(3) = amrLimit
		bcastIntegers(4) = frameOut
		bcastIntegers(5) = remeshInterval
		
		bcastReals(1) = dt
		bcastReals(2) = tfinal
		bcastReals(3) = shapeParam
		bcastReals(4) = rotationRate
		bcastReals(5) = vorticityMax
		bcastReals(6) = lat0
		bcastReals(7) = vortInitY
		bcastReals(8) = circulationTol
		bcastReals(9) = flowMapVarTol
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
	
	dt = bcastReals(1)
	tfinal = bcastReals(2)
	shapeParam = bcastReals(3)
	rotationRate = bcastReals(4)
	vorticityMax = bcastReals(5)
	lat0 = bcastReals(6)
	vortInitY = bcastReals(7)
	circulationTol = bcastReals(8)
	flowMapVarTol = bcastReals(9)
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