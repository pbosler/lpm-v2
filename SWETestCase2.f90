program SWETestCase2

use NumberKindsModule
use OutputWriterModule
use LoggerModule
use ParticlesModule
use EdgesModule
use FacesModule
use PolyMesh2dModule
use FieldModule
use MPISetupModule
use PSEDirectSumModule
use SphereGeomModule
use SphereSWEModule
use SphereSWESolverModule

implicit none

include 'mpif.h'

! mesh variables
type(SWEMesh) :: sphere
integer(kint) :: initNest
integer(kint) :: maxNest
integer(kint) :: meshSeed
integer(kint) :: amrLimit
integer(kint) :: faceKind
namelist /meshDefine/ faceKind, initNest, maxNest, amrLimit

! test case variables
real(kreal), parameter :: radius = 1.0_kreal
real(kreal), parameter :: u0 = 2.0_kreal * PI / 12.0_kreal
real(kreal), parameter :: g = 10000.0_kreal
real(kreal), parameter :: h0 = 0.0005_kreal
real(kreal), parameter :: omg = 0.0_kreal

!
! timestepping
!
type(SWESolver) :: solver
real(kreal) :: dt
real(kreal) :: t
real(kreal) :: tfinal
integer(kint) :: timeJ
integer(kint) :: nTimesteps
type(PSE) :: pseSetup

namelist /timestepping/ dt, tfinal

!
! i/o variables
!
character(len=MAX_STRING_LENGTH) :: outputDir
character(len=MAX_STRING_LENGTH) :: outputRoot
character(len=MAX_STRING_LENGTH) :: vtkFile
character(len=MAX_STRING_LENGTH) :: matlabFile
character(len=MAX_STRING_LENGTH) :: vtkRoot
character(len=MAX_STRING_LENGTH) :: meshString
integer(kint) :: frameOut
integer(kint) :: frameCounter
namelist /fileIO/ outputDir, outputRoot, frameOut

!
!	computing environment
!
type(Logger) :: exeLog
character(len=MAX_STRING_LENGTH) :: logstring
character(len=15) :: logKey = "SWE_TC2"
integer(kint), parameter :: logLevel = DEBUG_LOGGING_LEVEL
character(len=MAX_STRING_LENGTH) :: namelistFilename
integer(kint) :: mpiErrCode
real(kreal) :: programStart, programEnd

!--------------------------------
!	initialize : setup computing environment
!--------------------------------
call MPI_INIT(mpiErrCode)
call MPI_COMM_SIZE(MPI_COMM_WORLD, numProcs, mpiErrCode)
call MPI_COMM_RANK(MPI_COMM_WORLD, procRank, mpiErrCode)

call InitLogger(exeLog, procRank)

programStart = MPI_WTIME()

call ReadNamelistFile(procRank)

!--------------------------------
!	initialize : build mesh / problem setup
!--------------------------------

call New(sphere, meshSeed, initNest, maxNest, amrLimit, radius, omg, g)
call SetInitialRelVortOnMesh( sphere, TestCase2Vorticity )
call SetInitialDivergenceOnMesh( sphere, TestCase2Divergence )
call SetInitialVelocityOnMesh( sphere, TestCase2Velocity )
call SetInitialHOnMesh( sphere, TestCase2Height )
call SetInitialPotVortOnMesh(sphere)
call SetBottomHeightOnMesh(sphere, topography)
call AddTracers( sphere, 4, [1,1,1,1])
sphere%tracers(1)%name = "doubleDot1"
sphere%tracers(2)%name = "lapEta2"
sphere%tracers(3)%name = "doubleDot3"
sphere%tracers(4)%name = "lapEta3"

sphere%tracers(1)%N = sphere%mesh%particles%N
sphere%tracers(2)%N = sphere%mesh%particles%N
sphere%tracers(3)%N = sphere%mesh%particles%N
sphere%tracers(4)%N = sphere%mesh%particles%N


call New(pseSetup, sphere%mesh)
sphere%pseEps = pseSetup%eps

call New(solver, sphere, topography)
nTimesteps = floor(tfinal/dt)
t = 0.0_kreal
frameCounter = 1

if ( procRank == 0 ) then

	call LogStats(sphere, exeLog)

	if ( meshSeed == ICOS_TRI_SPHERE_SEED ) then
		write(meshString,'(A,I1,A)') '_icosTri', initNest, '_'
	elseif (meshSeed == CUBED_SPHERE_SEED ) then
		write(meshString,'(A,I1,A)') '_cubedSph', initNest, '_'
	endif
	
	write(vtkRoot,'(4A)') trim(outputDir), '/vtkOut/', trim(outputRoot), trim(meshString)
	write(vtkFile,'(A,I0.4,A)') trim(vtkRoot), 0, '.vtk'
	
	call OutputToVTK( sphere, vtkFile )
	
	call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logKey)//" t = ", t )	
endif

!--------------------------------
!	run : evolve the problem in time (currently no remeshing)
!--------------------------------

do timeJ = 0, nTimesteps - 1

	call Timestep( solver, sphere, dt, topography )
	
	t = real(timeJ + 1, kreal) * dt
	sphere%mesh%t = t
	
	
	if ( procRank == 0 .AND. mod(timeJ+1, frameOut) == 0 ) then
		write(vtkFile,'(A,I0.4,A)') trim(vtkRoot), frameCounter, '.vtk'
		
		call OutputToVTK( sphere, vtkFile )
		frameCounter = frameCounter + 1
		
		call LogStats(sphere, exeLog)
		
		call Logmessage(exeLog, TRACE_LOGGING_LEVEL, trim(logKey)//" t = ", t)
	endif
enddo

!--------------------------------
!	finalize : clean up
!--------------------------------

programEnd = MPI_WTIME()
write(logstring,'(A,F12.2,A)') "PROGRAM COMPLETE : elapsed time ", programEnd - programStart, " seconds."
call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logKey)//" ", logstring)

call Delete(pseSetup)
call Delete(solver)
call Delete(sphere)

call MPI_FINALIZE(mpiErrCode)

contains

subroutine ReadNamelistFile( rank )
	integer(kint), intent(in) :: rank
	!
	character(len=MAX_STRING_LENGTH) :: namelistFilename
	integer(kint), parameter :: initBCAST_intSize = 5
	integer(kint), parameter :: initBCAST_realSize = 2
	integer(kint) :: readStat
	integer(kint), dimension(initBCAST_intSize) :: bcastIntegers
	real(kreal), dimension(initBCAST_realSize) :: bcastReals
	integer(kint) :: mpiErrCode
	
	if ( COMMAND_ARGUMENT_COUNT() /=1 ) then
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
		
		if ( amrLimit == 0 ) maxNest = initNest
		
		bcastIntegers(1) = meshSeed
		bcastIntegers(2) = initNest
		bcastIntegers(3) = maxNest
		bcastIntegers(4) = amrLimit
		bcastIntegers(5) = frameOut
		
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
	
	dt = bcastReals(1)
	tfinal = bcastReals(2)
end subroutine

pure function TestCase2Vorticity( x, y, z)
!   Equation 94 from Williamson et al. 1992
	real(kreal) :: TestCase2Vorticity
	real(kreal), intent(in) :: x, y, z
	!
	real(kreal) :: lat
	lat = Latitude(x, y, z)
	TestCase2Vorticity = 2.0_kreal * u0 / radius * sin(lat)
end function

function TestCase2Height( x, y, z )
!   Equation 95 from Williamson et al. 1992
	real(kreal) :: TestCase2Height
	real(kreal), intent(in) :: x, y, z
	!
	real(kreal) :: lat
	lat = Latitude(x, y, z)
	TestCase2Height = h0 - ( radius * omg * u0 + 0.5_kreal * u0*u0 ) * ( sin(lat) * sin(lat) ) / g
end function

pure function topography( x, y, z )
	real(kreal) :: topography
	real(kreal), intent(in) :: x, y, z
	topography = 0.0_kreal
end function 

pure function TestCase2Divergence( x, y, z )
	real(kreal) :: TestCase2Divergence
	real(kreal), intent(in) :: x, y, z
	TestCase2Divergence = 0.0_kreal
end function

pure function TestCase2Velocity( x, y, z )
	real(kreal), dimension(3) :: TestCase2Velocity
	real(kreal), intent(in) :: x, y, z
	
	TestCase2Velocity(1) = -u0 * y
	TestCase2Velocity(2) =  u0 * x
	TestCase2Velocity(3) = 0.0_kreal 	
end function

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
