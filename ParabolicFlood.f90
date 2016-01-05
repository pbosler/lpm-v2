program ParabolicFloodWave

use NumberKindsModule
use LoggerModule 
use PlanarSWEModule
use ParticlesModule
use FieldModule
use MPISetupModule
use PSEDirectSumModule
use SWEPlaneSolverModule

implicit none

include 'mpif.h'

!
!	computing environment
!
type(Logger) :: exeLog
character(len=MAX_STRING_LENGTH) :: logstring
character(len=15) :: logKey = "ParFlood"
integer(kint), parameter :: logLevel = DEBUG_LOGGING_LEVEL
character(len=MAX_STRING_LENGTH) :: namelistFilename
integer(kint) :: mpiErrCode
real(kreal) :: programStart, programEnd
integer(kint) :: i

!
!	mesh variables
!
type(SWEMesh) :: mesh
integer(kint) :: meshSeed
integer(kint) :: faceKind
integer(kint) :: initNest
integer(kint) :: maxNest
integer(kint) :: amrLimit
real(kreal), parameter :: meshRadius = 6.0_kreal
real(kreal), parameter :: f0 = 0.0_kreal
real(kreal), parameter :: beta = 0.0_kreal
real(kreal), parameter :: g = 1.0_kreal
character(len=MAX_STRING_LENGTH) :: meshString

namelist /meshDefine/ faceKind, initNest, maxNest, amrLimit

!
!	problem variables
!
real(kreal) :: r0, h0, TT
real(kreal) :: vel(2)
namelist /parabolicFlood/ r0, h0

!
!	timestepping
!
type(SWESolver) :: timekeeper
real(kreal) :: dt
real(kreal) :: tfinal
real(kreal) :: t
integer(kint) :: timeJ, nTimesteps
type(PSE) :: pseSetup

namelist /timestepping/ dt, tfinal

!
!	I/O variables
!
character(len=MAX_STRING_LENGTH) :: outputDir, outputRoot, vtkFile, matlabFile, vtkRoot
integer(kint) :: frameOut, frameCounter
namelist /fileIO/ outputDir, outputRoot, frameOut

!--------------------------------
!	initialize : setup computing environment
!--------------------------------
call MPI_INIT(mpiErrCode)
call MPI_COMM_SIZE(MPI_COMM_WORLD, numProcs, mpiErrCode)
call MPI_COMM_RANK(MPI_COMM_WORLD, procRank, mpiErrCode)
programStart = MPI_WTIME()

call InitLogger(exeLog, procRank)

call ReadNamelistFile(procRank)

!--------------------------------
!	initialize : setup problem / build mesh
!--------------------------------

call New(mesh, meshSeed, initNest, maxNest, amrLimit, meshRadius, f0, beta, g)
call SetInitialVorticityOnMesh( mesh, initVorticity )
call SetInitialDivergenceOnMesh( mesh, initDivergence)
call SetInitialVelocityOnMesh(mesh, initVelocity )
call SetInitialHOnMesh( mesh, initH )
call SetBottomHeightOnMesh( mesh, bottomTopo)
call SetInitialPotVortOnMesh( mesh )

call New(timekeeper, mesh, bottomTopo)
nTimesteps = floor(tfinal/dt)
t = 0.0_kreal
frameCounter = 1
call New(pseSetup, mesh%mesh)

call AddTracers(mesh, 3, [2,1,1])
mesh%tracers(1)%name = "exactVelocity"
mesh%tracers(1)%units = "m/s"
mesh%tracers(2)%name = "exactH"
mesh%tracers(2)%units = "m"
mesh%tracers(3)%name = "exactDivergence"
mesh%tracers(3)%units = "1/s"

TT = r0 / sqrt(2.0_kreal * g * h0 )

do i = 1, mesh%mesh%particles%N
	call InsertVectorToField(mesh%tracers(1), exactVelocity(mesh%mesh%particles%x(i), mesh%mesh%particles%y(i),t, TT ) )
	call InsertScalarToField(mesh%tracers(2), exactH(mesh%mesh%particles%x(i), mesh%mesh%particles%y(i),t, TT ) )
	call InsertScalarToField(mesh%tracers(3), exactDivergence(mesh%mesh%particles%x(i), mesh%mesh%particles%y(i),t, TT ) )
enddo

if ( procRank == 0 ) then
	call LogStats(mesh, exeLog)
	
	if ( meshSeed == TRI_HEX_SEED ) then
		write(meshString, '(A,I1,A)') '_triHex', initNest, '_'
	elseif ( meshSeed == QUAD_RECT_SEED ) then
		write(meshString, '(A,I1,A)') '_quadRect', initNest, '_'
	endif

	write(vtkRoot,'(4A)') trim(outputDir), '/vtkOut/', trim(outputRoot), trim(meshString)
	write(vtkFile,'(A,I0.4,A)') trim(vtkRoot), 0, '.vtk'
	
	call OutputToVTKFile(mesh, vtkFile)
	call LogMessage(exelog, TRACE_LOGGING_LEVEL, trim(logkey)//" t = ", t)
endif

!--------------------------------
!	run : evolve the problem in time (currently no remeshing)
!--------------------------------
do timeJ = 0, nTimesteps - 1
	
	call Timestep( timekeeper, mesh, dt, bottomTopo)
	
	t = real(timeJ+1,kreal) * dt
	mesh%mesh%t = t
	
	do i = 1, mesh%mesh%particles%N
		vel = exactVelocity(mesh%mesh%particles%x(i), mesh%mesh%particles%y(i), t, TT)
		mesh%tracers(1)%xComp(i) = vel(1)
		mesh%tracers(1)%yComp(i) = vel(2)
		mesh%tracers(2)%scalar(i) = exactH(mesh%mesh%particles%x(i), mesh%mesh%particles%y(i), t, TT)
		mesh%tracers(3)%scalar(i) = exactDivergence(mesh%mesh%particles%x(i), mesh%mesh%particles%y(i), t, TT)
	enddo
	
	if ( procRank == 0 .AND. mod(timeJ+1, frameOut) == 0 ) then
		write(vtkFile,'(A,I0.4,A)') trim(vtkRoot), frameCounter, '.vtk'
		
		call OutputToVTKFile( mesh, vtkFile )
			
		frameCounter = frameCounter + 1
		call LogMessage(exelog, TRACE_LOGGING_LEVEL, trim(logkey)//" t = ", t)
	endif	
enddo

!--------------------------------
!	finalize : clean up
!--------------------------------

programEnd = MPI_WTIME()
write(logString,'(A,F4.1,A)') "at t = ", tfinal, ","
call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logString)//" max relVort err = ", maxval(abs(mesh%relVort%scalar)))
call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logString)//" max h err = ", &
	maxval(abs(mesh%h%scalar - mesh%tracers(2)%scalar))/maxval(abs(mesh%tracers(2)%scalar)))
call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logString)//" max div err = ",&
	 maxval(abs(mesh%divergence%scalar - mesh%tracers(3)%scalar))/maxval(abs(mesh%tracers(3)%scalar)))

write(logstring,'(A,F12.2,A)') "PROGRAM COMPLETE : elapsed time ", programEnd - programStart, " seconds."
call LogMessage(exelog, TRACE_LOGGING_LEVEL, trim(logkey)//" ", logString)

call Delete(pseSetup)
call Delete(timekeeper)
call Delete(mesh)

call MPI_FINALIZE(mpiErrCode)

contains

function initDivergence( x, y )
	real(kreal) :: initDivergence
	real(kreal), intent(in) :: x, y
	initDivergence = 0.0_kreal
end function 

function initVorticity( x, y )
	real(kreal) :: initVorticity
	real(kreal), intent(in) :: x, y
	initVorticity = 0.0_kreal
end function 

function initVelocity(x, y)
	real(kreal) :: initVelocity(2)
	real(kreal), intent(in) :: x, y
	initVelocity = 0.0_kreal
end function

function initH(x, y)
	real(kreal) :: initH
	real(kreal), intent(in) :: x, y
	initH = max(h0 * ( 1.0_kreal - (x * x + y * y)/r0/r0), 0.0_kreal)
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

function bottomTopo(x,y)
	real(kreal) :: bottomTopo
	real(kreal), intent(in) :: x, y
	bottomTopo = 0.0_kreal
end function

function DepthLaplacian( x, y )
	real(kreal) :: DepthLaplacian
	real(kreal), intent(in) :: x, y
	DepthLaplacian = 0.0_kreal
end function

function exactVelocity( x, y, t, TT)
	real(kreal) :: exactVelocity(2)
	real(kreal), intent(in) :: x, y, t, TT
	exactVelocity(1) = x * t / ( t*t + TT*TT )
	exactVelocity(2) = y * t / ( t*t + TT*TT )
end function

function exactDivergence(x, y, t, TT)
	real(kreal) :: exactDivergence
	real(kreal), intent(in) :: x, y, t, TT
	exactDivergence = 2.0_kreal * t / ( t*t + TT*TT)
end function

function exactH(x, y, t, TT)
	real(kreal) :: exactH
	real(kreal), intent(in) :: x, y, t, TT
	exactH = max( h0 * ( TT*TT / ( t*t + TT*TT ) - (x*x + y*y) * TT**4 / ((t*t + TT*TT)**2)/r0/r0 ), 0.0_kreal)
end function

subroutine ReadNamelistFile(rank)
	integer(kint), intent(in) :: rank
	!
	integer(kint), parameter :: initBCAST_intSize = 5
	integer(kint), parameter :: initBCAST_realSize = 4
	integer(kint) :: readStat
	integer(kint) :: bcastIntegers(initBCAST_intSize)
	real(kreal) :: bcastReals(initBCAST_realSize)
	integer(kint) :: mpiErrCode
	
	if ( COMMAND_ARGUMENT_COUNT() /= 1 ) then
		call LogMessage(exeLog,ERROR_LOGGING_LEVEL, trim(logKey), " expected namelist file name as 1st argument.")
		stop
	endif
	
	if ( rank == 0 ) then
		call GET_COMMAND_ARGUMENT(1, namelistFilename)
		
		open( unit=READ_UNIT, file=namelistFilename, status='OLD', action='READ', iostat=readStat)
			if (readStat /= 0 ) then	
				call LogMessage(exeLog, ERROR_LOGGING_LEVEL, trim(logKey), " cannot read namelist file.")
				stop
			endif
			read(READ_UNIT, nml=meshDefine)
			rewind(READ_UNIT)
			read(READ_UNIT, nml=parabolicFlood)
			rewind(READ_UNIT)
			read(READ_UNIT, nml=timestepping)
			rewind(READ_UNIT)
			read(READ_UNIT, nml=fileIO)		
			! namelist /meshDefine/ faceKind, initNest, maxNest, amrLimit	
			
			if ( faceKind == 3 ) then
				meshSeed = TRI_HEX_SEED
			elseif ( faceKind == 4 ) then
				meshSeed = QUAD_RECT_SEED
			else
				call LogMessage(exeLog, WARNING_LOGGING_LEVEL, trim(logKey)//" ReadNamelistFile WARNING :", &
					" invalid faceKind -- using triangles.")
				meshSeed = TRI_HEX_SEED
			endif
			
			if (amrLimit == 0) maxNest = initNest
			
			bcastIntegers(1) = meshSeed
			bcastIntegers(2) = initNest
			bcastIntegers(3) = maxNest
			bcastIntegers(4) = amrLimit
			bcastIntegers(5) = frameOut
			
			bcastReals(1) = dt
			bcastReals(2) = tfinal
			bcastReals(3) = r0
			bcastReals(4) = h0
		close(READ_UNIT)
	endif
	
	call MPI_BCAST(bcastIntegers, initBCAST_intSize, MPI_INTEGER, 0, MPI_COMM_WORLD, mpiErrCode)
	
	call MPI_BCAST(bcastReals, initBCAST_realSize, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpiErrCode)
	
	meshSeed = bcastIntegers(1)
	initNest = bcastIntegers(2)
	maxNest = bcastIntegers(3)
	amrLimit = bcastIntegers(4)
	frameOut = bcastIntegers(5)
	
	dt = bcastReals(1)
	tfinal = bcastReals(2)
	r0 = bcastReals(3)
	h0 = bcastReals(4)
end subroutine
end program