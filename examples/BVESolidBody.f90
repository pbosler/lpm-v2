program BVESolidBody

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

implicit none

include 'mpif.h'

! mesh variables
type(BVEMesh) :: sphere
integer(kint) :: initNest
integer(kint) :: maxNest
integer(kint) :: faceKind
integer(kint) :: meshSeed
integer(kint) :: amrLimit
real(kreal), parameter :: radius = 1.0_kreal
namelist /meshDefine/ faceKind, initNest, maxNest, amrLimit

! test case variables
real(kreal) :: omg = 2.0_kreal * PI
real(kreal), dimension(3) :: exactPos, exactVel
real(kreal) :: exactVor
real(kreal), allocatable, dimension(:) :: kineticEnergy, enstrophy

! timestepping
type(BVESolver) :: solver
real(kreal) :: dt
real(kreal) :: t
real(kreal) :: tfinal
integer(kint) :: nTimesteps
integer(kint) :: timeJ
namelist /timestepping/ dt, tfinal

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

! computing environment
type(Logger) :: exeLog
character(len=MAX_STRING_LENGTH) :: logString
integer(kint), parameter :: logLevel = DEBUG_LOGGING_LEVEL
character(len=18) :: logKey = "BVESolidBody"
integer(kint) :: mpiErrCode
real(kreal) :: timeStart, timeEnd
integer(kint) :: i

!--------------------------------
!	initialize : setup computing environment
!--------------------------------

call MPI_INIT(mpiErrCode)
call MPI_COMM_SIZE(MPI_COMM_WORLD, numProcs, mpiErrCode)
call MPI_COMM_RANK(MPI_COMM_WORLD, procRank, mpiErrCode)

call InitLogger(exeLog, procRank)

timeStart = MPI_WTIME()

call ReadNamelistFile( procRank )

!--------------------------------
!	initialize : build mesh / problem setup
!--------------------------------

call New( sphere, meshSeed, initNest, maxNest, amrLimit, radius, 0.0_kreal )
call SetInitialVorticityOnMesh(sphere, SolidBodyVorticity)
call SetStreamFunctionsOnMesh(sphere)
call SetVelocityOnMesh(sphere)

call New(solver, sphere)
nTimesteps = floor( tfinal / dt )
t = 0.0_kreal

allocate(kineticEnergy(nTimesteps+1))
allocate(enstrophy(nTimesteps+1))
kineticEnergy(1) = TotalKE(sphere)
enstrophy(1) = TotalEnstrophy(sphere)

call AddTracers(sphere, 3, [3, 3, 1])

sphere%tracers(1)%N = sphere%mesh%particles%N
sphere%tracers(2)%N = sphere%mesh%particles%N
sphere%tracers(3)%N = sphere%mesh%particles%N
do i = 1, sphere%tracers(1)%N	
	exactPos = ExactPosition( sphere%mesh%particles%x0(i), sphere%mesh%particles%y0(i), sphere%mesh%particles%z0(i) )
	sphere%tracers(1)%xComp(i) = sphere%mesh%particles%x(i) - exactPos(1)
	sphere%tracers(1)%yComp(i) = sphere%mesh%particles%y(i) - exactPos(2)
	sphere%tracers(1)%zComp(i) = sphere%mesh%particles%z(i) - exactPos(3)

	exactVel = SolidBodyVelocity( sphere%mesh%particles%x(i), sphere%mesh%particles%y(i), sphere%mesh%particles%z(i))
	sphere%tracers(2)%xComp(i) = sphere%velocity%xComp(i) - exactVel(1)
	sphere%tracers(2)%yComp(i) = sphere%velocity%yComp(i) - exactVel(2)
	sphere%tracers(2)%zComp(i) = sphere%velocity%zComp(i) - exactVel(3)
	
	exactVor = SolidBodyVorticity( sphere%mesh%particles%x(i), sphere%mesh%particles%y(i), sphere%mesh%particles%z(i))
	sphere%tracers(3)%scalar(i) = sphere%relVort%scalar(i) - exactVor
enddo
sphere%tracers(1)%name = "positionError"
sphere%tracers(2)%name = "velocityError"
sphere%tracers(3)%name = "vorticityError"
call MultiplyFieldByScalar(sphere%tracers(2), 1.0_kreal / (2.0_kreal * PI ))
call MultiplyFieldByScalar(sphere%tracers(3), 1.0_kreal / (2.0_kreal * PI ))

frameCounter = 0
if ( procRank == 0 ) then
	call LogStats(sphere, exeLog)
	
	if ( meshSeed == ICOS_TRI_SPHERE_SEED ) then
		write(meshString,'(A,I1,A)') '_icosTri', initNest, '_'
	elseif (meshSeed == CUBED_SPHERE_SEED ) then
		write(meshString,'(A,I1,A)') '_cubedSph', initNest, '_'
	endif
	
	write(vtkRoot,'(4A)') trim(outputDir), '/vtkOut/', trim(outputRoot), trim(meshString)
	write(vtkFile,'(A,I0.4,A)') trim(vtkRoot), frameCounter, '.vtk'
	
	call OutputToVTK( sphere, vtkFile )
	frameCounter = frameCounter + 1
	
	call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logKey)//" t = ", t )	
endif

!--------------------------------
!	run : evolve the problem in time (currently no remeshing)
!--------------------------------
do timeJ = 0, nTimesteps - 1
	
	call Timestep(solver, sphere, dt )
	
	t = real(timeJ + 1, kreal) * dt
	sphere%mesh%t = t
	
	do i = 1, sphere%tracers(1)%N
		exactPos = ExactPosition( sphere%mesh%particles%x0(i), sphere%mesh%particles%y0(i), sphere%mesh%particles%z0(i) )
		sphere%tracers(1)%xComp(i) = sphere%mesh%particles%x(i) - exactPos(1)
		sphere%tracers(1)%yComp(i) = sphere%mesh%particles%y(i) - exactPos(2)
		sphere%tracers(1)%zComp(i) = sphere%mesh%particles%z(i) - exactPos(3)
		
		exactVel = SolidBodyVelocity( sphere%mesh%particles%x(i), sphere%mesh%particles%y(i), sphere%mesh%particles%z(i))
		sphere%tracers(2)%xComp(i) = sphere%velocity%xComp(i) - exactVel(1)
		sphere%tracers(2)%yComp(i) = sphere%velocity%yComp(i) - exactVel(2)
		sphere%tracers(2)%zComp(i) = sphere%velocity%zComp(i) - exactVel(3)
		
		exactVor = SolidBodyVorticity( sphere%mesh%particles%x(i), sphere%mesh%particles%y(i), sphere%mesh%particles%z(i))
		sphere%tracers(3)%scalar(i) = sphere%relVort%scalar(i) - exactVor
	enddo
	call MultiplyFieldByScalar(sphere%tracers(2), 1.0_kreal / (2.0_kreal * PI ))
call MultiplyFieldByScalar(sphere%tracers(3), 1.0_kreal / (2.0_kreal * PI ))
	
	enstrophy(timeJ + 2) = TotalEnstrophy(sphere)
	kineticEnergy(timeJ + 2) = TotalKE(sphere)
	
	if ( procRank == 0 .AND. mod(timeJ + 1, frameOut) == 0 ) then
		write(vtkFile,'(A,I0.4,A)') trim(vtkRoot), frameCounter, '.vtk'
		
		call OutputToVTK( sphere, vtkFile )
		frameCounter = frameCounter + 1
		
		call Logmessage(exeLog, TRACE_LOGGING_LEVEL, trim(logKey)//" t = ", t)
	endif
enddo

! write final output
if ( procRank == 0 ) then
	write(matlabFile,'(5A)') trim(outputDir), '/', trim(outputRoot), trim(meshString), '.m'
	open(unit = WRITE_UNIT_1, file=matlabFile, status='REPLACE', action='WRITE')
		call WriteToMatlab( kineticEnergy, WRITE_UNIT_1, "kineticEnergy")
		call WriteToMatlab( enstrophy, WRITE_UNIT_1, "enstrophy")
	close(WRITE_UNIT_1)
endif

!--------------------------------
!	finalize : clean up
!--------------------------------

call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logKey)//" final velocity rel. err = ", &
		MaxMagnitude(sphere%tracers(2)))
call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logKey)//" final position err = ", MaxMagnitude(sphere%tracers(1)) )
call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logKey)//" final vorticity err = ", &
		maxval(abs(sphere%tracers(3)%scalar)) )

timeEnd = MPI_WTIME()
write(logstring,'(A,F12.2,A)') "PROGRAM COMPLETE : elapsed time ", timeEnd - timeStart, " seconds."
call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logKey)//" ", logstring)

deallocate(kineticEnergy)
deallocate(enstrophy)
call Delete(solver)
call Delete(sphere)

call MPI_FINALIZE(mpiErrCode)

contains

pure function ExactPosition( x0, y0, z0 )
	real(kreal), dimension(3) :: ExactPosition
	real(kreal), intent(in) :: x0, y0, z0
	real(kreal), dimension(3,3) :: rotMat
	
	rotMat(1,1) = cos( omg * t )
	rotMat(2,1) = sin( omg * t )
	rotMat(3,1) = 0.0_kreal
	rotMat(1,2) = -sin( omg * t )
	rotMat(2,2) = cos( omg * t )
	rotMat(3,2) = 0.0_kreal
	rotMat(1,3) = 0.0_kreal
	rotMat(2,3) = 0.0_kreal
	rotMat(3,3) = 1.0_kreal
	
	ExactPosition = MATMUL( rotMat, [x0, y0, z0])
end function

pure function SolidBodyVorticity( x, y, z)
	real(kreal) :: SolidBodyVorticity
	real(kreal), intent(in) :: x, y, z
	SolidBodyVorticity = 2.0_kreal * omg * z / radius
end function

pure function SolidBodyVelocity( x, y, z)
	real(kreal), dimension(3) :: SolidBodyVelocity
	real(kreal), intent(in) :: x, y, z
	SolidBodyVelocity(1) = - omg * y
	SolidBodyVelocity(2) =   omg * x
	SolidBodyVelocity(3) = 0.0_kreal
end function

subroutine ReadNamelistFile( rank )
	integer(kint), intent(in) :: rank
	!
	character(len=MAX_STRING_LENGTH) :: namelistFilename
	integer(kint), parameter :: initBcast_intSize = 5
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
		
		if ( amrLimit == 0 ) maxNest = initNest

!		namelist /meshDefine/ faceKind, initNest, maxNest, amrLimit, radius
!		namelist /timestepping/ dt, tfinal		
!		namelist /fileIO/ outputDir, outputRoot, frameOut

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