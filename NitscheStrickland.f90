program NitcheStrickland
!> @file NitscheStrickland.f90
!> Driver program for a planar SWE problem defined by a symmetric vortex and zero initial divergence.
!> @author Peter Bosler, Sandia National Laboratories, Center for Computing Research
!>
!> Demonstrates flow interaction with topography.@n
!> Demontstrates the use of @ref PlanarSWE, @ref SWEPlaneSolver, and @ref PSEDirectSum. @n
!> Driver program based on the setup of an isentropic gas dymanics problem from Nitsche and Strickland (2002).@n
!> Intial conditions are a symmetric vortex located at the origin, uniform fluid depth and zero initial divergence.
!> 
!> Users may alter the source code to change the bottom topography function (each change requires a new compilation).@n
!> 
!> Other input is supplied via namelist file.  
!>
!> @image html NS2002WithTopo.png "Flow interaction with topography. A vortex interacting with a mountain generates gravity waves."
!>
!> For a complete description of the problem that this application is based upon, see @n
!> [M. Nitsche and J. Strickland, Extension of the gridless vortex method into the compressible flow regime. <i> J. Turb. </i> 3 (2002)](http://www.math.unm.edu/~nitsche/webhome/preprints.html)
!> 
use NumberKindsModule
use LoggerModule 
use PlanarSWEModule
use ParticlesModule
use FieldModule
use MPISetupModule
use SWEPlaneSolverModule
use PSEDirectSumModule
use BIVARRemeshModule

implicit none

include 'mpif.h'

!
!	computing environment
!
type(Logger) :: exeLog
character(len=MAX_STRING_LENGTH) :: logstring
character(len=15) :: logKey = "NS2002"
integer(kint), parameter :: logLevel = DEBUG_LOGGING_LEVEL
character(len=MAX_STRING_LENGTH) :: namelistFilename
integer(kint) :: mpiErrCode
real(kreal) :: programStart, programEnd
!
!	mesh variables
!
type(SWEMesh) :: mesh
integer(kint) :: meshSeed
integer(kint) :: faceKind
integer(kint) :: initNest
integer(kint) :: maxNest
integer(kint) :: amrLimit
real(kreal), parameter :: meshRadius = 10.0_kreal
real(kreal) :: f0 
real(kreal) :: beta 
real(kreal), parameter :: g = 1.0_kreal
character(len=MAX_STRING_LENGTH) :: meshString

namelist /meshDefine/ faceKind, initNest, maxNest, amrLimit, f0, beta
!
!	timestepping
!
type(SWESolver) :: timekeeper
real(kreal) :: dt
real(kreal) :: tfinal
real(kreal) :: t
integer(kint) :: timeJ, nTimesteps
type(PSE) :: pseSetup

!
!	remeshing
!
integer(kint) :: remeshInterval
integer(kint) :: remeshCounter
type(SWEMesh) :: tempMesh

namelist /timestepping/ dt, tfinal, remeshInterval
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

print *, "hello from proc ", procRank

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

call New(pseSetup, mesh%mesh)
mesh%pseEps = pseSetup%eps

call New(timekeeper, mesh, bottomTopo)
nTimesteps = floor(tfinal/dt)
t = 0.0_kreal
frameCounter = 1


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
remeshCounter = 0
!print *, "proc ", procRank, ": starting time stepping loop."
do timeJ = 0, nTimesteps - 1
	!
	!	remesh
	!
	if (mod(timeJ+1, remeshInterval) == 0) then
		call MPI_BARRIER(MPI_COMM_WORLD, mpiErrCode)
!		print *, "proc ", procRank, ": entering remesh."
		remeshCounter = remeshCounter + 1
		call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logkey)//" remesh ", remeshCounter)
		call New(tempMesh, meshSeed, initNest, maxNest, amrLimit, meshRadius, f0, beta, g)
		call SetFieldN(tempMesh)
		
		if (procRank == 0 ) then
			call OutputToVTKFile(tempMesh, trim(outputDir)//"/vtkOut/preInterp.vtk")
		endif
		
		print *, "proc ", procRank, ": new mesh built."
		
		call SetBottomHeightOnMesh(tempMesh, bottomTopo)
		call DirectRemeshPlanarSWE(mesh, tempMesh, .FALSE.)
		
		if (procRank == 0 ) then
			call OutputToVTKFile(tempMesh, trim(outputDir)//"/vtkOut/postInterp.vtk")
		endif
		
		call Copy(mesh, tempMesh)
		
!		call StartSection(exeLog, "post copy")
!		call LogStats(mesh, exeLog)
!		call EndSection(exeLog)
		
		call Delete(tempMesh)
		
		call Delete(timekeeper)

		call Delete(pseSetup)		
		call New(pseSetup, mesh%mesh)
		mesh%pseEps = pseSetup%eps
		
		call New(timekeeper, mesh, bottomTopo)
!		call LogStats(timekeeper, exeLog)
	endif
	
	call Timestep( timekeeper, mesh, dt, bottomTopo)
!	if (timeJ + 1 == remeshInterval) then
!		call LogStats(timekeeper, exeLog)
!		call StartSection(exeLog, "post timestep")
!		call LogStats(mesh, exeLog)
!		call EndSection(exeLog)
!	endif
	
	t = real(timeJ+1,kreal) * dt
	mesh%mesh%t = t
	
	if ( procRank == 0 .AND. mod(timeJ+1, frameOut) == 0 ) then
		write(vtkFile,'(A,I0.4,A)') trim(vtkRoot), frameCounter, '.vtk'
		
		call OutputToVTKFile( mesh, vtkFile )
			
		frameCounter = frameCounter + 1
		call LogMessage(exelog, TRACE_LOGGING_LEVEL, trim(logkey)//" t = ", mesh%mesh%t)
	endif	
enddo



!--------------------------------
!	finalize : clean up
!--------------------------------

programEnd = MPI_WTIME()

write(logstring,'(A,F12.2,A)') "PROGRAM COMPLETE : elapsed time ", programEnd - programStart, " seconds."
call LogMessage(exelog, TRACE_LOGGING_LEVEL, trim(logkey)//" ", logString)

call Delete(pseSetup)
call Delete(timekeeper)
call Delete(mesh)

call MPI_FINALIZE(mpiErrCode)

contains

!> @brief Defines the initial divergence.
!> 
!> Conforms to the numberkindsmodule::scalarFnOf2DSpace interface.
!> 
!> @param[in] x
!> @param[in] y
!> @return initial divergence, @f$ \delta(x,y,0) @f$
function initDivergence( x, y )
	real(kreal) :: initDivergence
	real(kreal), intent(in) :: x, y
	initDivergence = 0.0_kreal
end function 

!> @brief Defines the initial vorticity.
!> 
!> Conforms to the numberkindsmodule::scalarFnOf2DSpace interface.
!> 
!> @param[in] x
!> @param[in] y
!> @return initial vorticiy, @f$ \zeta(x,y,0) @f$
function initVorticity( x, y)
	real(kreal) :: initVorticity
	real(kreal), intent(in) :: x, y
	!
	real(kreal) :: r2
	real(kreal), parameter :: b = 0.5_kreal
	r2 = x*x + y*y
	initVorticity = (3.0_kreal * sqrt(r2) /( r2 + ZERO_TOL * ZERO_TOL) - 2.0_kreal * b * sqrt(r2)) * r2 * exp(-b * r2)
end function 

!> @brief Defines the initial velocity.
!> 
!> Conforms to the numberkindsmodule::vectorFnOf2DSpace interface.
!> 
!> @param[in] x
!> @param[in] y
!> @return initial velocity, @f$ \vec{u}(x,y,0) @f$
function initVelocity(x, y)
	real(kreal) :: initVelocity(2)
	real(kreal), intent(in) :: x, y
	!
	real(kreal) :: utheta, r2, theta
	real(kreal), parameter :: b = 0.5_kreal
	
	r2 = x * x + y * y
	utheta = r2 * exp( - b * r2)
	theta = atan2(y,x)
	
	initVelocity(1) = - utheta * sin(theta)
	initVelocity(2) =   utheta * cos(theta)
end function 

!> @brief Defines the initial fluid depth.
!> 
!> Conforms to the numberkindsmodule::scalarFnOf2DSpace interface.
!> 
!> @param[in] x
!> @param[in] y
!> @return initial depth, @f$ h(x,y,0) @f$
function initH(x, y)
	real(kreal) :: initH
	real(kreal), intent(in) :: x, y
	initH = 1.0_kreal - bottomTopo(x,y)
end function

!> @brief Defines the bottom topography.
!> 
!> Conforms to the numberkindsmodule::scalarFnOf2DSpace interface.
!> 
!> @param[in] x
!> @param[in] y
!> @return height of bottom, @f$ h_B(x,y) @f$
function bottomTopo(x,y)
	real(kreal) :: bottomTopo
	real(kreal), intent(in) :: x, y
	!bottomTopo = 0.0_kreal
	bottomTopo = 0.5_kreal * exp( -5.0_kreal * (x-0.75_kreal) * (x-0.75_kreal) - 5.0_kreal * y * y)
end function

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

!> @brief Reads a namelist file, which must be specified on the command line at run-time execution as the first argument,
!> to define the user-specified variables for this driver program.
!> 
!> Only MPI rank 0 reads the file; it then broadcasts the relevant data to all other ranks.
!>
!> @param[in] rank MPI rank
subroutine ReadNamelistFile(rank)
	integer(kint), intent(in) :: rank
	!
	integer(kint), parameter :: initBCAST_intSize = 6
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
			bcastIntegers(6) = remeshInterval
			
			bcastReals(1) = dt
			bcastReals(2) = tfinal
			bcastReals(3) = f0
			bcastReals(4) = beta
		close(READ_UNIT)
	endif
	
	call MPI_BCAST(bcastIntegers, initBCAST_intSize, MPI_INTEGER, 0, MPI_COMM_WORLD, mpiErrCode)
	
	call MPI_BCAST(bcastReals, initBCAST_realSize, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpiErrCode)
	
	meshSeed = bcastIntegers(1)
	initNest = bcastIntegers(2)
	maxNest = bcastIntegers(3)
	amrLimit = bcastIntegers(4)
	frameOut = bcastIntegers(5)
	remeshInterval = bcastIntegers(6)
	
	dt = bcastReals(1)
	tfinal = bcastReals(2)
	f0 = bcastReals(3)
	beta = bcastReals(4)
end subroutine

end program