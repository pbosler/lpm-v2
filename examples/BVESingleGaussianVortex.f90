program BVESingleGaussianVortexDriver
!> @file BVESingleGaussianVortex.f90
!> Driver program for the problem of a single Gaussian vortex on a rotating sphere.
!> @author Peter Bosler, Sandia National Laboratories Center for Computing Research
!>
!> Driver program for the problem of a single Gaussian vortex on a rotating sphere.@n
!> Demonstrates the solution of an inviscid incompressible flow on the sphere using @ref SphereBVE,
!> @ref SphereBVESolver, @ref Refinement, and @ref SSRFPACKRemesh.
!>
!> @image html BVEGaussVortSol.png "Late-time vorticity distribution of an initially Gaussian vortex on a rotating sphere (coarse resolution)."
!>
!> 
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
use SSRFPACKInterfaceModule
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
namelist /meshDefine/ faceKind, initNest, amrLimit

! test case variables
real(kreal), allocatable, dimension(:) :: kineticEnergy, enstrophy
real(kreal), save :: GAUSS_CONST = 0.0_kreal
real(kreal) :: shapeParam
real(kreal) :: rotRate
real(kreal) :: vortStrength
real(kreal) :: vortInitLat
real(kreal) :: vortInitLon
real(kreal), dimension(3) :: vortCenter
real(kreal) :: flowMapVarTol
real(kreal) :: circulationTol
namelist /gaussVort/ shapeParam, rotRate, vortStrength, vortInitLat, vortInitLon, flowMapVarTol, circulationTol

! remeshing
type(BVERemeshSource) :: remesh
integer(kint) :: remeshInterval
integer(kint) :: remeshCounter
type(BVEMesh) :: tempSphere

! timestepping
type(BVESolver) :: solver
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
character(len=28) :: logKey = "BVEGaussianVortex"
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

call SetSigmaFlag(1)

!
!	initialize mesh and spatial fields
!
t = 0.0_kreal
call New( sphere, meshSeed, initNest, maxNest, amrLimit, radius, rotRate )

vortCenter = radius * [ cos(vortInitLon) * cos(vortInitLat), sin(vortInitLon) * cos(vortInitLat), sin(vortInitLat)]
call SetInitialVorticityOnMesh( sphere, GaussianVortexVorticity)
GAUSS_CONST = SetGaussConst(sphere)
call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logKey)//" GAUSS_CONST = ", GAUSS_CONST)
call SetInitialVorticityOnMesh( sphere, GaussianVortexVorticity)
call LogMessage(exeLog, DEBUG_LOGGING_LEVEL, trim(logKey)//" ", "uniform mesh ready...")

AMR = (maxNest > initNest)
if ( AMR ) then 
	call LogMessage(exeLog, DEBUG_LOGGING_LEVEL, trim(logKey)//" ", "starting adaptive refinement...")
	call New(refinement, sphere%mesh%faces%N_Max)
	call SetAbsoluteTolerances( sphere, circulationTol, flowMapVarTol)
	call LogMessage(exeLog, DEBUG_LOGGING_LEVEL, trim(logKey)//" ", "AMR ready.")

	do i = 1, amrLimit
		write(logString,"(A,I2)") 'amr loop ', i
		call StartSection(exeLog, trim(logString) )
		call IterateMeshRefinementOneVariableAndFlowMap( refinement, sphere%mesh, sphere%relVort, ScalarIntegralRefinement, &
			circulationTol, "circulation refinement", flowMapVarTol, nParticlesBefore, nParticlesAfter)
		
		call SetInitialVorticityOnMesh(sphere, GaussianVortexVorticity)
		GAUSS_CONST = SetGaussConst(sphere)
		call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logKey)//" GAUSS_CONST = ", GAUSS_CONST)
		call SetInitialVorticityOnMesh(sphere, GaussianVortexVorticity)	
		call EndSection(exeLog)
	enddo
	
	call Delete(refinement)
	call LoadBalance(sphere%mpiParticles, sphere%mesh%particles%N, numProcs)
endif

call SetStreamFunctionsOnMesh(sphere)
call SetVelocityOnMesh(sphere)

call AddTracers(sphere, 1, [1])
sphere%tracers(1)%name = "initLat"
sphere%tracers(1)%units = "radians"
do i = 1, sphere%mesh%particles%N
	vec3 = PhysCoord(sphere%mesh%particles, i)
	call InsertScalarToField(sphere%tracers(1), Latitude(vec3))
enddo

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

!$acc enter data create(solver)
call New(solver, sphere)



nTimesteps = floor(tfinal / dt)
allocate(kineticEnergy(nTimesteps+1))
allocate(enstrophy(nTimesteps+1))
kineticEnergy(1) = TotalKE(sphere)
enstrophy(1) = TotalEnstrophy(sphere)
!--------------------------------
!	run : evolve the problem in time 
!--------------------------------

call LogMessage(exeLog, DEBUG_LOGGING_LEVEL, trim(logkey)//" ", "starting timestepping loop.")
do timeJ = 0, nTimesteps - 1

	if ( mod(timeJ+1, remeshInterval) == 0 ) then
		call New(remesh, sphere)
	
		call New( tempSphere, meshSeed, initNest, maxNest, amrLimit, radius, rotRate )
		call AddTracers(tempSphere, 1, [1])
		tempSphere%tracers(1)%name = "initLat"
		tempSphere%tracers(1)%units = "radians"
		
		call LagrangianRemeshBVEWithVorticityFunction(remesh, sphere, tempSphere, AMR, GaussianVortexVorticity, &
					scalarIntegralRefinement, circulationTol, "circulation refinement", &
					RefineFlowMapYN = .TRUE., flowMapVarTol = flowMapVarTol, nLagTracers = 1, tracerFn1 = InitLatTracer )
		
		call Copy(sphere, tempSphere)
		
		remeshCounter = remeshCounter + 1
		
		call Delete(tempSphere)
		call Delete(remesh)
		
		call Delete(solver)
		!$acc exit data delete(solver)
		
		
		call New(solver, sphere)
	endif
	
	call Timestep(solver, sphere, dt)
	
	t = real(timeJ +1, kreal) * dt
	sphere%mesh%t = t
	
	enstrophy(timeJ+2) = TotalEnstrophy(sphere)
	kineticEnergy(timeJ+2) = TotalKE(sphere)
	
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
call Delete(sphere)
if (AMR) call Delete(refinement)

call MPI_FINALIZE(mpiErrCode)

contains

!> @brief Returns the latitude of a Lagrangian coordinate.  
!> 
!> Conforms to numberkindsmodule::scalarFnOf3DSpace interface.
pure function InitLatTracer( x0, y0, z0 )
	real(kreal) :: InitLatTracer
	real(kreal), intent(in) :: x0, y0, z0
	InitLatTracer = Latitude( x0, y0, z0 )
end function

!> @brief Sets absolute tolerances for AMR criteria based on an initially uniform mesh with corresponding vorticity distribution.
!> @param[in] aBVEMesh @ref SphereBVE mesh
!> @param[inout] circTol On input, relative tolerance (between 0 and 1).  On output, absolute tolerance for circulation magnitude about a face.
!> @param[inout] lagVarTol On input, relative tolerance (between 0 and 1).  On output, absolute tolerance for Lagrangian variation magnitude on each face.
subroutine SetAbsoluteTolerances( aBVEMesh, circTol, lagVarTol )
	type(BVEMesh), intent(in) :: aBVEMesh
	real(kreal), intent(inout) :: circTol
	real(kreal), intent(inout) :: lagVarTol
	!
	circTol = circTol * MaxCirculationMagnitudePerFace(aBVEMesh)
	lagVarTol = lagVarTol * MaxLagrangianVariationPerFace(aBVEMesh%mesh)
end subroutine

!> @brief Vorticity distribution function, used to define initial conditions and for indirect vorticity interpolation
!> to @f$ t = 0 @f$.
!>
!> Conforms to the numberkindsmodule::scalarFnOf3DSpace interface.
!> 
!> @param[in] x
!> @param[in] y
!> @param[in] z 
!> @return initial vorticity value at at location (x,y,z)
function GaussianVortexVorticity( x, y, z )
	real(kreal) :: GaussianVortexVorticity
	real(kreal), intent(in) :: x, y, z
	GaussianVortexVorticity = vortStrength * exp( - shapeParam * shapeParam * &
		(radius * radius - x * vortCenter(1) - y * vortCenter(2) - z * vortCenter(3)) ) - GAUSS_CONST
end function

!> @brief Defines a constant to ensure that the total integral of vorticity is zero over the whole sphere.
!> 
!> @param[in] aBVEMesh @ref SphereBVE mesh
!> @return C such that 
!> @f[
!> 		\int_S \zeta(\vec{x})\,dA = \int_{S} \zeta_G(\vec{x})\,dA - C = 0,
!> @f]
!> where @f$ \zeta_G @f$ is a Gaussian distribution on a sphere.
function SetGaussConst( aBVEMesh )
	real(kreal) :: SetGaussConst
	type(BVEMesh), intent(in) :: aBVEMesh
	SetGaussConst = sum( aBVEMesh%relVort%scalar(1:aBVEMesh%mesh%particles%N) * &
						 aBVEMesh%mesh%particles%area(1:aBVEMesh%mesh%particles%N), &
						 MASK=aBVEMesh%mesh%particles%isActive(1:aBVEMesh%mesh%particles%N) ) / &
						 (4.0_kreal * PI * aBVEMesh%radius )
end function

!> @brief Reads a namelist file, which must be specified on the command line at run-time execution as the first argument,
!> to define the user-specified variables for this driver program.
!> 
!> Only MPI rank 0 reads the file; it then broadcasts the relevant data to all other ranks.
!>
!> @param[in] rank MPI rank
subroutine ReadNamelistFile( rank )
	integer(kint), intent(in) :: rank
	!
	character(len=MAX_STRING_LENGTH) :: namelistFilename
	integer(kint), parameter :: initBcast_intSize = 6
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
		bcastReals(3) = shapeParam
		bcastReals(4) = rotRate
		bcastReals(5) = vortStrength
		bcastReals(6) = vortInitLon
		bcastReals(7) = vortInitLat
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
	
	meshSeed = bcastIntegers(1)
	initNest = bcastIntegers(2)
	maxNest = bcastIntegers(3)
	amrLimit = bcastIntegers(4)
	frameOut = bcastIntegers(5)
	remeshInterval = bcastIntegers(6)
	
	dt = bcastReals(1)
	tfinal = bcastReals(2)
	shapeParam = bcastReals(3)
	rotRate = bcastReals(4)
	vortStrength = bcastReals(5)
	vortInitLon = bcastReals(6)
	vortInitLat = bcastReals(7)
	circulationTol = bcastReals(8)
	flowMapVarTol = bcastReals(9)
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