program CollidingDipolesDriver
!> @file CollidingDipoles.f90
!> Driver program for the problem of two colliding Lamb dipoles in the plane.
!> @author Peter Bosler, Sandia National Laboratories Center for Computing Research
!>
!> Driver program for the problem of two colliding Lamb dipoles in the plane.@n
!> Demonstrates the solution of an inviscid incompressible flow in the plane using @ref PlanarIncompressible,
!> @ref PlanarIncompressibleSolver, @ref Refinement, and @ref BIVARRemesh.
!>
!> @image html collidingDipolesInitCond.png "Initial vorticity distribution for colliding dipoles problem."
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
use PlaneGeomModule
use PlanarIncompressibleModule
use PlanarIncompressibleSolverModule
use RefinementModule
use BIVARRemeshModule

implicit none

include 'mpif.h'

! mesh variables
type(PlaneMeshIncompressible) :: plane 
integer(kint) :: initNest
integer(kint) :: maxNest
integer(kint) :: amrLimit
integer(kint) :: faceKind
integer(kint) :: meshSeed
real(kreal) :: meshRadius

! AMR variables
type(RefineSetup) :: refinement
logical(klog) :: useAMR
logical(klog) :: doFlowMapRefine
integer(kint) :: nParticlesBefore
integer(kint) :: nParticlesAfter

namelist /meshDefine/ faceKind, initNest, maxNest, amrLimit, meshRadius

! test case variables
real(kreal), allocatable, dimension(:) :: kineticEnergy, enstrophy
real(kreal) :: radius1
real(kreal) :: initX1
real(kreal) :: initY1
real(kreal) :: initStrength1
real(kreal) :: radius2
real(kreal) :: initX2
real(kreal) :: initY2
real(kreal) :: initStrength2
real(kreal) :: circulationTol
real(kreal) :: flowMapVarTol
real(kreal), parameter :: LAMB_K0 = 3.8317_kreal

namelist /twoDipoles/ radius1, initX1, initY1, initStrength1, radius2, initX2, initY2, initStrength2, &
					  circulationTol, flowMapVarTol

! remeshing variables
integer(kint) :: remeshInterval
integer(kint) :: resetLagParamInterval
integer(kint) :: remeshCounter
type(PlaneMeshIncompressible) :: tempMesh
type(PlaneMeshIncompressible) :: refMesh

! timestepping variables
type(PlaneSolver) :: solver
real(kreal) :: t
real(kreal) :: dt
real(kreal) :: tfinal
integer(kint) :: timeJ
integer(kint) :: nTimesteps

namelist /timestepping/ dt, tfinal, remeshInterval, resetLagParamInterval

! i/o variables
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
logical(klog), save :: logInit = .FALSE.
character(len=MAX_STRING_LENGTH) :: logString
integer(kint), parameter :: logLevel = DEBUG_LOGGING_LEVEL
character(len=28) :: logKey = "twoDipoles"
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

!
!	initialize mesh and spatial fields
!
call New( plane, initNest, maxNest, meshSeed, amrLimit, meshRadius)
call SetInitialVorticityOnMesh( plane, TwoDipolesVorticity )

useAMR = ( amrLimit > 0 )
if ( useAMR ) then

	doFlowMapRefine = .FALSE.
	
	call New(refinement, plane%mesh%faces%N_Max)
	
	call SetAbsoluteTolerances( plane, circulationTol, flowMapVarTol )
	call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logKey)//" circulation tol = ", circulationTol)
	call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logKey)//" flowMapVarTol = ", flowMapVarTol)
	
	do i = 1, amrLimit
		call IterateMeshRefinementOneVariableAndFlowMap( refinement, plane%mesh, plane%vorticity, &
			ScalarIntegralRefinement, circulationTol, "circulation refinement", flowMapVarTol, &
			nParticlesBefore, nParticlesAfter)
		call SetInitialVorticityOnMesh(plane, TwoDipolesVorticity)
	enddo
	
	call LoadBalance( plane%mpiParticles, plane%mesh%particles%N, numProcs)
	
	call Delete(refinement)
endif

call AddTracers( plane, 1, [2])
plane%tracers(1)%name = "initX"
call StoreLagParamAsTracer( plane, 1 )
call SetStreamFunctionOnMesh( plane )
call SetVelocityOnMesh( plane )

!
!	output initial state
!
frameCounter = 0
t = 0.0_kreal
if ( procRank == 0 ) then
	call LogStats(plane, exeLog)
	
	call MeshSeedString( meshString, meshSeed)
	
	if ( useAMR ) then
		write( meshString, '(A,I1,A,I1,A)') trim(meshString)//"_AMR_", initNest, 'to', initNest+amrLimit, '_'
	else
		write( meshString, '(A,I1,A)') trim(meshString)//'_', initNest, '_'
	endif
	
	write(vtkRoot,'(4A)') trim(outputDir), '/vtkOut/', trim(outputRoot), trim(meshString)
	write(vtkFile,'(A,I0.4,A)') trim(vtkRoot), frameCounter, '.vtk'
	write(vtkRoot,'(4A)') trim(outputDir), '/vtkOut/', trim(outputRoot), trim(meshString)
	write(vtkFile,'(A,I0.4,A)') trim(vtkRoot), frameCounter, '.vtk'
	
	call OutputToVTK(plane, vtkFile )
	
	frameCounter = frameCounter + 1
	
	call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logkey)//" t = ", plane%mesh%t)
endif

!
!	initialize timestepping
!
call New(solver, plane)
nTimesteps = floor(tfinal/dt)
allocate(kineticEnergy(nTimesteps+1))
allocate(enstrophy(nTimesteps+1))
kineticEnergy(1) = TotalKE(plane)
enstrophy(1) = TotalEnstrophy(plane)
remeshCounter = 0

!
!	initialize remeshing
!
!call New(tempMesh, initNest, maxNest, meshSeed, amrLimit, meshRadius)
!call AddTracers( tempMesh, 1, [2])
!tempMesh%tracers(1)%name = "initX"
!
!call New( refMesh, initNest, maxNest, meshSeed, amrLimit, meshRadius )
!call AddTracers( refMesh, 1, [2])
!refMesh%tracers(1)%name = "initX"

!--------------------------------
!	run : evolve the problem in time 
!--------------------------------

do timeJ = 0, nTimesteps - 1
	!
	! remesh
	!
	if ( mod( timeJ+1, remeshInterval) == 0 ) then
		remeshCounter = remeshCounter + 1
		!
		!	choose remeshing procedure
		!
		if ( remeshCounter < resetLagParamInterval ) then
			!
			!	remesh to t = 0
			!
			call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logKey)//" remeshing to t = 0, remeshCounter = ", remeshCounter)
			
			call New(tempMesh, initNest, maxNest, meshSeed, amrLimit, meshRadius)
			call AddTracers( tempMesh, 1, [2])

			call LagrangianRemeshPlanarIncompressibleWithVorticityFunction( tempMesh, plane, TwoDipolesVorticity, &
					ScalarIntegralRefinement, circulationTol, "circulation refinement", RefineFlowMapYN = doFlowMapRefine, &
					flowMapVarTol = flowMapVarTol )
			call StoreLagParamAsTracer( tempMesh, 1)
		
			call Copy( plane, tempMesh)
			call Delete(tempMesh)
		elseif ( remeshCounter == resetLagParamInterval ) then
			!
			!	remesh to t = 0, keep old mesh as reference defined for t = t_remesh
			!
			call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logKey)//" remeshing to t = 0, remeshCounter = ", remeshCounter)
			
			t = plane%mesh%t
			call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logKey)//" setting new ref_t = ", t)
			
			call New( refMesh, initNest, maxNest, meshSeed, amrLimit, meshRadius )
			call AddTracers( refMesh, 1, [2])
			refMesh%tracers(1)%name = "initX"
			
			call LagrangianRemeshPlanarIncompressibleWithVorticityFunction( refMesh, plane, TwoDipolesVorticity, &
					ScalarIntegralRefinement, circulationTol, "circulation refinement", RefineFlowMapYN = doFlowMapRefine, &
					flowMapVarTol = flowMapVarTol )
			call StoreLagParamAsTracer(refMesh, 1)
			
			call Copy(plane, refMesh)
			refMesh%mesh%t = t		
			
			call ResetLagrangianParameter(refMesh%mesh%particles)
			call ResetLagrangianParameter(plane%mesh%particles)
		elseif ( remeshCounter > resetLagParamInterval .AND. mod(remeshCounter, resetLagParamInterval) == 0 ) then
			!
			!	remesh to previous reference time, then keep current mesh as new reference
			!
			call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logKey)//" remeshing to t = ", refMesh%mesh%t)
			t = plane%mesh%t
			call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logKey)//" setting new ref_t = ", t)
				
			call New(tempMesh, initNest, maxNest, meshSeed, amrLimit, meshRadius )
			call AddTracers( tempMesh, 1, [2])
			tempMesh%tracers(1)%name = "initX"
			
			call LogMessage(exeLog, DEBUG_LOGGING_LEVEL, trim(logKey)//" calling remesh routine = ", t)
			
			call LagrangianRemeshPlanarIncompressibleToReferenceMesh( tempMesh, plane, refMesh, ScalarIntegralRefinement, &
					circulationTol, "circulation refinement", RefineFlowMapYN = doFlowMapRefine, &
					flowMapVarTol = flowMapVarTol)
			
			call LogMessage(exeLog, DEBUG_LOGGING_LEVEL, trim(logKey)//" returned from remesh = ", t)
			
			call Copy(plane, tempMesh)
			call Copy(refMesh, tempMesh)
			refMesh%mesh%t = t
			
			call ResetLagrangianParameter(refMesh%mesh%particles)
			call ResetLagrangianParameter(plane%mesh%particles)
			call Delete(tempMesh)
		else
			!
			!	remesh to reference time
			!
			call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logKey)//" remeshing to t = ", refMesh%mesh%t)
			call New(tempMesh, initNest, maxNest, meshSeed, amrLimit, meshRadius )
			call AddTracers( tempMesh, 1, [2])
			tempMesh%tracers(1)%name = "initX"
			
			call LagrangianRemeshPlanarIncompressibleToReferenceMesh( tempMesh, plane, refMesh, ScalarIntegralRefinement, &
					circulationTol, "circulation refinement", RefineFlowMapYN = doFlowMapRefine, &
					flowMapVarTol = flowMapVarTol)
					
			call Copy(plane, tempMesh)
			call Delete(tempMesh)
		endif
		
		
		call Delete(solver)
		call New(solver, plane)
	endif
	
	!
	! advance time step
	!
	call Timestep(solver, plane, dt)
	plane%mesh%t = real( timeJ + 1, kreal) * dt
	
	kineticEnergy(timeJ + 2) = TotalKE(plane)
	enstrophy(timeJ + 2) = TotalEnstrophy(plane)
	
	if ( procRank == 0 .AND. mod(timeJ + 1, frameOut) == 0 ) then
		write(vtkFile,'(A,I0.4,A)') trim(vtkRoot), frameCounter, '.vtk'
		call OutputToVTK(plane, vtkFile)
		frameCounter = frameCounter + 1
		call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logKey)//" t = ", plane%mesh%t)
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
call Delete(refMesh)
call Delete(solver)
call Delete(plane)

call MPI_FINALIZE(mpiErrCode)

contains

!> @brief Stores the Lagrangian coordinate of each particle as a passive tracer.
!> Helpful for visualization of transport
!>
!> @param[inout] aPlane @ref PlanarIncompressible mesh
!> @param[in] tracerID index to `PlanarIncompressible%%tracers(:)` array to store the Lagrangian coordinate
subroutine StoreLagParamAsTracer( aPlane, tracerID )
	type(PlaneMeshIncompressible), intent(inout) :: aPlane
	integer(kint), intent(in) :: tracerID
	!
	integer(kint) :: i
	
	do i = 1, aPlane%mesh%particles%N
		aPlane%tracers(tracerID)%xComp(i) = aPlane%mesh%particles%x0(i)
		aPlane%tracers(tracerID)%yCOmp(i) = aPlane%mesh%particles%y0(i)
	enddo
end subroutine

!> @brief Vorticity distribution function, used to define initial conditions and for indirect vorticity interpolation
!> to @f$ t = 0 @f$.
!>
!> Conforms to the numberkindsmodule::scalarFnOf2DSpace interface.
!> 
!> @param[in] x
!> @param[in] y
!> @return vorticity value at initial time at location (x,y)
function TwoDipolesVorticity( x, y ) 
	real(kreal) :: TwoDipolesVorticity
	real(kreal), intent(in) :: x, y
	
	TwoDipolesVorticity = LambDipole( x, y, radius1, initStrength1, initX1, initY1 ) + &
						  LambDipole( x, y, radius2, initStrength2, initX2, initY2 )
end function

!> @brief Vorticity function associated with a single Lamb dipole at location `(xcent, ycent)`
!> 
!> @param[in] x x-coordinate of output vorticity
!> @param[in] y y-coordinate of output vorticity
!> @param[in] lambR radius of dipole
!> @param[in] U0 strength parameter 
!> @param[in] xcent x-coordinate of dipole center
!> @param[in] ycent y-coordinate of dipole center
function LambDipole( x, y, lambR, U0, xcent, ycent)
	real(kreal) ::LambDipole
	real(kreal), intent(in) :: x, y, xCent, yCent, lambR, U0
	!
	real(kreal) :: r, k, sintheta, denom
	real(kreal) :: BESSJ0, BESSJ1
	external BESSJ1, BESSJ0

	r = sqrt( (x - xcent) * (x - xcent) + (y - ycent) * (y - ycent) )

	if ( r > lambR .OR. r < ZERO_TOL) then
		LambDipole = 0.0_kreal
	else
		k = LAMB_K0 / lambR

		sintheta = y / r

		denom = BESSJ0( LAMB_K0 )

		LambDipole = -2.0_kreal * U0 * k * BESSJ1( k * r) * sintheta / denom
	endif
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
	integer(kint), parameter :: initBcast_intSize = 7
	integer(kint), parameter :: initBcast_realSize = 13
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
			read(READ_UNIT, nml=TwoDipoles)
			rewind(READ_UNIT)
			read(READ_UNIT, nml=timestepping)
			rewind(READ_UNIT)
			read(READ_UNIT, nml=fileIO)
		close(READ_UNIT)
		
		if ( faceKind == 3 ) then
			meshSeed = TRI_HEX_SEED
		elseif ( faceKind == 4) then
			meshSeed = QUAD_RECT_SEED
		else
			call LogMessage(exeLog, WARNING_LOGGING_LEVEL, trim(logkey)//" ReadNamelistFile WARNING : ", &
				" invalid faceKind -- using triangles.")
			meshSeed = QUAD_RECT_SEED
		endif
		
		bcastIntegers(1) = initNest
		bcastIntegers(2) = maxNest
		bcastIntegers(3) = amrLimit
		bcastIntegers(4) = frameOut
		bcastIntegers(5) = remeshInterval
		bcastIntegers(6) = meshSeed
		bcastIntegers(7) = resetLagParamInterval
		
		bcastReals(1) = dt
		bcastReals(2) = tfinal
		bcastReals(3) = radius1
		bcastReals(4) = initX1
		bcastReals(5) = initY1
		bcastReals(6) = initStrength1
		bcastReals(7) = radius2
		bcastReals(8) = initX2
		bcastReals(9) = initY2
		bcastReals(10) = initStrength2	
		bcastReals(11) = circulationTol
		bcastReals(12) = flowMapVarTol
		bcastReals(13) = meshRadius
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
	meshSeed = bcastIntegers(6)
	resetLagParamInterval = bcastIntegers(7)
	
	dt = bcastReals(1)
	tfinal = bcastReals(2)
	radius1 = bcastReals(3)
	initX1 = bcastReals(4)
	initY1 = bcastReals(5)
	initStrength1 = bcastReals(6)
	radius2 = bcastReals(7)
	initX2 = bcastReals(8)
	initY2 = bcastReals(9)
	initStrength2 = bcastReals(10)
	circulationTol = bcastReals(11)
	flowMapVarTol = bcastReals(12)
	meshRadius = bcastReals(13)
end subroutine


!> @brief Initializes a @ref Logger for this executable program.
!> 
!> Output is controlled by message priority level and MPI rank.
!> 
!> @param[in] aLog @ref Logger to initialize
!> @param[in] rank MPI rank
subroutine InitLogger(aLog,rank)
	type(Logger), intent(out) :: aLog
	integer(kint), intent(in) :: rank
	write(logKey,'(A,A,I0.2,A)') trim(logKey),'_',rank,' : '
	if ( rank == 0 ) then
		call New(aLog,logLevel)
	else
		call New(aLog,ERROR_LOGGING_LEVEL)
	endif
	logInit = .TRUE.
end subroutine

end program


!>    @brief This function calculates the first kind modified Bessel function
!>     of integer order N, for any REAL X. 
!>
!>     We use here the classical
!>     recursion formula, when X > N. For X < N, the Miller's algorithm
!>     is used to avoid overflows.
!>
!>     REFERENCE:
!>     C.W.CLENSHAW, CHEBYSHEV SERIES FOR MATHEMATICAL FUNCTIONS,
!>     MATHEMATICAL TABLES, VOL.5, 1962.
!> 
!> @param[in] n order of Bessel function
!> @param[in] x 
!> @return Bessel function output, @f$ B_n(x) @f$
FUNCTION BESSJ (N,X)

      PARAMETER (IACC = 40,BIGNO = 1.D10, BIGNI = 1.D-10)
      REAL *8 X,BESSJ,BESSJ0,BESSJ1,TOX,BJM,BJ,BJP,SUM
      IF (N.EQ.0) THEN
      BESSJ = BESSJ0(X)
      RETURN
      ENDIF
      IF (N.EQ.1) THEN
      BESSJ = BESSJ1(X)
      RETURN
      ENDIF
      IF (X.EQ.0.) THEN
      BESSJ = 0.
      RETURN
      ENDIF
      print*,'here'
      TOX = 2./X
      IF (X.GT.FLOAT(N)) THEN
      BJM = BESSJ0(X)
      BJ  = BESSJ1(X)
      DO 11 J = 1,N-1
      BJP = J*TOX*BJ-BJM
      BJM = BJ
      BJ  = BJP
   11 CONTINUE
      BESSJ = BJ
      ELSE
      M = 2*((N+INT(SQRT(FLOAT(IACC*N))))/2)
      BESSJ = 0.
      JSUM = 0
      SUM = 0.
      BJP = 0.
      BJ  = 1.
      DO 12 J = M,1,-1
      BJM = J*TOX*BJ-BJP
      BJP = BJ
      BJ  = BJM
      IF (ABS(BJ).GT.BIGNO) THEN
      BJ  = BJ*BIGNI
      BJP = BJP*BIGNI
      BESSJ = BESSJ*BIGNI
      SUM = SUM*BIGNI
      ENDIF
      IF (JSUM.NE.0) SUM = SUM+BJ
      JSUM = 1-JSUM
      IF (J.EQ.N) BESSJ = BJP
   12 CONTINUE
      SUM = 2.*SUM-BJ
      BESSJ = BESSJ/SUM
      ENDIF
      RETURN
      END


!>  @brief   This function calculates the First Kind Bessel Function of
!>     order 0, for any real number X. 
!> 
!> The polynomial approximation by
!> series of Chebyshev polynomials is used for 0<X<8 and 0<8/X<1.
!> 
!>     REFERENCES:
!>     M.ABRAMOWITZ,I.A.STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS, 1965.
!>     C.W.CLENSHAW, NATIONAL PHYSICAL LABORATORY MATHEMATICAL TABLES,
!>     VOL.5, 1962.
!>
!> @param[in] x 
!> @return Bessel function output, @f$ B_0(x) @f$
      FUNCTION BESSJ0 (X)
      REAL *8 X,BESSJ0,AX,FR,FS,Z,FP,FQ,XX



      REAL *8 Y,P1,P2,P3,P4,P5,R1,R2,R3,R4,R5,R6  &
               ,Q1,Q2,Q3,Q4,Q5,S1,S2,S3,S4,S5,S6
      DATA P1,P2,P3,P4,P5 /1.D0,-.1098628627D-2,.2734510407D-4, &
      -.2073370639D-5,.2093887211D-6 /
      DATA Q1,Q2,Q3,Q4,Q5 /-.1562499995D-1,.1430488765D-3, &
      -.6911147651D-5,.7621095161D-6,-.9349451520D-7 /
      DATA R1,R2,R3,R4,R5,R6 /57568490574.D0,-13362590354.D0, &
      651619640.7D0,-11214424.18D0,77392.33017D0,-184.9052456D0 /
      DATA S1,S2,S3,S4,S5,S6 /57568490411.D0,1029532985.D0, &
      9494680.718D0,59272.64853D0,267.8532712D0,1.D0 /
      IF(X.EQ.0.D0) GO TO 1
      AX = ABS (X)
      IF (AX.LT.8.) THEN
      Y = X*X
      FR = R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6))))
      FS = S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6))))
      BESSJ0 = FR/FS
      ELSE
      Z = 8./AX
      Y = Z*Z
      XX = AX-.785398164
      FP = P1+Y*(P2+Y*(P3+Y*(P4+Y*P5)))
      FQ = Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)))
      BESSJ0 = SQRT(.636619772/AX)*(FP*COS(XX)-Z*FQ*SIN(XX))
      ENDIF
      RETURN
    1 BESSJ0 = 1.D0
      RETURN
      END
! ---------------------------------------------------------------------------

!>  @brief   This function calculates the First Kind Bessel Function of
!>     order 1, for any real number X. 
!> 
!> The polynomial approximation by
!> series of Chebyshev polynomials is used for 0<X<8 and 0<8/X<1.
!> 
!>     REFERENCES:
!>     M.ABRAMOWITZ,I.A.STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS, 1965.
!>     C.W.CLENSHAW, NATIONAL PHYSICAL LABORATORY MATHEMATICAL TABLES,
!>     VOL.5, 1962.
!>
!> @param[in] x 
!> @return Bessel function output, @f$ B_1(x) @f$
      FUNCTION BESSJ1 (X)
      REAL *8 X,BESSJ1,AX,FR,FS,Z,FP,FQ,XX
      REAL *8 Y,P1,P2,P3,P4,P5,P6,R1,R2,R3,R4,R5,R6  &
               ,Q1,Q2,Q3,Q4,Q5,S1,S2,S3,S4,S5,S6
      DATA P1,P2,P3,P4,P5 /1.D0,.183105D-2,-.3516396496D-4,  &
      .2457520174D-5,-.240337019D-6 /,P6 /.636619772D0 /
      DATA Q1,Q2,Q3,Q4,Q5 /.04687499995D0,-.2002690873D-3,   &
      .8449199096D-5,-.88228987D-6,.105787412D-6 /
      DATA R1,R2,R3,R4,R5,R6 /72362614232.D0,-7895059235.D0, &
      242396853.1D0,-2972611.439D0,15704.48260D0,-30.16036606D0 /
      DATA S1,S2,S3,S4,S5,S6 /144725228442.D0,2300535178.D0, &
      18583304.74D0,99447.43394D0,376.9991397D0,1.D0 /

      AX = ABS(X)
      IF (AX.LT.8.) THEN
      Y = X*X
      FR = R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6))))
      FS = S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6))))
      BESSJ1 = X*(FR/FS)
      ELSE
      Z = 8./AX
      Y = Z*Z
      XX = AX-2.35619491
      FP = P1+Y*(P2+Y*(P3+Y*(P4+Y*P5)))
      FQ = Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)))
      BESSJ1 = SQRT(P6/AX)*(COS(XX)*FP-Z*SIN(XX)*FQ)*SIGN(S6,X)
      ENDIF
      RETURN
      END