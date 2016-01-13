program SpherePSEConvTest

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
use PSEDirectSumModule

implicit none

include 'mpif.h'

!
! mesh variables
!
type(PolyMesh2d) :: sphere
integer(kint) :: initNest
integer(kint) :: maxNest
integer(kint) :: faceKind
integer(kint) :: meshSeed
integer(kint) :: amrLimit
real(kreal) :: radius = 1.0_kreal
namelist /meshDefine/ faceKind, initNest, amrLimit, radius
!
! pse variables
!
type(PSE) :: pseSetup
real(kreal) :: power
!
! test case variables
!
type(Field) :: constScalar
type(Field) :: constScalarInterp
type(Field) :: constScalarInterpErr
type(Field) :: constScalarLap
type(Field) :: harmonic
type(Field) :: harmonicInterp
type(Field) :: harmonicInterpError
type(Field) :: harmonicLapExact
type(Field) :: harmonicLapPSE
type(Field) :: harmonicLapError
integer(kint), parameter :: nLat = 181
integer(kint), parameter :: nLon = 360
real(kreal) :: dlam
real(kreal) :: lats(nLat)
real(kreal) :: lons(nLon)
real(kreal) :: constData(nLat, nLon)
real(kreal) :: constInterp(nLat, nLon)
real(kreal) :: constLap(nLat, nLon)
real(kreal) :: harmData(nLat, nLon)
real(kreal) :: harmInterp(nLat, nLon)
real(kreal) :: harmLap(nLat, nLon)
!
! i/o
!
character(len=MAX_STRING_LENGTH) :: outputDir
character(len=MAX_STRING_LENGTH) :: outputRoot
character(len=MAX_STRING_LENGTH) :: vtkFile
character(len=MAX_STRING_LENGTH) :: matlabFile
character(len=MAX_STRING_LENGTH) :: meshString
namelist /fileIO/ outputDir, outputRoot
integer(kint) :: writeStat
!
! computing environment / general
!
type(MPISetup) :: mpiParticles
type(MPISetup) :: mpiLongitudes
type(Logger) :: exeLog
character(len=MAX_STRING_LENGTH) :: logString
integer(kint), parameter :: logLevel = DEBUG_LOGGING_LEVEL
character(len=28) :: logKey = "SpherePSEConvTest"
integer(kint) :: mpiErrCode
real(kreal) :: timeStart, timeEnd
integer(kint) :: i, j
real(kreal), dimension(3) :: vec3


!--------------------------------
!	initialize : setup mesh and computing environment
!--------------------------------


call MPI_INIT(mpiErrCode)
call MPI_COMM_SIZE(MPI_COMM_WORLD, numProcs, mpiErrCode)
call MPI_COMM_RANK(MPI_COMM_WORLD, procRank, mpiErrCode)

call InitLogger(exeLog, procRank)

call ReadNamelistFile( procRank )

timeStart = MPI_WTIME()

!
!	mesh and spatial field setup
!
call New(sphere, meshSeed, initNest, maxNest, amrLimit, radius)
call New(constScalar, 1, sphere%particles%N_Max, "constantScalar")
call New(harmonic, 1, sphere%particles%N_Max, "sphereHarm54")

call SetScalarFieldOnMesh( sphere, constScalar, ConstantScalarFn )
call SetScalarFieldOnMesh( sphere, harmonic, SphericalHarmonicFn )

dlam = 360.0_kreal / nLon
do i = 1, nLat
	lats(i) = - 0.5_kreal * PI + (i-1) * dlam * DEG_2_RAD
enddo
do j = 1, nLon
	lons(j) = (j-1) * dlam * DEG_2_RAD
enddo


!
!	write output to vtk/paraview
!
if ( procRank == 0 ) then

	open( unit=WRITE_UNIT_1, file=vtkFile, status='REPLACE', action='WRITE', iostat=writeStat)
		if ( writeStat /= 0 ) then
			call LogMessage(exeLog, ERROR_LOGGING_LEVEL, trim(logkey)//" OutputToVTK ERROR writing to file = ", trim(vtkFile))
			return
		endif
	
		call WriteVTKPoints( sphere%particles, WRITE_UNIT_1)
		call WriteFacesToVTKPolygons( sphere%faces, WRITE_UNIT_1)
	
		call WriteVTKPointDataSectionHeader(WRITE_UNIT_1, sphere%particles%N)
		call WriteFieldToVTKPointData( constScalar, WRITE_UNIT_1)
		call WriteFieldToVTKPointData( harmonic, WRITE_UNIT_1 )
	
		call WriteFaceAreaToVTKCellData( sphere%faces, sphere%particles, WRITE_UNIT_1)
	close(WRITE_UNIT_1)


!
!	write output to matlab
!
	open( unit=WRITE_UNIT_1, file=vtkFile, status='REPLACE', action='WRITE', iostat=writeStat)
		if ( writeStat /= 0 ) then
			call LogMessage(exeLog, ERROR_LOGGING_LEVEL, trim(logkey)//" OutputToMatlab ERROR writing to file = ", trim(matlabFile))
			return
		endif
		
		call WriteToMatlab(lons, WRITE_UNIT_1, "lons")
		call WriteToMatlab(lats, WRITE_UNIT_1, "lats")
		

	close(WRITE_UNIT_1)
endif
!
!----------------
! PROGRAM END
!----------------
!
programEnd = MPI_WTIME()

write(logstring,'(A,F12.2,A)') "PROGRAM COMPLETE : elapsed time ", programEnd - programStart, " seconds."
call LogMessage(exelog, TRACE_LOGGING_LEVEL, trim(logkey)//" ", logString)


call Delete(constScalar)
call Delete(harmonic)
call Delete(sphere)

contains

subroutine SetScalarFieldOnMesh(sphere, aField, scalarFn)
	type(PolyMesh2d), intent(in) :: sphere
	type(Field), intent(inout) :: aField
	procedure(scalarFnOf3DSpace) :: scalarFn
	!
	integer(kint) :: i
	
	aField%N = sphere%particles%N
	do i = 1, sphere%particles%N
		aField%scalar(i) = scalarFn( sphere%particles%x(i), sphere%particles%y(i), sphere%particles%z(i) )
	enddo
end subroutine

pure function ConstantScalarFn( x, y, z )
	real(kreal) :: ConstantScalarFn
	real(kreal), intent(in) :: x, y, z
	ConstantScalarFn = 2.0_kreal
end function

pure function SphericalHarmonicFn( x, y, z )
	real(kreal) :: SphericalHarmonicFn
	real(kreal), intent(in) :: x, y, z
	!
	real(kreal) :: lat, lon
	
	lat = Latitude( x, y, z )
	lon = Longitude( x, y, z)
	
	SphericalHarmonicFn = 3.0_kreal * sqrt(35.0_kreal) * cos( 4.0_kreal * lon ) * sin( lat ) *  &
		(-1.0_kreal + sin( lat ) * sin( lat ) )**2
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
	integer(kint), parameter :: initBcast_intSize = 4
	integer(kint), parameter :: initBcast_realSize = 1
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
		
		if (meshSeed == ICOS_TRI_SPHERE_SEED) then
			if ( initNest == maxNest ) then
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
	
		write(vtkFile,'(5A)') trim(outputDir), '/vtkOut/', trim(outputRoot), trim(meshString), '.vtk'
		write(matlabFile,'(4)') trim(outputDir), trim(outputRoot), trim(meshString), '.m'
		
		bcastIntegers(1) = meshSeed
		bcastIntegers(2) = initNest
		bcastIntegers(3) = maxNest
		bcastIntegers(4) = amrLimit
		
		bcastReals(1) = radius
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
	
	radius = bcastReals(1)
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