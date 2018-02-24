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
namelist /meshDefine/ faceKind, initNest, amrLimit, radius, psePower
!
! pse variables
!
type(PSE) :: pseSetup
real(kreal) :: psePower
!
! test case variables
!
type(Field) :: constScalar
type(Field) :: constScalarInterp
type(Field) :: constScalarInterpErr
type(Field) :: constScalarLap
real(kreal), parameter :: constVal = 2.0_kreal
type(Field) :: harmonic
type(Field) :: harmonicInterp
type(Field) :: harmonicInterpError
type(Field) :: harmonicLapExact
type(Field) :: harmonicLapPSE
type(Field) :: harmonicLapError
type(Field) :: streamFnScalar
integer(kint), parameter :: nLat = 181
integer(kint), parameter :: nLon = 360
real(kreal) :: dlam
real(kreal) :: lats(nLat)
real(kreal) :: lons(nLon)
real(kreal) :: constData(nLat, nLon)
real(kreal) :: constInterp(nLat, nLon)
real(kreal) :: constInterpErr(nLat,nLon)
real(kreal) :: constLap(nLat, nLon)
real(kreal) :: harmData(nLat, nLon)
real(kreal) :: harmInterp(nLat, nLon)
real(kreal) :: harmInterpError(nLat, nLon)
real(kreal) :: harmLap(nLat, nLon)
real(kreal) :: harmLapExact(nLat,nLon)
real(kreal) :: harmLapError(nLat,nLon)
real(kreal) :: stream(nLat,nLon)
real(kreal) ::unifLinfConst
real(kreal) ::unifL2Const
real(kreal) ::unifLinfConstLap
real(kreal) ::unifL2ConstLap
real(kreal) ::unifLinfHarm
real(kreal) ::unifL2Harm
real(kreal) ::unifLinfHarmLap
real(kreal) ::unifL2HarmLap
real(kreal) ::particlesLinfConst
real(kreal) ::particlesL2Const
real(kreal) ::particlesLinfConstLap
real(kreal) ::particlesL2ConstLap
real(kreal) ::particlesLinfHarm
real(kreal) ::particlesL2Harm
real(kreal) ::particlesLinfHarmLap
real(kreal) ::particlesL2HarmLap
real(kreal) :: harmDenom
real(kreal) :: harmLapDenom

!
! i/o
!
character(len=MAX_STRING_LENGTH) :: testTitle
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
integer(kint) :: i, j, k
real(kreal), dimension(3) :: xi, xj
real(kreal) :: xx, yy, zz, pseKin, lapKernel, greensKernel
real(kreal), parameter :: fourPi = 4.0_kreal * PI

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
if ( procRank == 0 ) then
	write(testTitle,'(A,X,F9.6)') trim(testTitle), MaxEdgeLength(sphere%edges, sphere%particles)
endif
call StartSection(exeLog, trim(testTitle))

call New(constScalar, 1, sphere%particles%N_Max, "constantScalar")
call New(harmonic, 1, sphere%particles%N_Max, "harmonicExact")
call New(harmonicLapExact, 1, sphere%particles%N_Max, "harmLapExact")
call New(streamFnScalar,1,sphere%particles%N_Max, "streamFn")
call SetScalarFieldOnMesh( sphere, constScalar, ConstantScalarFn )
call SetScalarFieldOnMesh( sphere, harmonic, SphericalHarmonicFn )
call SetScalarFieldOnMesh( sphere, harmonicLapExact, ExactHarmonicLaplacian )

call New(constScalarLap, 1, sphere%particles%N, "constScalarLap")
call New(harmonicLapError, 1, sphere%particles%N_Max, "harmLapError")
call New(harmonicLapPSE, 1, sphere%particles%N_Max, "harmLapPSE")

call New(constScalarInterp, 1, sphere%particles%N, "constScalarInterp")
call New(harmonicInterp, 1, sphere%particles%N, "harmonicInterp")
call New(harmonicInterpError, 1, sphere%particles%N, "harmonicInterpError")


dlam = 360.0_kreal / nLon
do i = 1, nLat
	lats(i) = - 0.5_kreal * PI + (i-1) * dlam * DEG_2_RAD
enddo
do j = 1, nLon
	lons(j) = (j-1) * dlam * DEG_2_RAD
enddo

do j = 1, nLon 
	do i = 1, nLat
		xx = radius * cos( lons(j) ) * cos(lats(i))
		yy = radius * sin( lons(j) ) * cos(lats(i))
		zz = radius * sin( lats(i) )		
		constData(i,j) = constVal
		harmData(i,j) = SphericalHarmonicFn( xx, yy, zz)
		harmLapExact(i,j) = ExactHarmonicLaplacian( xx, yy, zz)
	enddo
enddo

!--------------------------------
!	run : perform interpolation tests
!--------------------------------
call New( mpiLongitudes, nLon, numProcs )
call New( mpiParticles, sphere%particles%N, numProcs)

call New(pseSetup, sphere, psepower)

!
!	interpolation to uniform lat-lon grid
!
constLap = 0.0_kreal
harmLap = 0.0_kreal
stream = 0.0_kreal
do j = mpiLongitudes%indexStart(procRank), mpiLongitudes%indexEnd(procRank)
	do i = 1, nLat
		xx = radius * cos( lons(j) ) * cos(lats(i))
		yy = radius * sin( lons(j) ) * cos(lats(i))
		zz = radius * sin( lats(i) )
		constInterp(i,j) = PSESphereInterpolateScalar( pseSetup, sphere, constScalar, [xx, yy, zz])
		harmInterp(i,j) = PSESphereInterpolateScalar( pseSetup, sphere, harmonic, [xx, yy, zz] )
		
		do k = 1, sphere%particles%N
			if ( sphere%particles%isActive(k) ) then
				pseKin = SphereDistance( PhysCoord(sphere%particles, k), [xx, yy, zz]) / pseSetup%eps
				lapKernel = bivariateLaplacianKernel8( pseKin ) / pseSetup%eps**2
				greensKernel = log( radius * radius - xx * sphere%particles%x(k) - yy * sphere%particles%y(k) - &
					zz * sphere%particles%z(k) ) / fourPi
							
				constLap(i,j) = constLap(i,j) + lapKernel * ( constScalar%scalar(k) - constData(i,j) ) * &
					sphere%particles%area(k) / pseSetup%eps**2
				harmLap(i,j) = harmLap(i,j) + lapKernel * ( harmonic%scalar(k) - harmData(i,j) ) * &
					sphere%particles%area(k) / pseSetup%eps**2
				stream(i,j) = stream(i,j) + greensKernel * harmonicLapExact%scalar(k) * sphere%particles%area(k)
			endif
		enddo
	enddo
enddo

do i = 0, numProcs - 1
	call MPI_BCAST( constInterp(:, mpiLongitudes%indexStart(i):mpiLongitudes%indexEnd(i)), &
		nLat * mpiLongitudes%messageLength(i), MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode )
	call MPI_BCAST( harmInterp(:, mpiLongitudes%indexStart(i):mpiLongitudes%indexEnd(i)), &
		nLat * mpiLongitudes%messageLength(i), MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode )
	call MPI_BCAST( constLap(:, mpiLongitudes%indexStart(i):mpiLongitudes%indexEnd(i)), &
		nLat * mpiLongitudes%messageLength(i), MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode )
	call MPI_BCAST( harmLap(:, mpiLongitudes%indexStart(i):mpiLongitudes%indexEnd(i)), &
		nLat * mpiLongitudes%messageLength(i), MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode )
	call MPI_BCAST( stream(:, mpiLongitudes%indexStart(i):mpiLongitudes%indexEnd(i)), &
		nLat * mpiLongitudes%messageLength(i), MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode )
enddo

constInterpErr = abs( constInterp - constVal )
harmInterpError = abs( harmInterp - harmData )
harmLapError = abs( harmLap - harmLapExact )

!
!	interpolation to particles & fields
!
constScalarLap%N = sphere%particles%N
harmonicLapPSE%N = sphere%particles%N
harmonicLapError%N = sphere%particles%N
harmonicInterp%N = sphere%particles%N
harmonicInterpError%N = sphere%particles%N
constScalarInterp%N = sphere%particles%N 
streamFnScalar%N = sphere%particles%N

call SetFieldToZero( constScalarLap )
call SetFieldToZero( harmonicLapPSE )
call SetFieldToZero( streamFnScalar )

do i = mpiParticles%indexStart(procRank), mpiParticles%indexEnd(procRank)
	xi = PhysCoord(sphere%particles, i)
	
	constScalarInterp%scalar(i) = PSESphereInterpolateScalar( pseSetup, sphere, constScalar, xi )
	harmonicInterp%scalar(i) = PSESphereInterpolateScalar( pseSetup, sphere, harmonic, xi)
	
	do j = 1, sphere%particles%N
		if ( sphere%particles%isActive(j) ) then
			xj = PhysCoord(sphere%particles, j)
			pseKin = SphereDistance( xi, xj ) / pseSetup%eps
			lapKernel = bivariateLaplacianKernel8( pseKin ) / pseSetup%eps**2
			greensKernel = log( radius * radius - sum(xi*xj) ) / fourPi

			if ( j == i ) greensKernel = 0.0_kreal
			
			constScalarLap%scalar(i) = constScalarLap%scalar(i) + lapKernel * ( constScalar%scalar(j) - constScalar%scalar(i)) * &
				sphere%particles%area(j) / pseSetup%eps**2
			harmonicLapPSE%scalar(i) = harmonicLapPSE%scalar(i) + lapKernel * ( harmonic%scalar(j) - harmonic%scalar(i) ) * &
				sphere%particles%area(j) / pseSetup%eps**2
			
			streamFnScalar%scalar(i) = streamFnScalar%scalar(i) + greensKernel * harmonicLapExact%scalar(j) * &
				sphere%particles%area(j)
		endif
	enddo
enddo

do i = 0, numProcs - 1
	call MPI_BCAST( constScalarInterp%scalar( mpiParticles%indexStart(i): mpiParticles%indexEnd(i)), &
		mpiParticles%messageLength(i), MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)
	call MPI_BCAST( constScalarLap%scalar( mpiParticles%indexStart(i): mpiParticles%indexEnd(i)), &
		mpiParticles%messageLength(i), MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)
	call MPI_BCAST( harmonicInterp%scalar( mpiParticles%indexStart(i): mpiParticles%indexEnd(i)), &
		mpiParticles%messageLength(i), MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)
	call MPI_BCAST( harmonicLapPSE%scalar( mpiParticles%indexStart(i): mpiParticles%indexEnd(i)), &
		mpiParticles%messageLength(i), MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)
	call MPI_BCAST( streamFnScalar%scalar( mpiParticles%indexStart(i): mpiParticles%indexEnd(i)), &
		mpiParticles%messageLength(i), MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)		
enddo

harmonicInterpError%scalar = abs( harmonicInterp%scalar - harmonic%scalar )
harmonicLapError%scalar = abs( harmonicLapPSE%scalar - harmonicLapExact%scalar )

!
!	compute error norms
!
particlesLinfConst = maxval(abs( constScalarInterp%scalar - constVal )) / abs( constVal)
particlesLinfConstLap = maxval( abs( constScalarLap%scalar ))
particlesLinfHarm = maxval( harmonicInterpError%scalar ) / maxval(abs(harmonic%scalar))
particlesLinfHarmLap = maxval( harmonicLapError%scalar ) / maxval(abs(harmonicLapExact%scalar))

harmDenom = 0.0_kreal
harmLapDenom = 0.0_kreal
particlesL2Const = 0.0_kreal
particlesL2ConstLap = 0.0_kreal
particlesL2Harm = 0.0_kreal
particlesL2HarmLap = 0.0_kreal 
do k = mpiParticles%indexStart(procRank), mpiParticles%indexEnd(procRank)
	if ( sphere%particles%isActive(k) ) then
		harmDenom = harmDenom + harmonic%scalar(k)**2 * sphere%particles%area(k)
		harmLapDenom = harmLapDenom + harmonicLapExact%scalar(k)**2 * sphere%particles%area(k)
		
		particlesL2Const = particlesL2Const + ( constScalarInterp%scalar(k) - constVal )**2 * sphere%particles%area(k)
		particlesL2ConstLap = particlesL2ConstLap + constScalarLap%scalar(k)**2 * sphere%particles%area(k)
		
		particlesL2Harm = particlesL2Harm + harmonicInterpError%scalar(k)**2 * sphere%particles%area(k)
		particlesL2HarmLap = particlesL2HarmLap + harmonicLapError%scalar(k)**2 * sphere%particles%area(k)
	endif
enddo

call MPI_Reduce( [harmDenom, harmLapDenom, particlesL2Const, particlesL2ConstLap, particlesL2Harm, particlesL2HarmLap], & 
	[harmDenom, harmLapDenom, particlesL2Const, particlesL2ConstLap, particlesL2Harm, particlesL2HarmLap],&
	6, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, mpiErrCode )

particlesL2Const = particlesL2Const / ( 4.0_kreal * PI * radius * radius * constVal * constVal )
particlesL2Harm = particlesL2Harm / harmDenom
particlesL2HarmLap = particlesL2HarmLap / harmLapDenom

unifL2Const = 0.0_kreal
unifL2ConstLap = 0.0_kreal
unifL2Harm = 0.0_kreal
unifL2HarmLap = 0.0_kreal
do j = mpiLongitudes%indexStart(procRank), mpiLongitudes%indexEnd(procRank)
	do i = 1, nLat
		unifL2Const = unifL2Const + ( constInterp(i,j) - constVal )**2 * cos(lats(i))
		unifL2Harm = unifL2Harm + ( harmInterp(i,j) - harmData(i,j) )**2 * cos(lats(i))
		unifL2ConstLap = unifL2ConstLap + constLap(i,j) ** 2 * cos(lats(i))
		unifL2HarmLap = unifL2HarmLap + ( harmLap(i,j) - harmLapExact(i,j) )**2 * cos(lats(i))
	enddo
enddo

call MPI_Reduce( [unifL2Const, unifL2ConstLap, unifL2Harm, unifL2HarmLap], &
				 [unifL2Const, unifL2ConstLap, unifL2Harm, unifL2HarmLap], 4, MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
				 MPI_COMM_WORLD, mpiErrCode )

unifLinfConst = maxval( abs( constInterp - constVal )) / abs( constVal )
unifLinfConstLap = maxval( abs( constLap ) )
unifLinfHarm = maxval( abs( harmInterp - harmData )) / maxval(abs(harmData))
unifLinfHarmLap = maxval( abs( harmLap - harmLapExact)) / maxval(abs(harmLapExact))

unifL2Const = unifL2Const / ( 4.0_kreal * PI * radius * radius * constVal * constVal )
unifL2Harm = unifL2Harm / harmDenom
unifL2HarmLap = unifL2HarmLap / harmLapDenom 

call StartSection(exeLog, "Particles to Uniform Grid approximations")
call LogMessage(exeLog, TRACE_LOGGING_LEVEL, "constInterpLinf = ", unifLinfConst )
call LogMessage(exeLog, TRACE_LOGGING_LEVEL, "constInterpL2 = ", unifL2Const )
call LogMessage(exeLog, TRACE_LOGGING_LEVEL, "constLapLinf = ", unifLinfConstLap )
call LogMessage(exeLog, TRACE_LOGGING_LEVEL, "constLapL2 = ", unifL2ConstLap )
call LogMessage(exeLog, TRACE_LOGGING_LEVEL, "harmInterpLinf = ", unifLinfHarm )
call LogMessage(exeLog, TRACE_LOGGING_LEVEL, "harmInterpL2 = ", unifL2Harm )
call LogMessage(exeLog, TRACE_LOGGING_LEVEL, "harmLapLinf = ", unifLinfHarmLap )
call LogMessage(exeLog, TRACE_LOGGING_LEVEL, "harmLapL2 = ", unifL2HarmLap )
call EndSection(exeLog)

call StartSection(exeLog, "Particles to Particles approximations")
call LogMessage(exeLog, TRACE_LOGGING_LEVEL, "constInterpLinf = ", particlesLinfConst )
call LogMessage(exeLog, TRACE_LOGGING_LEVEL, "constInterpL2 = ", particlesL2Const )
call LogMessage(exeLog, TRACE_LOGGING_LEVEL, "constLapLinf = ", particlesLinfConstLap )
call LogMessage(exeLog, TRACE_LOGGING_LEVEL, "constLapL2 = ", particlesL2ConstLap )
call LogMessage(exeLog, TRACE_LOGGING_LEVEL, "harmInterpLinf = ", particlesLinfHarm )
call LogMessage(exeLog, TRACE_LOGGING_LEVEL, "harmInterpL2 = ", particlesL2Harm )
call LogMessage(exeLog, TRACE_LOGGING_LEVEL, "harmLapLinf = ", particlesLinfHarmLap )
call LogMessage(exeLog, TRACE_LOGGING_LEVEL, "harmLapL2 = ", particlesL2HarmLap )
call EndSection(exeLog)

harmonicInterpError%scalar = harmonicInterpError%scalar / maxval(abs(harmonic%scalar))
harmonicLapError%scalar = harmonicLapError%scalar / maxval(abs(harmonicLapExact%scalar))


!--------------------------------
!	finalize : write output, clean up
!--------------------------------
if ( procRank == 0 ) then
	!
	!	write output to vtk/paraview
	!
	open( unit=WRITE_UNIT_1, file=vtkFile, status='REPLACE', action='WRITE', iostat=writeStat)
		if ( writeStat /= 0 ) then
			call LogMessage(exeLog, ERROR_LOGGING_LEVEL, trim(logkey)//" OutputToVTK ERROR writing to file = ", trim(vtkFile))
		else	
			call WriteVTKPoints( sphere%particles, WRITE_UNIT_1)
			call WriteFacesToVTKPolygons( sphere%faces, WRITE_UNIT_1)
	
			call WriteVTKPointDataSectionHeader(WRITE_UNIT_1, sphere%particles%N)
		
			call WriteFieldToVTKPointData( constScalar, WRITE_UNIT_1)
			call WriteFieldToVTKPointData( constScalarInterp, WRITE_UNIT_1)
			call WriteFieldToVTKPointData( constScalarLap, WRITE_UNIT_1)
		
			call WriteFieldToVTKPointData( harmonic, WRITE_UNIT_1 )
			call WriteFieldToVTKPointData( harmonicInterp, WRITE_UNIT_1 )
			call WriteFieldToVTKPointData( harmonicInterpError, WRITE_UNIT_1 )
		
			call WriteFieldToVTKPointData( harmonicLapPSE, WRITE_UNIT_1 )	
			call WriteFieldToVTKPointData( harmonicLapExact, WRITE_UNIT_1 )
			call WriteFieldToVTKPointData( harmonicLapError, WRITE_UNIT_1 )
		
			call WriteFieldToVTKPointData( streamFnScalar, WRITE_UNIT_1 )
	
			call WriteFaceAreaToVTKCellData( sphere%faces, sphere%particles, WRITE_UNIT_1)
		endif
	close(WRITE_UNIT_1)
	!
	!	write output to matlab
	!
	open( unit=WRITE_UNIT_1, file=matlabFile, status='REPLACE', action='WRITE', iostat=writeStat)
		if ( writeStat /= 0 ) then
			call LogMessage(exeLog, ERROR_LOGGING_LEVEL, trim(logkey)//" OutputToMatlab ERROR writing to file = ", trim(matlabFile))
		else
			call WriteToMatlab(lons, WRITE_UNIT_1, "lons")
			call WriteToMatlab(lats, WRITE_UNIT_1, "lats")

			call WriteToMatlab( constData, WRITE_UNIT_1, "const")
			call WriteToMatlab( constLap, WRITE_UNIT_1, "constLap")
			call WriteToMatlab( constInterp, WRITE_UNIT_1, "constInterp")
		
			call WriteToMatlab( harmData, WRITE_UNIT_1, "harm54")		
			call WriteToMatlab( harmInterp, WRITE_UNIT_1, "harmInterp")
			call WriteToMatlab( harmLap, WRITE_UNIT_1, "harmLap")
			call WriteToMatlab( harmLapExact, WRITE_UNIT_1, "harmLapExact")
		
			call WriteToMatlab( stream, WRITE_UNIT_1, "streamFn")
		endif
	close(WRITE_UNIT_1)
endif
!
!----------------
! PROGRAM END
!----------------
!
timeEnd = MPI_WTIME()

write(logstring,'(A,F12.2,A)') "PROGRAM COMPLETE : elapsed time ", timeEnd - timeStart, " seconds."
call LogMessage(exelog, TRACE_LOGGING_LEVEL, trim(logkey)//" ", logString)

call EndSection(exeLog)

call Delete(harmonicLapError)
call Delete(harmonicLapExact)
call Delete(harmonicLapPSE)
call Delete(pseSetup)
call Delete(mpiParticles)
call Delete(mpiLongitudes)
call Delete(constScalar)
call Delete(harmonic)
call Delete(sphere)

call MPI_FINALIZE(mpiErrCode)

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
	ConstantScalarFn = constVal
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

pure function ExactHarmonicLaplacian( x, y, z )
	real(kreal) :: ExactHarmonicLaplacian
	real(kreal), intent(in) :: x, y, z
	ExactHarmonicLaplacian = -30.0_kreal * SphericalHarmonicFn(x,y,z)
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
				write(meshString, '(A,I1)') '_icosTri', initNest
				write(testTitle,'(A,I1,A)') "PSE Convergence Test, Icos. Tri., nest = ", initNest, " mesh size = "
			else
				write(meshString, '(2(A,I1))') '_icosTriAMR', initNest, 'to', maxNest
			endif
		elseif (meshSeed == CUBED_SPHERE_SEED ) then
			if ( initNest == maxNest ) then
				write(meshString, '(A,I1)') '_cubedSphere', initNest
				write(testTitle,'(A,I1,A)') "PSE Convergence Test, Cubed Sphere, nest = ", initNest, " mesh size = "
			else
				write(meshString, '(2(A,I1))') '_cubedSphereAMR', initNest, 'to', maxNest
			endif
		endif
	
		write(vtkFile,'(4A)') trim(outputDir), trim(outputRoot), trim(meshString), '.vtk'
		write(matlabFile,'(5A)') trim(outputDir), '/', trim(outputRoot), trim(meshString), '.m'
		
		bcastIntegers(1) = meshSeed
		bcastIntegers(2) = initNest
		bcastIntegers(3) = maxNest
		bcastIntegers(4) = amrLimit
		
		bcastReals(1) = radius
		bcastReals(2) = psePower
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
	psePower = bcastReals(2)
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