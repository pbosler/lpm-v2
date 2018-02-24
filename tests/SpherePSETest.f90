program TestPSESphere

use NumberKindsModule
use OutputWriterModule
use LoggerModule
use ParticlesModule
use FacesModule
use PolyMesh2dModule
use FieldModule
use MPISetupModule
use PSEDirectSumModule
use SphereGeomModule

implicit none

include 'mpif.h'

! mesh variables
type(PolyMesh2d) :: sphere
integer(kint) :: initNest
integer(kint) :: maxNest
integer(kint) :: meshSeed
integer(kint) :: amrLimit
real(kreal) :: radius

! field tests
type(Field) :: constantScalar
type(Field) :: estGradConst
type(Field) :: estLapConst
type(Field) :: sphHarm
type(Field) :: estGradSphHarm
type(Field) :: estLapSphHarm
type(Field) :: exactGradHarm
type(Field) :: lapSphHarmError
type(Field) :: gradHarmError
real(kreal) :: maxGradMag
real(kreal) :: maxLapMag

! pse interpolation/approximation
type(PSE) :: pseSetup

! uniform grid
real(kreal) :: lons(360), lats(181)
real(kreal) :: interp(181,360)
real(kreal) :: harmInterp(181,360)
real(kreal) :: exactHarm(181,360)
real(kreal) :: xInterp(3)
integer(kint), parameter :: nLat = 181
integer(kint), parameter :: nLon = 360

! environment / machine variables
type(MPISetup) :: particlesMPI
type(MPISetup) :: unifGridMPI
integer(kint) :: mpiErrCode, narg
character(len=100) :: arg
type(Logger) :: exeLog
integer(kint), parameter :: logLevel = DEBUG_LOGGING_LEVEL
character(len=MAX_STRING_LENGTH) :: logString
character(len=15) :: logKey = "spherePSE"
real(kreal) :: programStart, programEnd

! output variables
character(len=100) :: matlabFilename
character(len=100) :: paraviewFilename

! general
integer(kint) :: i, j
real(kreal) :: xi(3), vecA(3), vecB(3)
!
!----------------
! PROGRAM START
!----------------
!
call MPI_INIT(mpiErrCode)
call MPI_COMM_SIZE(MPI_COMM_WORLD, numProcs, mpiErrCode)
call MPI_COMM_RANK(MPI_COMM_WORLD, procRank, mpiErrCode)

!print *, "Hello from proc ", procRank, " of ", numProcs

programStart = MPI_WTIME()

call InitLogger(exeLog, procRank)

if ( procRank == 0 ) then
	narg = IARGC()
	if ( narg /= 2 ) then
		stop "usage : spherePSETest.exe meshType initNest"
	else
		call GETARG(1, arg)
		read(arg, *) meshSeed
		if ( meshSeed == 3 ) then
			meshSeed = ICOS_TRI_SPHERE_SEED
		elseif (meshSeed == 4) then
			meshSeed = CUBED_SPHERE_SEED
		else
			stop "invalid arg1: meshType must be 3 or 4"
		endif
		call GETARG(2, arg)
		read(arg, *) initNest
	endif
endif
call MPI_BCAST(meshSeed, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpiErrCode)
call MPI_BCAST(initNest, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpiErrCode)

if ( meshSeed == ICOS_TRI_SPHERE_SEED ) then
	write(logString,'(A,I1,A)') "building icos tri mesh to initNest = ", initNest, "..."
	write(paraviewFilename,'(A,I1,A)') "icosTriPSE_", initNest, ".vtk"
	write(matlabFilename,'(A,I1,A)') "icosTriPSE_", initNest, ".m"
else
	write(logString,'(A,I1,A)') "building cubed sphere mesh to initNest = ", initNest, "..."
	write(paraviewFilename,'(A,I1,A)') "cubedSpherePSE_", initNest, ".vtk"
	write(matlabFilename,'(A,I1,A)') "cubedSpherePSE_", initNest, ".m"
endif
call LogMessage(exeLog, TRACE_LOGGING_LEVEL, logKey, logString )


!
!----------------
! build mesh
!----------------
!
maxNest = initNest
amrLimit = 0
radius = 1.0_kreal
call New( sphere, meshSeed, initNest, maxNest, amrLimit, radius )
call LogStats( sphere, exeLog )

call New( particlesMPI, sphere%particles%N, numProcs)

do j = 1, 360
	lons(j) = (j-1) * DEG_2_RAD
enddo
do i = 1, 181
	lats(i) = -PI/2.0_kreal + (i-1) * DEG_2_RAD
enddo
do j = 1, 360
	do i = 1, 181
		exactHarm(i,j) = 3.0_kreal * sqrt(35.0_kreal) * cos( 4.0_kreal * lons(j)) * sin(lats(i)) * &
						 (-1.0_kreal + sin(lats(i))*sin(lats(i)))**2
	enddo
enddo

call New( unifGridMPI, 360, numProcs )

!
!----------------
! set up source field data
!----------------
!
call New(constantScalar, 1, sphere%particles%N, "const", "n/a")

call New( sphHarm, 1, sphere%particles%N, "sphHarm54", "n/a")
call New( exactGradHarm, 3, sphere%particles%N, "gradSphHarm", "n/a")

do i = 1, sphere%particles%N
	xi = [sphere%particles%x(i), sphere%particles%y(i), sphere%particles%z(i)]
	call InsertScalarToField( constantScalar, 1.0_kreal )
	call InsertScalarToField( sphHarm, SphereHarmonic54( xi ))
	call InsertVectorToField( exactGradHarm, HarmGradient(xi) ) 
enddo



!
!----------------
! do PSE approximations
!----------------
!
call New(pseSetup, sphere, 1.5_kreal )

do j = unifGridMPI%indexStart(procRank), unifGridMPI%indexEnd(procRank)
	do i = 1, 181
		xInterp = [ radius * cos(lons(j)) * cos(lats(i)), radius * sin(lons(j)) * cos(lats(i)), radius * sin(lats(i))]
		interp(i,j) = PSESphereInterpolateScalar( pseSetup, sphere, constantScalar, xInterp )
		harmInterp(i,j) = PSESphereInterpolateScalar( pseSetup, sphere, sphHarm, xInterp)
	enddo
enddo
do i = 0, numProcs - 1
	call MPI_BCAST( interp(:, unifGridMPI%indexStart(i):unifGridMPI%indexEnd(i)), 181 * unifGridMPI%messageLength(i), &
				    MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)
	call MPI_BCAST( harmInterp(:, unifGridMPI%indexStart(i):unifGridMPI%indexEnd(i)), 181 * unifGridMPI%messageLength(i),&
					MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)
enddo

call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logkey), " ... PSE interpolation complete.")
call MPI_BARRIER(MPI_COMM_WORLD, mpiErrCode)

call New( estGradConst, 3, sphere%particles%N, "estGrad-const", "n/a")
call PSESphereGradientAtParticles( pseSetup, sphere, constantScalar, estGradConst, particlesMPI )

call New( estLapConst, 1, sphere%particles%N, "estLap-const", "n/a")
call PSESphereLaplacianAtParticles( pseSetup, sphere, constantScalar, estLapConst, particlesMPI)

call New( estGradSphHarm, 3, sphere%particles%N, "estGrad-sphHarm", "n/a")
call PSESphereGradientAtParticles(pseSetup, sphere, sphHarm, estGradSphHarm, particlesMPI)

call New( estLapSphHarm, 1, sphere%particles%N, "estLap-sphHarm", "n/a")
call PSESphereLaplacianAtParticles(pseSetup, sphere, sphHarm, estLapSphHarm, particlesMPI)

call New( lapSphHarmError, 1, sphere%particles%N, "errorLapSphHarm", "n/a")
call New( gradHarmError, 1, sphere%particles%N, "errorGradSphHarm", "n/a")
maxGradMag = MaxMagnitude( exactGradHarm )
maxLapMag = 30.0_kreal * maxval( abs(sphHarm%scalar) )
do i = 1, sphere%particles%N
	call InsertScalarToField( lapSphHarmError, abs( estLapSphHarm%scalar(i) + 30.0_kreal * &
				SphereHarmonic54( [sphere%particles%x(i), sphere%particles%y(i), sphere%particles%z(i)] ) )/maxLapMag )
	vecA = [estGradSphHarm%xComp(i), estGradSphHarm%yComp(i), estGradSphHarm%zComp(i)]
	vecB = [exactGradHarm%xComp(i), exactGradHarm%yComp(i), exactGradHarm%zComp(i) ]
	call InsertScalarToField( gradHarmError, sqrt(sum( (vecB- vecA)*(vecB-vecA)))/maxGradMag )			
enddo

call LogMessage(exeLog,TRACE_LOGGING_LEVEL,trim(logkey)//" interp err = ", &
		maxval(abs(harmInterp-exactHarm))/maxval(abs(exactHarm)))
call LogMessage(exeLog,TRACE_LOGGING_LEVEL,trim(logkey)//" gradient err = ", maxval(gradHarmError%scalar))
call LogMessage(exeLog,TRACE_LOGGING_LEVEL,trim(logkey)//" lap err = ", maxval(abs(lapSphHarmError%scalar)) )

!
!----------------
! write output
!----------------
!
if ( procRank == 0 ) then
	open(unit=WRITE_UNIT_1, file=paraviewFilename, action='WRITE',status='REPLACE')
		call WriteVTKPoints( sphere%particles, WRITE_UNIT_1 )
		call WriteFacesToVTKPolygons( sphere%faces, WRITE_UNIT_1)
		! write all point data before first cell data
		call WriteVTKPointDataSectionHeader(WRITE_UNIT_1, sphere%particles%N)
		call WriteVTKLagCoords( sphere%particles, WRITE_UNIT_1)
		call WriteFieldToVTKPointData( constantScalar, WRITE_UNIT_1)
		call WriteFieldToVTKPointData( estGradConst, WRITE_UNIT_1)
		call WriteFieldToVTKPointData( estLapConst, WRITE_UNIT_1)
		call WriteFieldToVTKPointData( sphHarm, WRITE_UNIT_1)
		call WriteFieldToVTKPointData( estGradSphHarm, WRITE_UNIT_1)
		call WriteFieldToVTKPointData( estLapSphHarm, WRITE_UNIT_1)
		call WriteFieldToVTKPointData( lapSphHarmError, WRITE_UNIT_1)
		call WriteFieldToVTKPointData( gradHarmError, WRITE_UNIT_1)
		call WriteFieldToVTKPointData( exactGradHarm, WRITE_UNIT_1)
		! now only write cell data
		call WriteFaceAreaToVTKCellData(sphere%faces, sphere%particles, WRITE_UNIT_1)	
		call WriteFieldToVTKCellData( constantScalar, WRITE_UNIT_1, sphere%faces)
	close(WRITE_UNIT_1)
	if ( numProcs == 1 ) then
		open(unit=WRITE_UNIT_2, file=matlabFilename, action='WRITE',status='REPLACE')
			call WriteToMatlab( lons, WRITE_UNIT_2, "lons")
			call WriteToMatlab( lats, WRITE_UNIT_2, "lats")
			call WriteToMatlab( interp, WRITE_UNIT_2, "constScalarInterp")
		close(WRITE_UNIT_2)
	endif
endif

if ( procRank == 1 ) then
	open(unit=WRITE_UNIT_2, file=matlabFilename, action='WRITE',status='REPLACE')
		call WriteToMatlab( lons, WRITE_UNIT_2, "lons")
		call WriteToMatlab( lats, WRITE_UNIT_2, "lats")
		call WriteToMatlab( interp, WRITE_UNIT_2, "constScalarInterp")
		call WriteToMatlab( harmInterp, WRITE_UNIT_2, "sphHarmInterp")
		call WriteToMatlab( exactHarm, WRITE_UNIT_2, "exactHarm")
	close(WRITE_UNIT_2)
endif

!
!----------------
! PROGRAM END
!----------------
!
programEnd = MPI_WTIME()

write(logstring,'(A,F12.2,A)') "PROGRAM COMPLETE : elapsed time ", programEnd - programStart, " seconds."
call LogMessage(exelog, TRACE_LOGGING_LEVEL, trim(logkey)//" ", logString)

call Delete(gradHarmError)
call Delete(lapSphHarmError)
call Delete(estLapSphHarm)
call Delete(estGradSphHarm)
call Delete(sphHarm)
call Delete(estGradConst)
call Delete(constantScalar)
call Delete(pseSetup)
call Delete(unifGridMPI)
call Delete(particlesMPI)
call Delete(sphere)
call Delete(exeLog)
call MPI_Finalize(mpiErrCode)

contains

pure function SphereHarmonic54(xyz)
	real(kreal) :: SphereHarmonic54
	real(kreal), intent(in) :: xyz(3)
	SphereHarmonic54 = 3.0_kreal * sqrt(35.0_kreal) * cos( 4.0_kreal * Longitude(xyz)) * sin(Latitude(xyz)) * &
					   (-1.0_kreal + sin(Latitude(xyz)) * sin(Latitude(xyz)))**2
end function

pure function HarmGradient( xyz ) 
	real(kreal) :: HarmGradient(3)
	real(kreal), intent(in) :: xyz(3)
	!
	real(kreal) :: lon, lat, u, v
	
	lat = Latitude(xyz)
	lon = Longitude(xyz)
	
	u = -12.0_kreal*sqrt(35.0_kreal)*sin(4.0_kreal * lon)*(-1.0_kreal + sin(lat)*sin(lat))**2 * tan(lat)
	v =  12.0_kreal*sqrt(35.0_kreal)*cos(4.0_kreal * lon)*cos(lat)*sin(lat)*sin(lat)*(-1.0_kreal + sin(lat)*sin(lat)) + &
		 3.0_kreal*sqrt(35.0_kreal) * cos(4.0_kreal * lon)*cos(lat)*(-1.0_kreal + sin(lat)*sin(lat))**2
	
	HarmGradient(1) =-u * sin(lon) - v * sin(lat)*cos(lon)
	HarmGradient(2) = u * cos(lon) - v * sin(lat)*sin(lon)
	HarmGradient(3) = v * cos(lat)
end function

subroutine InitLogger(aLog,rank)
! Initialize a logger for this processor
	type(Logger), intent(out) :: aLog
	integer(kint), intent(in) :: rank
	write(logKey,'(A,A,I0.2,A)') trim(logKey),'_',rank,' : '
	if ( rank == 0 ) then
		call New(aLog,logLevel)
	else
		call New(aLog,ERROR_LOGGING_LEVEL)
	endif
end subroutine

end program