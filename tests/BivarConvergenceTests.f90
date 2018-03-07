program PlanarInterpConvergence

use NumberKindsModule
use UtilitiesModule
use LoggerModule
use PolyMesh2dModule
use ParticlesModule
use EdgesModule
use FacesModule
use STDIntVectorModule
use FieldModule
use BIVARModule

implicit none

type(Logger) :: exeLog

type(PolyMesh2d) :: triMesh
integer(kint), parameter :: meshSeed = TRI_HEX_SEED
integer(kint) :: initNest
integer(kint) :: maxNest
integer(kint) :: amrLimit
real(kreal) :: ampFactor

type(Field) :: scalar
type(Field) :: estGrad, exactGrad, gradError
type(Field) :: est2ndPartials, exact2ndPartials, partialsError
type(Field) :: estLap, exactLap, lapError

real(kreal), parameter :: b = 3.0_kreal
real(kreal), parameter :: xc = 0.0_kreal, yc = 0.0_kreal
real(kreal), parameter :: maxAbsLap = 36.0_kreal
integer(kint), parameter :: nestStart = 5, nestEnd = 5, testSize = nestEnd - nestStart + 1
integer(kint) :: nestCtr
real(kreal) :: estGradError(testSize), estLapError(testSize), interpError(testSize), meshSize(testSize)
real(kreal) :: testStart, testEnd, maxGradMag

integer(kint) :: i, j

integer(kint) :: nTri, nL
integer(kint), allocatable :: ipt(:), ipl(:), iwl(:), iwp(:)
real(kreal), allocatable :: wk(:), partials(:)

integer(kint), parameter :: nn = 501

real(kreal), parameter :: dx = 8.0_kreal / real( nn - 1, kreal)
real(kreal), parameter :: xmin = -4.0_kreal
real(kreal), parameter :: xmax = 4.0_kreal
real(kreal), parameter :: ymin = xmin
real(kreal), parameter :: ymax = xmax

real(kreal) :: x(nn), y(nn)
real(kreal) :: interpScalar(nn,nn), exactScalar(nn,nn)
real(kreal) :: xVecA(3), xVecB(3)
integer(kint) :: inTri(nn,nn)

character(len=MAX_STRING_LENGTH) :: logstring
character(len=MAX_STRING_LENGTH) :: filename
!
!----------------
! PROGRAM START
!----------------
!
call New(exeLog, DEBUG_LOGGING_LEVEL)


nestCtr = 1

do initNest = nestStart, nestEnd

	write(logstring,'(2(A,I3),A)') "test ", nestCtr, ", of ", nestEnd-nestStart +1, "..."

	call LogMessage(exeLog, TRACE_LOGGING_LEVEL,"Interpolation Convergence : ", logString)
	call cpu_time(testStart)
	!
	! build a mesh
	!
	maxNest = initNest
	amrLimit = 0
	ampFactor = 3.0_kreal
	call New(triMesh, meshSeed, initNest, maxNest, amrLimit, ampFactor)
	! define a scalar field on the mesh
	call New(scalar, 1, triMesh%particles%N, "gaussScalar", "n/a")
	call New(estGrad,2, triMesh%particles%N, "estGradient", "n/a")
	call New(exactGrad,2,triMesh%particles%N,"exactGradient","n/a")
	call New(gradError,1,triMesh%particles%N, "gradError", "n/a")
	call New(est2ndPartials,3,triMesh%particles%N, "estPartials","n/a")
	call New(exact2ndPartials,3,triMesh%particles%N,"exactPartials","n/a")
	call New(partialsError,3,triMesh%particles%N, "partialsError","n/a")
	call New(estLap, 1, triMesh%particles%N, "estLaplacian","n/a")
	call New(exactLap,1, triMesh%particles%N,"exactLaplacian", "n/a")
	call New(lapError,1, triMesh%particles%N,"lapError","n/a")

	do i = 1, triMesh%particles%N
		call InsertScalarToField( scalar, Gaussian( [triMesh%particles%x(i), triMesh%particles%y(i)], b) )
		call InsertVectorToField(exactGrad, GaussGrad( [triMesh%particles%x(i), triMesh%particles%y(i)], b))
		call InsertVectorToField(exact2ndPartials, Gauss2ndDerivs([triMesh%particles%x(i), triMesh%particles%y(i)], b))
		call InsertScalarToField(exactLap, GaussLap([triMesh%particles%x(i), triMesh%particles%y(i)],b) )
	enddo

	!
	! define uniform mesh to serve as interpolation target
	!
	do i = 1, nn
		x(i) = xmin + dx * (i-1)
		y(i) = ymin + dx * (i-1)
	enddo

	allocate(ipt(6*triMesh%particles%N - 15))  ! triangle vertices : triangle I has vertices at
											   !  	indices ipt(3*I-2), ipt(3*I-1), and ipt(3*I), for I = 1,...,nTri.
	allocate(ipl(6*triMesh%particles%N))	   ! border edges and triangles : border edge I has
											   ! 	endpoint 1, endpoint 2, and triangle index stored at
											   !	ipl(3*I-2), ipl(3*I-1), and ipl(3*I), for I = 1,...,nL.
	allocate(iwl(18*triMesh%particles%N))	! workspace only
	allocate(iwp(triMesh%particles%N))		! workspace only
	allocate(wk(8*triMesh%particles%N))		! workspace only
	!
	!	triangulate particles
	!
	call idtang( triMesh%particles%N, triMesh%particles%x, triMesh%particles%y, &
				 nTri, ipt, nL, ipl, iwl, iwp, wk)
	!
	!	locate triangles that contain destination points
	!
	inTri = 0
	do j = 1, nn
		do i = 1, nn
			call idlctn( triMesh%particles%N, triMesh%particles%x, triMesh%particles%y, nTri, ipt, nL, ipl, &
						 x(j), y(i), inTri(i,j), iwl, wk )
		enddo
	enddo
	!
	!	estimate partial derivatives at particles
	!
	allocate(partials(5*triMesh%particles%N)) ! partial derivatives : partials of of particle I stored at
											  !		partials(5*I-4:5*I)
	call idpdrv( triMesh%particles%N, triMesh%particles%x, triMesh%particles%y, scalar%scalar, nTri, ipt, partials, wk)
	do i = 1, triMesh%particles%N
		call InsertVectorToField( estGrad, partials(5*i-4:5*i-3) )
		call InsertVectorToField( est2ndPartials, partials(5*i-2:5*i))
		call InsertScalarToField( estLap, partials(5*i-2) + partials(5*i))
	enddo

	do i = 1, triMesh%particles%N
		partials(5*i-4:5*i-3) = GaussGrad( [triMesh%particles%x(i), triMesh%particles%y(i)], b)
		partials(5*i-2:5*i) = Gauss2ndDerivs( [triMesh%particles%x(i), triMesh%particles%y(i)], b)
	enddo

	!
	! 	interpolate the scalar
	!
	interpScalar = 0.0_kreal
	do j = 1, nn
		do i = 1, nn
			exactScalar(i,j) = Gaussian( [x(j),y(i)], b)
			if ( pointIsOutsideMesh(triMesh, [ x(j), y(i), 0.0_kreal ]) ) then
				interpScalar(i,j) = 0.0_kreal
			else
				call idptip( triMesh%particles%N, triMesh%particles%x, triMesh%particles%y, scalar%scalar, &
						 nTri, ipt, nL, ipl, partials, inTri(i,j), x(j), y(i), interpScalar(i,j))
			endif
		enddo
	enddo


	!
	! calculate derivative error at particles
	!
	do i = 1, triMesh%particles%N
		xVecA = 0.0_kreal
		xVecB = 0.0_kreal
		xVecA(1:2) = [estGrad%xComp(i), estGrad%yComp(i)]
		xVecB(1:2) = GaussGrad( [triMesh%particles%x(i), triMesh%particles%y(i)], b)
		call InsertScalarToField( gradError, sqrt(sum( xVecA - xVecB)))
		call InsertScalarToField( lapError, estLap%scalar(i) - exactLap%scalar(i))
		xVecA = [est2ndPartials%xComp(i), est2ndPartials%yComp(i), est2ndPartials%zComp(i)]
		xVecB = Gauss2ndDerivs( [triMesh%particles%x(i), triMesh%particles%y(i)], b)
		call InsertVectorToField( partialsError, xVecA - xVecB )
	enddo

	maxGradMag = MaxMagnitude(exactGrad)
	meshSize(nestCtr) = MaxEdgeLength(triMesh%edges, triMesh%particles)
	estGradError(nestCtr) = maxval(gradError%scalar)/maxGradMag
	estLapError(nestCtr) = maxval(abs(lapError%scalar))/maxAbsLap
	interpError(nestCtr) = maxval(abs(interpScalar-exactScalar))

	write(6,'(4A24)') "dx", "gradErr-particles", "lapErr-particles", "interp error"
	write(6,'(4F24.10)') meshSize(nestCtr), estGradError(nestCtr), estLapError(nestCtr), interpError(nestCtr)

	write(filename,'(A,I1,A)') 'BivarTestTriMesh', initNest, '.m'
	open(unit=WRITE_UNIT_1,file=filename,status='REPLACE',action='WRITE')
		call WriteParticlesToMatlab(triMesh%particles, WRITE_UNIT_1)
		call WriteFieldToMatlab(scalar, WRITE_UNIT_1)
		call WriteFieldToMatlab(estGrad, WRITE_UNIT_1)
		call WriteFieldToMatlab(exactGrad, WRITE_UNIT_1)
		call WriteFieldToMatlab(estLap, WRITE_UNIT_1)
		call WriteFieldToMatlab(exactLap, WRITE_UNIT_1)
		call WriteFieldToMatlab(gradError, WRITE_UNIT_1)
		call WriteFieldToMatlab(lapError, WRITE_UNIT_1)

		write(WRITE_UNIT_1,'(A)',advance='NO') "xi = ["
		do i = 1, nn - 1
			write(WRITE_UNIT_1,'(F18.12,A)', advance='NO') x(i), ", "
		enddo
		write(WRITE_UNIT_1,'(F18.12,A)') x(nn), "];"
		write(WRITE_UNIT_1, '(A)') "yi = xi; "

		write(WRITE_UNIT_1,'(A)',advance='NO') "interp = ["
		do i = 1, nn - 1
			do j = 1, nn - 1
				write(WRITE_UNIT_1,'(F18.12,A)',advance='NO') interpScalar(i,j), ", "
			enddo
			write(WRITE_UNIT_1,'(F18.12,A)') interpScalar(i,nn), "; ..."
		enddo
		do j = 1, nn - 1
			write(WRITE_UNIT_1,'(F18.12,A)',advance='NO') interpScalar(nn,j), ", "
		enddo
		write(WRITE_UNIT_1,'(F18.12,A)') interpScalar(nn,nn), "];"

	close(WRITE_UNIT_1)

	call cpu_time(testEnd)
	write(6,'(A,I8,A,F12.2,A)') "nParticles = ", triMesh%particles%N, ": elapsed time = ", testEnd-testStart, " seconds."

	!
	! reset for next iteration
	!
	deallocate(partials)
	deallocate(wk)
	deallocate(iwl)
	deallocate(iwp)
	deallocate(ipl)
	deallocate(ipt)
	call Delete(lapError)
	call Delete(exactLap)
	call Delete(estLap)
	call Delete(partialsError)
	call Delete(exact2ndPartials)
	call Delete(est2ndPartials)
	call Delete(gradError)
	call Delete(exactGrad)
	call Delete(estGrad)
	call Delete(scalar)
	call Delete(triMesh)

	call LogMessage(exeLog,TRACE_LOGGING_LEVEL, "test complete for initNest = ", initNest)
	nestCtr = nestCtr + 1
enddo
!
!write(6,'(4A24)') "dx", "gradErr-particles", "lapErr-particles", "interp error"
!do i = 1, testSize
!	write(6,'(4F24.10)') meshSize(i), estGradError(i), estLapError(i), interpError(i)
!enddo
!
!if ((meshSize(1) /= 0.375) .OR. estGradError(1) > 0.3953 .OR. interpError(1) > 0.0003134) then
!    testPass = .FALSE.
!    call LogMessage(exeLog, ERROR_LOGGING_LEVEL, "test result: ", "regression failed at index 1")
!endif
!
!if ((meshSize(2) /= 0.1875) .OR. estGradError(2) > 0.2091 .OR. interpError(2) > 0.0003134) then
!    testPass = .FALSE.
!    call LogMessage(exeLog, ERROR_LOGGING_LEVEL, "test result: ", "regression failed at index 2")
!endif

if ((meshSize(1) /= 0.09375) .OR. estGradError(1) > 0.10728 .OR. interpError(1) > 0.0003134) then
    testPass = .FALSE.
    call LogMessage(exeLog, ERROR_LOGGING_LEVEL, "test result: ", "regression failed at index 3")
endif

if (testPass) call LogMessage(exeLog, TRACE_LOGGING_LEVEL, "test result: ", "PASSED.")


contains

function Gaussian( xy, b )
	real(kreal) :: Gaussian
	real(kreal), intent(in) :: xy(2)
	real(kreal), intent(in) :: b

	Gaussian = exp( - b * b * ( (xy(1)-xc)*(xy(1)-xc) + (xy(2)-yc)*(xy(2)-yc)))
end function

function GaussGrad(xy, b)
	real(kreal) :: GaussGrad(2)
	real(kreal), intent(in) :: xy(2)
	real(kreal), intent(in) :: b
	GaussGrad(1) = xy(1)-xc
	GaussGrad(2) = xy(2)-yc
	GaussGrad = -2.0_kreal * GaussGrad * b * b * exp( - b * b * ( (xy(1)-xc)*(xy(1)-xc) + (xy(2)-yc)*(xy(2)-yc) ))
end function

function GaussLap(xy, b)
	real(kreal) :: GaussLap
	real(kreal), intent(in) :: xy(2)
	real(kreal), intent(in) :: b
	GaussLap = ( b*b * ( (xy(1)-xc)*(xy(1)-xc) + (xy(2)-yc)*(xy(2)-yc) ) - 1.0_kreal ) * 4.0_kreal * b * b * &
			exp( - b * b * ( (xy(1)-xc)*(xy(1)-xc) + (xy(2)-yc)*(xy(2)-yc) ))
end function


function Gauss2ndDerivs( xy, b)
	real(kreal) :: Gauss2ndDerivs(3)
	real(kreal), intent(in) :: xy(2)
	real(kreal), intent(in) :: b
	Gauss2ndDerivs(1) = 2.0_kreal * b * b *(2.0_kreal * b*b *(xy(1)-xc)*(xy(1)-xc) - 1.0_kreal)
	Gauss2ndDerivs(2) = 4.0_kreal * b**4 * (xy(1)-xc)*(xy(2)-yc)
	Gauss2ndDerivs(3) = 2.0_kreal * b * b *(2.0_kreal * b*b * (xy(2)-yc)*(xy(2)-yc) - 1.0_kreal)
	Gauss2ndDerivs = Gauss2ndDerivs * exp( - b * b * ( (xy(1)-xc)*(xy(1)-xc) + (xy(2)-yc)*(xy(2)-yc) ))
end function

end program
