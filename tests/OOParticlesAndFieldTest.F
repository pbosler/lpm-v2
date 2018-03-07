program OOParticlesAndFieldTest

use NumberKindsModule
use LoggerModule
use UtilitiesModule
use ParticlesOOModule
use FieldOOModule

implicit none

integer(kint) :: i, j
type(Logger) :: exeLog
type(Field) :: scalarField
type(Field) :: vectorField
type(Particles) :: aParticles

real(kreal), parameter :: dx = 0.1_kreal
real(kreal), parameter :: xmin = -5.0_kreal
real(kreal), parameter :: xmax = 5.0_kreal
real(kreal), parameter :: ymin = xmin
real(kreal), parameter :: ymax = xmax
integer(kint), parameter :: nn = 101

real(kreal), dimension(2) :: vec

character(len=56) :: logkey = "ParticlesAndFieldTest"
character(len=256) :: logstring

call New(exeLog, DEBUG_LOGGING_LEVEL)
write(logstring,'(A,I4,A,I4,A)') "part 1: build a particle set from an ", nn, " x ", nn, " uniform grid."
call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logkey)//" ", trim(logstring))
call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logkey)//" ", "part 2: define a scalar field on the particles")
call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logkey)//" ", "part 3: define a vector field on the particles")
call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logkey)//" ", "part 4: output particles and fields to matlab m file.")

call aParticles%init(nn*nn, PLANAR_GEOM)

do i=1,nn
    vec(1) = xmin + dx * (i-1)
    do j=1,nn
        vec(2) = ymin + dx *(j-1)
        call aParticles%insert(vec, vec)
    enddo
enddo

call aparticles%logstats(exeLog)

call scalarField%init(1, aParticles%N, "sineWave")
do i=1, aParticles%N
    call scalarField%insertScalar(sineWave([aParticles%x(i), aParticles%y(i)]))
enddo
call scalarField%logStats(exeLog)

call vectorField%init(2, aParticles%N, "gradSineWave")
do i=1, aParticles%N
    call vectorField%insertVector(sineWaveGrad([aParticles%x(i), aParticles%y(i)]))
enddo
call vectorField%logStats(exeLog)

open(unit=WRITE_UNIT_1, file='ooParticlesAndFieldTest.m', action='write', status='replace')
call aParticles%writeMatlab(WRITE_UNIT_1)
call scalarField%writeMatlab(WRITE_UNIT_1)
call vectorField%writeMatlab(WRITE_UNIT_1)
close(WRITE_UNIT_1)

if (testPass) call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logkey)//" result: ", "PASSED.")

contains

pure function sineWave(xy)
    real(kreal) :: sineWave
    real(kreal), dimension(2), intent(in) :: xy
    sineWave = sin(2.0_kreal*PI* xy(1)/5.0_kreal) * sin(2.0_kreal*PI*xy(2)/5.0_kreal)
end function

pure function sineWaveGrad(xy)
    real(kreal), dimension(2) :: sineWaveGrad
    real(kreal), dimension(2), intent(in) :: xy

    sinewaveGrad(1) = 2.0_kreal*PI/5.0_kreal * cos(2.0_kreal*PI*xy(1)/5.0_kreal) * sin(2.0_kreal*PI*xy(2)/5.0_kreal)
    sinewaveGrad(2) = 2.0_kreal*PI/5.0_kreal * sin(2.0_kreal*PI*xy(1)/5.0_kreal) * cos(2.0_kreal*PI*xy(2)/5.0_kreal)
end function


end program
