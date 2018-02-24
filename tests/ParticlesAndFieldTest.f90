program ParticlesAndFieldTester

use NumberKindsModule
use LoggerModule
use ParticlesModule
use FieldModule

implicit none

integer(kint) :: i, j
real(kreal), parameter :: dx = 0.1;
real(kreal), parameter :: xmin = -5.0_kreal
real(kreal), parameter :: xmax = 5.0_kreal
real(kreal), parameter :: ymin = xmin
real(kreal), parameter :: ymax = xmax
integer(kint), parameter :: nn = 101

type(Logger) :: exeLog
type(Particles) :: aParticles
type(Field) :: scalarField
type(Field) :: vectorField

real(kreal) :: vec(2)
print *, "ParticlesAndFieldTest : 1. build a particle set from a ", nn," x ", nn, " uniform planar mesh."
print *, "ParticlesAndFieldTest : 2. define a scalar field on that particle set."
print *, "ParticlesAndFieldTest : 3. define a vector field on that particle set."
print *, "ParticlesAndFieldTest : 4. test Matlab output."

call New(exeLog, DEBUG_LOGGING_LEVEL )

open(unit=WRITE_UNIT_1, file="particlesFieldTest.m", STATUS='REPLACE', action="WRITE")

call New(aParticles, nn * nn, PLANAR_GEOM )

do i = 1, nn
	vec(1) = xmin + dx * (i-1)
	do j = 1, nn
		vec(2) = ymin + dx * (j-1)
		call InsertParticle(aParticles, vec, vec )
	enddo
enddo

call LogStats(aParticles, exeLog )
call WriteParticlesToMatlab( aParticles, WRITE_UNIT_1 )

call New(scalarField, 1, aParticles%N, "sineWaves" )

do i = 1, aParticles%N
	call InsertScalarToField( scalarField, SineWaves([aParticles%x(i), aParticles%y(i)]) )
enddo

call LogStats(scalarField, exeLog )
call WriteFieldToMatlab(scalarField, WRITE_UNIT_1)

call New(vectorField, 2, aParticles%N, "sineWavesGrad")

do i = 1, aParticles%N
	call InsertVectorToField( vectorField, SineWavesGrad([aParticles%x(i), aParticles%y(i)]) )
enddo

call LogStats(vectorField, exeLog)
call WriteFieldToMatlab(vectorField, WRITE_UNIT_1)

write(WRITE_UNIT_1,*) "figure(1);clf; scatter(x,y,8,scalarField_sineWaves);"
write(WRITE_UNIT_1,*) "sineWavesGrad = vectorField_sineWavesGrad';"
write(WRITE_UNIT_1,*) "figure(2);clf; quiver(x,y,sineWavesGrad(1,:),sineWavesGrad(2,:));"

call Delete(vectorField)
call Delete(scalarField)
call Delete(aParticles)
close(WRITE_UNIT_1)
call Delete(exeLog)

contains 

function SineWaves( xy )
	real(kreal) :: SineWaves
	real(kreal), intent(in) :: xy(2)
	SineWaves = sin( 2.0_kreal * PI * xy(1) / 5.0_kreal) * sin(2.0_kreal * PI * xy(2) / 5.0_kreal )
end function

function SineWavesGrad( xy ) 
	real(kreal) :: SineWavesGrad(2)
	real(kreal), intent(in) :: xy(2)
	SineWavesGrad(1) = 2.0_kreal * PI / 5.0_kreal * cos(2.0_kreal * PI * xy(1) / 5.0_kreal) * &
		sin(2.0_kreal * PI * xy(2) / 5.0_kreal)
	SineWavesGrad(2) = 2.0_kreal * PI / 5.0_kreal * sin(2.0_kreal * PI * xy(1) / 5.0_kreal) * &
	    cos(2.0_kreal * PI * xy(2) / 5.0_kreal)
end function 

end program