program CubicGLLTest

#include "LpmConfig.h"

use NumberKindsModule
use UtilitiesModule
use CubicGLLModule
use LoggerModule
use SphereGeomModule
use PlaneGeomModule

implicit none

real(kreal), dimension(3,4) :: planeVerts, sphereVerts
real(kreal), dimension(3) :: elXy
real(kreal), dimension(2) :: rcrd, dvec

character(len=56) :: logkey = "CubicGLLTest:"
character(len=MAX_STRING_LENGTH) :: logstring
integer(kint) :: i
real(kreal) :: xi, eta
type(Logger) :: exeLog

testPass = .TRUE.

planeVerts = dzero
planeVerts(1:2,1) = [4.0_kreal, 3.0_kreal]
planeVerts(1:2,2) = [2.0_kreal, 2.0_kreal]
planeVerts(1:2,3) = [5.0_kreal, 1.0_kreal]
planeVerts(1:2,4) = [8.0_kreal, 3.0_kreal]

elXy = bilinearMap(planeVerts, -one, one)
if ( sqrt(sum( (elxy-planeVerts(:,2))*(elxy-planeVerts(:,1)))) > ZERO_TOL ) then
    write(logstring,'(2(A,G18.8),A)') "bilinearMap(-1,1) = (", elXy(1), ", ", elXy(1), "); correct solution is (4,3)"
    call LogMessage(exeLog, ERROR_LOGGING_LEVEL, trim(logkey)//" ", logstring)
    testPass = .FALSE.
endif

elXy = bilinearMap(planeVerts, -one, -one)
if ( sqrt(sum( (elxy-planeVerts(:,2))*(elxy-planeVerts(:,2)))) > ZERO_TOL ) then
    write(logstring,'(2(A,G18.8),A)') "bilinearMap(-1,-1) = (", elXy(1), ", ", elXy(1), "); correct solution is (2,2)"
    call LogMessage(exeLog, ERROR_LOGGING_LEVEL, trim(logkey)//" ", logstring)
    testPass = .FALSE.
endif

elXy = bilinearMap(planeVerts, one, -one)
if ( sqrt(sum( (elxy-planeVerts(:,3))*(elxy-planeVerts(:,3)))) > ZERO_TOL ) then
    write(logstring,'(2(A,G18.8),A)') "bilinearMap(1,-1) = (", elXy(1), ", ", elXy(1), "); correct solution is (5,1)"
    call LogMessage(exeLog, ERROR_LOGGING_LEVEL, trim(logkey)//" ", logstring)
    testPass = .FALSE.
endif

elXy = bilinearMap(planeVerts, one, one)
if ( sqrt(sum( (elxy-planeVerts(:,4))*(elxy-planeVerts(:,4)))) > ZERO_TOL ) then
    write(logstring,'(2(A,G18.8),A)') "bilinearMap(1,1) = (", elXy(1), ", ", elXy(1), "); correct solution is (8,3)"
    call LogMessage(exeLog, ERROR_LOGGING_LEVEL, trim(logkey)//" error: ", logstring)
    testPass = .FALSE.
endif

elxy = [19.0_kreal, 9.0_kreal, dzero] * 0.25_kreal
rcrd = invertBilinear(planeVerts, elXy)
if (sqrt(sum(rcrd*rcrd)) > ZERO_TOL) then
    write(logstring,'(A,2G15.8,A)') 'invertBilinear([19/4,9/4,0]) = (', rcrd, "); correct solution is (0,0)"
    call LogMessage(exeLog,ERROR_LOGGING_LEVEL,trim(logkey)//" error: ", logstring)
    testPass = .FALSE.
endif

call RANDOM_SEED()
do i=1, 100
    call RANDOM_NUMBER(xi)
    call RANDOM_NUMBER(eta)
    xi = 2.0_kreal * xi - one
    eta = 2.0_kreal * eta - one
    elXy = bilinearMap(planeVerts, xi, eta)
    rcrd = invertBilinear(planeVerts, elXy)
    dvec(1) = rcrd(1) - xi
    dvec(2) = rcrd(2) - eta
    if (sqrt(sum(dvec*dvec)) > ZERO_TOL) then
        write(logstring,'(2(A,G15.8),2(A,2G15.8),a)') 'invertBilinear([', xi, ',',eta,  ']) = (', rcrd, &
               '); correct solution is (', xi, eta, ')'
        call LogMessage(exeLog, ERROR_LOGGING_LEVEL, trim(logKey)// " error: ", logstring)
        testPass = .FALSE.
    endif
enddo
if (testPass) call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logkey)//" result: ", "PASSED")

end program
