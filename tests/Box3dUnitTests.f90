program Box3dTest

use NumberKindsModule
use LoggerModule
use Box3dModule
use SphereGeomModule
use OutputWriterModule
use UtilitiesModule

#include "LpmConfig.h"

implicit none

type(Box3d) :: unitBox
type(Box3d), dimension(:), allocatable :: kids
type(Box3d) :: box1
type(Box3d) :: box2

integer(kint), parameter :: logLevel = DEBUG_LOGGING_LEVEL
character(len=28) :: logKey = "Box3dTest"

real(kreal), dimension(3), parameter :: queryPt= one / 2.0_kreal
integer(kint) :: i, boxid

type(Logger) :: exeLog
character(len=256) :: logString

call InitLogger(exeLog, procRank)

call logMessage(exeLog, TRACE_LOGGING_LEVEL, "Box3d: ", "Unit Test Start")

call unitBox%init(-one, one, -one, one, -one, one)
call unitBox%logStats(exeLog)

if (unitBox%containsPoint(queryPt)) then
    call logMessage(exeLog, TRACE_LOGGING_LEVEL, "unit box contains ", queryPt)
else
    call logMessage(exeLog, ERROR_LOGGING_LEVEL, "unit box does not contain ", queryPt)
endif

allocate(kids(8))

kids = unitBox%bisectAll()
do i=1,8
    call kids(i)%logStats(exeLog)
    if (kids(i)%containsPoint(queryPt)) then
        boxid = i
        call logMessage(exeLog, TRACE_LOGGING_LEVEL, "(0.5, 0.5, 0.5) found in box ", boxid)
    endif
enddo

if (boxid /= 8) then
    call logMessage(exeLog, ERROR_LOGGING_LEVEL, "(0.5, 0.5, 0.5) not found ", " in box 8.")
endif

deallocate(kids)

if (testPass) then
    call logMessage(exeLog, TRACE_LOGGING_LEVEL, "result: ", "PASSED")
else
    call logMessage(exeLog, ERROR_LOGGING_LEVEL, "result: ", "FAIL")
endif

contains

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
