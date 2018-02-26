program CubicMeshTester

use NumberKindsModule
use LoggerModule
use PolyMesh2dModule
use ParticlesOOModule
use EdgesOOModule
use FacesOOModule
use STDIntVectorModule


implicit none

integer(kint) :: i, j, k

type(Logger) :: exeLog

integer(kint) :: initNest

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call New(exeLog, DEBUG_LOGGING_LEVEL)
call getInput()
call LogMessage(exeLog, TRACE_LOGGING_LEVEL, "CubicMeshTest: initNest = ", initNest)


if (testPass) call LogMessage(exeLog, TRACE_LOGGING_LEVEL, "CubicMeshTest result: ", "PASSED.")

contains

subroutine getInput()
    integer(kint) :: narg
    character(len=100) :: arg
    narg = IARGC()
    if (narg==0) then
        initNest = 1
    else 
        call GETARG(1, arg)
        read(arg, *) initNest
    endif
end subroutine

end program