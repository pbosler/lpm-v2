program OOMeshTest

use NumberKindsModule
use UtilitiesModule
use LoggerModule
use PolyMeshOOModule

implicit none

type(PolyMesh2d), pointer :: mesh => null()
integer(kint) :: initNest
character(len=56) :: mesh_type, mfile

type(Logger) :: exeLog

testPass = .TRUE.

call New(exeLog, DEBUG_LOGGING_LEVEL)
call getInput(initNest)

!-------------------------------------
!   Test Planar mesh, linear triangles
!-------------------------------------
mesh_type = "planar_tri"
write(mfile,'(A,I1,A)') trim(mesh_type), initNest, '.m'
call StartSection(exeLog, "TEST START", mesh_type)
call runTest(mesh, mesh_type, initNest, mfile)
call EndSection(exeLog)
!
!!-------------------------------------
!!   Test Planar mesh, linear quads
!!-------------------------------------
!!
mesh_type = "planar_quad"
write(mfile,'(A,I1,A)') trim(mesh_type), initNest, '.m'
call StartSection(exeLog, "TEST START", mesh_type)
call runTest(mesh, mesh_type, initNest, mfile)
call EndSection(exeLog)
!!
!!-------------------------------------
!!   Test Planar mesh, cubic quads
!!-------------------------------------
mesh_type = "planar_cubic_quad"
write(mfile,'(A,I1,A)') trim(mesh_type), initNest, '.m'
call StartSection(exeLog, "TEST START", mesh_type)
call runTest(mesh, mesh_type, initNest, mfile)
call EndSection(exeLog)
!
!!-------------------------------------
!!   Test Spherical mesh, linear triangles
!!-------------------------------------
mesh_type = "icos_tri_sphere"
write(mfile,'(A,I1,A)') trim(mesh_type), initNest, '.m'
call StartSection(exeLog, "TEST START", mesh_type)
call runTest(mesh, mesh_type, initNest, mfile)
call EndSection(exeLog)
!!
!!!-------------------------------------
!!!   Test Spherical mesh, linear quads
!!!-------------------------------------
mesh_type = "cubed_sphere"
write(mfile,'(A,I1,A)') trim(mesh_type), initNest, '.m'
call StartSection(exeLog, "TEST START", mesh_type)
call runTest(mesh, mesh_type, initNest,mfile)
call EndSection(exeLog)


if (testPass) call LogMessage(exeLog, TRACE_LOGGING_LEVEL, "Test Result: ", "PASSED")

contains

subroutine runTest(ptr, meshType, initNest, oname)
    type(PolyMesh2d), pointer :: ptr
    character(len=*), intent(in) :: meshType
    integer(kint), intent(in) :: initNest
    character(len=*), intent(in) :: oname
    !
    integer(kint), parameter :: amrLimit = 0
    real(kreal), parameter :: amp = 1.0_kreal
    integer(kint) :: i

    allocate(ptr)
    call ptr%init(meshType, initNest, initNest, amrLimit, amp)
    call ptr%logStats(exeLog)

    open(unit=WRITE_UNIT_1, file=oname, status='replace', action='write')
    call ptr%writeMatlab(WRITE_UNIT_1)
    close(WRITE_UNIT_1)

    if (areaTest(ptr)) then
        print *, "AreaTest: pass"
    else
        call LogMessage(exeLog, ERROR_LOGGING_LEVEL, trim(mesh_type)//" areaTest ", "FAIL")
        call LogMessage(exeLog, ERROR_LOGGING_LEVEL, trim(mesh_type)//" faceArea ", surfAreaFaces(ptr))
        call LogMessage(exeLog, ERROR_LOGGING_LEVEL, trim(mesh_type)//" particleArea ", surfAreaParticles(ptr))
    endif
    deallocate(ptr)
    nullify(ptr)
end subroutine

function areaTest(ptr)
    logical(klog) :: areaTest
    type(PolyMesh2d), pointer :: ptr
    areaTest = (abs((surfAreaFaces(ptr) - surfAreaParticles(ptr))) < ZERO_TOL)
end function

function surfAreaFaces(ptr)
    real(kreal) :: surfAreaFaces
    type(PolyMesh2d), pointer :: ptr
    !
    integer(kint) :: i

    surfAreaFaces = dzero
    do i=1,ptr%faces%n
        !if (.not. ptr%faces%hasChildren(i)) then
            surfAreaFaces = surfAreaFaces + ptr%faces%area(i)
        !endif
    enddo
end function

function surfAreaParticles(ptr)
    real(kreal) :: surfAreaParticles
    type(PolyMesh2d), pointer :: ptr
    !
    integer(kint) :: i, j, nv, nc

    nv = size(ptr%faces%vertices,1)
    nc = size(ptr%faces%centerParticles,1)

    surfAreaParticles = dzero
    do i = 1, ptr%faces%n
        if (.not. ptr%faces%hasChildren(i)) then
            do j=1,nv
                surfAreaParticles = surfAreaParticles + ptr%particles%weight(ptr%faces%vertices(j,i))
            enddo
            do j=1,nc
                surfAreaParticles = surfAreaParticles + ptr%particles%weight(ptr%faces%centerParticles(j,i))
            enddo
        endif
    enddo
end function

subroutine getInput(initNest)
    integer(kint), intent(out) :: initNest
    !
    integer(kint) :: narg
    character(len=128) :: arg
    narg = IARGC()
    if (narg == 0) then
        initNest = 0
    else
        call GETARG(1,arg)
        read(arg,*) initNest
    endif
end subroutine

end program
