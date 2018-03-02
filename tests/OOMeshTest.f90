program OOMeshTest

use NumberKindsModule
use LoggerModule
use PolyMeshOOModule

implicit none

type(PolyMesh2d), pointer :: mesh => null()
integer(kint) :: initNest
character(len=56) :: mesh_type

type(Logger) :: exeLog


call New(exeLog, DEBUG_LOGGING_LEVEL)
call getInput(initNest)

!-------------------------------------
!   Test Planar mesh, linear triangles
!-------------------------------------
mesh_type = "planar_tri"
call StartSection(exeLog, "TEST START", mesh_type)
call runTest(mesh, mesh_type, initNest)
call EndSection(exeLog)

!-------------------------------------
!   Test Planar mesh, linear quads
!-------------------------------------

mesh_type = "planar_quad"
call StartSection(exeLog, "TEST START", mesh_type)
call runTest(mesh, mesh_type, initNest)
call EndSection(exeLog)

!-------------------------------------
!   Test Planar mesh, cubic quads
!-------------------------------------
!mesh_type = "planar_cubic_quad"
!call StartSection(exeLog, "TEST START", mesh_type)
!call runTest(mesh, mesh_type, initNest)
!call EndSection(exeLog)

!-------------------------------------
!   Test Spherical mesh, linear triangles
!-------------------------------------
mesh_type = "icos_tri_sphere"
call StartSection(exeLog, "TEST START", mesh_type)
call runTest(mesh, mesh_type, initNest)
call EndSection(exeLog)

!-------------------------------------
!   Test Spherical mesh, linear quads
!-------------------------------------
mesh_type = "cubed_sphere"
call StartSection(exeLog, "TEST START", mesh_type)
call runTest(mesh, mesh_type, initNest)
call EndSection(exeLog)


contains

subroutine runTest(ptr, meshType, initNest)
    type(PolyMesh2d), pointer :: ptr
    character(len=*), intent(in) :: meshType
    integer(kint), intent(in) :: initNest
    !
    integer(kint), parameter :: amrLimit = 0
    real(kreal), parameter :: amp = 1.0_kreal
    
    allocate(ptr)
    call ptr%init(meshType, initNest, initNest, amrLimit, amp)
    
    call ptr%logStats(exeLog)
    
    deallocate(ptr)
    nullify(ptr)
end subroutine

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