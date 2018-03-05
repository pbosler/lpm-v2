program OOPlaneMestTest

use NumberKindsModule
use LoggerModule
use PolyMeshOOModule
use FieldOOModule

implicit none

type(PolyMesh2d), pointer :: mesh=>null()
type(Field), pointer :: scalarField => null()
type(Field), pointer :: vectorField => null()
integer(kint) ::  initNest, i, j
character(len=56) :: mesh_type, vtkFile

type(Logger) :: exeLog
character(len=56) :: logkey = "OOPlaneMestTest"
character(len=MAX_STRING_LENGTH) :: logstring

call New(exeLog, DEBUG_LOGGING_LEVEL)
call getInput(initNest)

call LogMessage(exeLog, TRACE_LOGGING_LEVEL,trim(logkey)//" part 1: build each planar mesh to initNest = ", initNest)
call LogMessage(exeLog, TRACE_LOGGING_LEVEL,trim(logkey), " part 2: define scalar and vector fields on mesh.")
call LogMessage(exeLog, TRACE_LOGGING_LEVEL,trim(logkey), " part 3: demonstrate quadrature on mesh.")
call LogMessage(exeLog, TRACE_LOGGING_LEVEL,trim(logKey), " part 4: output to vtk xml file.")

!-------------------------------------
!   Test Planar mesh, linear triangles
!-------------------------------------
mesh_type = "planar_tri"
mesh_type = "planar_quad"
mesh_type = "planar_cubic_quad"
write(vtkFile,'(A,I1,A)') trim(mesh_type), initNest, '.vtp'
call StartSection(exeLog,'Planar Tri Test.')

allocate(mesh)
allocate(scalarField)
allocate(vectorField)
call mesh%init(mesh_type, initNest, initNest, 0, 1.0_kreal)
call mesh%logStats(exeLog)

call scalarField%init(1, mesh%particles%N, "sineWave")
do i=1, mesh%particles%N
    call scalarField%insertScalar(sineWave([mesh%particles%x(i), mesh%particles%y(i)]))
enddo    

call vectorField%init(2, mesh%particles%N, "gradSineWave")
do i=1, mesh%particles%n
    call vectorField%insertVector(sineWaveGrad([mesh%particles%x(i), mesh%particles%y(i)]))
enddo

call LogMessage(exeLog, TRACE_LOGGING_LEVEL, "sineWave integral = ", mesh%integrateScalar(scalarField))

open(unit=WRITE_UNIT_1, file=vtkFile, action='write', status='replace')
call mesh%writeVTKSerialStartXML(WRITE_UNIT_1)
write(WRITE_UNIT_1,'(A)') '     <PointData Scalars="sineWave" Vectors="gradSineWave">'
call scalarField%writeVtkPointDataXML(WRITE_UNIT_1)
call vectorField%writeVtkPointDataXML(WRITE_UNIT_1)
write(WRITE_UNIT_1,'(A)') '     </PointData>'
call mesh%writeVTKSerialEndXML(WRITE_UNIT_1)
close(WRITE_UNIT_1)

write(vtkFile,'(A,I1,A)') trim(mesh_type), initNest, '.m'
open(unit=WRITE_UNIT_1, file=vtkFile, action='write', status='replace')
call mesh%writeMatlab(WRITE_UNIT_1)
call scalarField%writeMatlab(WRITE_UNIT_1)
call vectorField%writeMatlab(WRITE_UNIT_1)
close(WRITE_UNIT_1)

deallocate(vectorField)
deallocate(scalarField)
deallocate(mesh)
nullify(mesh)
call EndSection(exeLog)
!-------------------------------------
!   Test Planar mesh, linear quads
!-------------------------------------
mesh_type = "planar_quad"


!-------------------------------------
!   Test Planar mesh, cubic quads
!-------------------------------------
mesh_type = "planar_cubic_quad"

testPass = .FALSE.


if (testPass) call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logkey)//" result:", " PASSED.")

contains

subroutine getInput(initNest)
    integer(kint), intent(out) :: initNest
    !
    integer(kint) :: argc
    character(len=MAX_STRING_LENGTH) :: argv
    argc = IARGC()
    if (argc == 0) then
        initNest = 0
    else
        call GETARG(1,argv)
        read(argv,*) initNest
    endif
end subroutine

pure function sineWave(xy)
    real(kreal) :: sineWave
    real(kreal), dimension(2), intent(in) :: xy
    sineWave = sin(2.0_kreal*PI* xy(1)*1.0_kreal) * sin(2.0_kreal*PI*xy(2)*1.0_kreal)
end function

pure function sineWaveGrad(xy)
    real(kreal), dimension(2) :: sineWaveGrad
    real(kreal), dimension(2), intent(in) :: xy

    sinewaveGrad(1) = 2.0_kreal*PI*1.0_kreal * cos(2.0_kreal*PI*xy(1)*1.0_kreal) * sin(2.0_kreal*PI*xy(2)*1.0_kreal)
    sinewaveGrad(2) = 2.0_kreal*PI*1.0_kreal * sin(2.0_kreal*PI*xy(1)*1.0_kreal) * cos(2.0_kreal*PI*xy(2)*1.0_kreal)
end function

end program 