program OOPlaneMestTest

use NumberKindsModule
use LoggerModule
use OutputWriterModule
use UtilitiesModule
use PolyMeshOOModule
use FieldOOModule

implicit none

type(PolyMesh2d), pointer :: mesh=>null()
type(Field), pointer :: scalarField => null()
type(Field), pointer :: vectorField => null()
integer(kint) ::  initNest, i, j, doOutput
character(len=56) :: mesh_type, vtkFile

real(kreal), parameter :: intCosBump = 16.0_kreal / PI**2
real(kreal), parameter :: triIntExact = 1.4741613930426135683_kreal
real(kreal), parameter :: dx = 2.0_kreal / 50.0_kreal


type(Logger) :: exeLog
character(len=56) :: logkey = "OOPlaneMestTest"
character(len=MAX_STRING_LENGTH) :: logstring

call New(exeLog, DEBUG_LOGGING_LEVEL)
call getInput(initNest, doOutput)

call LogMessage(exeLog, TRACE_LOGGING_LEVEL,trim(logkey)//" part 1: build each planar mesh to initNest = ", initNest)
call LogMessage(exeLog, TRACE_LOGGING_LEVEL,trim(logkey), " part 2: define scalar and vector fields on mesh.")
call LogMessage(exeLog, TRACE_LOGGING_LEVEL,trim(logkey), " part 3: demonstrate quadrature on mesh.")
call LogMessage(exeLog, TRACE_LOGGING_LEVEL,trim(logkey), " part 4: demonstrate interpolation on a cubic mesh.")
call LogMessage(exeLog, TRACE_LOGGING_LEVEL,trim(logKey), " part 5: output to vtk xml file.")

!-------------------------------------
!   Test Planar mesh, linear triangles
!-------------------------------------
mesh_type = "planar_tri"
write(vtkFile,'(A,I1,A)') trim(mesh_type), initNest, '.vtp'
call StartSection(exeLog,'Planar Tri Test.')
call runTest(mesh, mesh_type, initNest, vtkFile)
call EndSection(exeLog)

!-------------------------------------
!   Test Planar mesh, linear quads
!-------------------------------------
mesh_type = "planar_quad"
write(vtkFile,'(A,I1,A)') trim(mesh_type), initNest, '.vtp'
call StartSection(exeLog,'Planar Quad Test.')
call runTest(mesh, mesh_type, initNest, vtkFile)
call EndSection(exeLog)


!-------------------------------------
!   Test Planar mesh, cubic quads
!-------------------------------------


mesh_type = "planar_cubic_quad"
write(vtkFile,'(A,I1,A)') trim(mesh_type), initNest, '.vtp'
call StartSection(exeLog,'Planar Cubic Quad Test.')
call runTest(mesh, mesh_type, initNest, vtkFile)
call EndSection(exeLog)

if (testPass) call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logkey)//" result:", " PASSED.")

contains

subroutine runTest(mesh_ptr, mesh_type, initNest, oname)
    type(PolyMesh2d), pointer, intent(inout) :: mesh_ptr
    character(len=*), intent(in) :: mesh_type
    integer(kint), intent(in) :: initNest
    character(len=*), intent(in) :: oname
    !
    real(kreal) :: computedIntegral
    real(kreal), dimension(3) :: interpTargetPt
    integer(kint) :: i, j
    real(kreal) :: interpErr
    real(kreal), dimension(51) :: unifX, unifY
    real(kreal), dimension(51,51) :: interp, exact

    allocate(mesh_ptr)
    allocate(scalarField)
    allocate(vectorField)
    call mesh_ptr%init(mesh_type, initNest, initNest, 0, 1.0_kreal)
    call mesh_ptr%logStats(exeLog)

    call scalarField%init(1, mesh_ptr%particles%N, "scalar")
    do i=1, mesh_ptr%particles%N
        call scalarField%insertScalar(cosineBump([mesh_ptr%particles%x(i), mesh_ptr%particles%y(i)]))
    enddo

    call vectorField%init(2, mesh_ptr%particles%N, "scalarGradient")
    do i=1, mesh_ptr%particles%n
        call vectorField%insertVector(bumpGrad([mesh_ptr%particles%x(i), mesh_ptr%particles%y(i)]))
    enddo

    computedIntegral = mesh_ptr%integrateScalar(scalarField)
    interpErr = dzero
    if (mesh_ptr%faceKind == QUAD_CUBIC_PANEL) then
        do j=1,51
            unifx(j) = -one + dx*(j-1)
        enddo
        do i=1, 51
            unify(i) = -one + dx*(i-1)
        enddo
        do j=1,51
            do i=1,51
                exact(i,j) = cosineBump([unifX(j), unifY(i)])
            enddo
        enddo
        interpTargetPt = dzero
        do j=1,51
            do i=1,51
                interpTargetPt(1) = unifx(j)
                interpTargetPt(2) = unify(i)
                interp(i,j) = mesh_ptr%interpolateScalar(interpTargetPt, scalarField)
            enddo
        enddo
        interpErr = maxval(abs(exact-interp))
        call LogMessage(exeLog, TRACE_LOGGING_LEVEL, "interp error = ", interpErr)
    endif

    call LogMessage(exeLog, TRACE_LOGGING_LEVEL, "scalar integral = ", computedIntegral)
    if (mesh_ptr%faceKind == TRI_PANEL) then
        call LogMessage(exeLog, TRACE_LOGGING_LEVEL, "mesh size = ", mesh_ptr%edges%avgLength(mesh_ptr%particles))
        call LogMessage(exeLog, TRACE_LOGGING_LEVEL, "abs(integralErr) = ", abs(computedIntegral-triIntExact))
    elseif (mesh_ptr%faceKind == QUAD_PANEL .or. mesh_ptr%faceKind==QUAD_CUBIC_PANEL) then
        call LogMessage(exeLog, TRACE_LOGGING_LEVEL, "mesh size = ", mesh_ptr%edges%avgLength(mesh_ptr%particles))
        call LogMessage(exeLog, TRACE_LOGGING_LEVEL, "abs(integralErr) = ", abs(computedIntegral-intCosBump))
    endif

    if (doOutput>0) then
        open(unit=WRITE_UNIT_1, file=vtkFile, action='write', status='replace')
        call mesh_ptr%writeVTKSerialStartXML(WRITE_UNIT_1)
        write(WRITE_UNIT_1,'(A)') '     <PointData Scalars="sineWave" Vectors="gradSineWave">'
        call scalarField%writeVtkPointDataXML(WRITE_UNIT_1)
        call vectorField%writeVtkPointDataXML(WRITE_UNIT_1)
        write(WRITE_UNIT_1,'(A)') '     </PointData>'
        call mesh_ptr%writeVTKSerialEndXML(WRITE_UNIT_1)
        close(WRITE_UNIT_1)
        write(vtkFile,'(A,I1,A)') trim(mesh_type), initNest, '.m'

        if (mesh_ptr%faceKind == QUAD_CUBIC_PANEL) then
            open(unit=WRITE_UNIT_1, file=vtkFile, action='write', status='replace')
            call mesh%writeMatlab(WRITE_UNIT_1)
            call scalarField%writeMatlab(WRITE_UNIT_1)
            call vectorField%writeMatlab(WRITE_UNIT_1)
            call writeToMatlab(unifX, WRITE_UNIT_1, "unifx")
            call writeToMatlab(unify, WRITE_UNIT_1, "unify")
            call writeToMatlab(exact, WRITE_UNIT_1, "exact")
            call writeToMatlab(interp, WRITE_UNIT_1, "interp")
            close(WRITE_UNIT_1)
        endif
    endif



    deallocate(vectorField)
    deallocate(scalarField)
    deallocate(mesh)
    nullify(mesh)
    nullify(vectorField)
    nullify(scalarField)
end subroutine

subroutine getInput(initNest, doOutput)
    integer(kint), intent(out) :: initNest, doOutput
    !
    integer(kint) :: argc
    character(len=MAX_STRING_LENGTH) :: argv
    argc = IARGC()
    if (argc == 0) then
        initNest = 0
        doOutput = 1
    elseif (argc==1) then
        call GETARG(1,argv)
        read(argv,*) initNest
        doOutput = 0
    else
        call GETARG(1,argv)
        read(argv,*) initNest
        call GETARG(2,argv)
        read(argv,*) doOutput
    endif
end subroutine

pure function cosineBump(xy)
    real(kreal) :: cosineBump
    real(kreal), intent(in) :: xy(2)
    cosineBump = cos(0.5_kreal * PI * xy(1)) * cos(0.5_kreal * PI * xy(2))
end function

pure function bumpGrad(xy)
    real(kreal) :: bumpGrad(2)
    real(kreal), intent(in) :: xy(2)
    bumpGrad(1) = -0.5_kreal * PI * sin(0.5_kreal*PI*xy(1)) * sin(0.5_kreal*PI*xy(2))
    bumpGrad(2) = -0.5_kreal *PI * cos(0.5_kreal*PI*xy(1))*sin(0.5_kreal*PI*xy(2))
end function

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
