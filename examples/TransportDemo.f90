program transportDemo

use NumberKindsModule
use UtilitiesModule
use LoggerModule
use MPIReplicatedDataModule
use MeshedParticlesModule
use RK4Module

implicit none

include 'mpif.h'

type(Logger) :: exeLog
character(MAX_STRING_LENGTH) :: logString
integer(kint) :: mpiErrCode

character(len=56) :: mesh_type
integer(kint) :: initNest
integer(kint) :: maxNest
integer(kint) :: amrLimit
real(kreal), parameter :: domainRadius=2.0_kreal
type(MeshedParticles) :: plane
type(MeshedParticles), pointer :: tempPlane => null()
integer(kint) :: i
real(kreal), dimension(2) :: xy
logical(klog) :: doOutput = .FALSE.

type(MPIReplicatedData) :: mpiParticles
type(MPIReplicatedData) :: mpiFaces
type(TransportStepper) :: time
real(kreal) :: dt, tfinal
integer(kint) :: nTimesteps, timeI, frameCounter

real(kreal) :: l2Err, linfErr, l2denom, ev, err
real(kreal) :: timeStart, timeEnd

character(len=128) :: fname

call MPI_INIT(mpiErrCode)
call MPI_COMM_SIZE(MPI_COMM_WORLD, numProcs, mpiErrCode)
call MPI_COMM_RANK(MPI_COMM_WORLD, procRank, mpiErrCode)

call New(exeLog, DEBUG_LOGGING_LEVEL)

timeStart = MPI_WTIME()

call getInput(initNest, doOutput)
maxNest = initNest
amrLimit = 0


mesh_type='planar_cubic_quad'
call plane%init(mesh_type, initNest, maxNest, amrLimit, domainRadius)
!print *, "returned from init"
if (procRank==0) call plane%logStats(exeLog)

plane%density%N = plane%particles%N
call plane%density%setToConst(one)
call plane%set2dVelocityWithFunction(rigidRotation, dzero)
call plane%addTracers(1)
plane%tracers(1)%name = "threeTracers"
do i=1, plane%particles%N
    xy = [plane%particles%x(i), plane%particles%y(i)]
    call plane%tracers(1)%insertScalar(ThreeTracers(xy))
enddo

if (procRank == 0 .and. doOutput) then
    write(fname, '(A,I1,A,I0.4,A)') "transport_"//trim(mesh_type), initNest, "_", 0, '.vtp'
    open(unit=WRITE_UNIT_1, file=fname, action='write', status='replace')
    call plane%writeVTK(WRITE_UNIT_1)
    close(WRITE_UNIT_1)
endif

call mpiParticles%init(plane%particles%N, numProcs)
call mpiFaces%init(plane%faces%N, numProcs)
call time%init(plane, mpiParticles, mpiFaces, .FALSE.)
tfinal = one
nTimesteps = 100
dt = 0.01_kreal

frameCounter = 1
do timeI = 1, nTimesteps
    if (mod(timeI,25) == 0) then ! remesh
        !call LogMessage(exeLog, TRACE_LOGGING_LEVEL, "commencing remesh, t = ", plane%t)
        allocate(tempPlane)
        !print *, "allocated temp plane"
        
        call tempPlane%init(mesh_type, initNest, maxNest, amrLimit, domainRadius)
        tempPlane%density%n = tempPlane%particles%N
        call tempPlane%density%setToConst(one)
        call tempPlane%set2dVelocityWithFunction(rigidRotation, plane%t)
        call tempPlane%setLagCoordFromOther(plane)
        tempPlane%t = plane%t
        call tempPlane%addTracers(1)
        tempPlane%tracers(1)%name = "threeTracers"
        do i=1, tempPlane%particles%N
            xy = [tempPlane%particles%x0(i), tempPlane%particles%y0(i)]
            call tempPlane%tracers(1)%insertScalar(ThreeTracers(xy))
        enddo

        !print *, "data on temp plane set."
        
        !if (procrank==0) call tempPlane%logStats(exeLog)
        
        
        if (procRank==0) call LogMessage(exeLog, TRACE_LOGGING_LEVEL, "remesh complete, t = ", plane%t)
        
        call plane%copy(tempPlane)
        deallocate(tempPlane)
        nullify(tempPlane)
        
        if (procRank == 0 .and. doOutput) then
            write(fname, '(A,I1,A,I0.4,A)') "transport_"//trim(mesh_type), initNest, "_", frameCounter, '.vtp'
            open(unit=WRITE_UNIT_1, file=fname, action='write', status='replace')
            call plane%writeVTK(WRITE_UNIT_1)
            close(WRITE_UNIT_1)
            frameCounter = frameCounter + 1
        endif
    endif

    call time%stepTransport2d(plane, mpiParticles, mpiFaces, dt, rigidRotation)
    
    if (mod(timeI, 25) == 0 .and. procRank == 0) call LogMessage(exeLog, TRACE_LOGGING_LEVEL, "time = ", plane%t)

    if (procRank == 0 .and. doOutput) then
        write(fname, '(A,I1,A,I0.4,A)') "transport_"//trim(mesh_type), initNest, "_", frameCounter, '.vtp'
        open(unit=WRITE_UNIT_1, file=fname, action='write', status='replace')
        call plane%writeVTK(WRITE_UNIT_1)
        close(WRITE_UNIT_1)
        frameCounter = frameCounter + 1
    endif
enddo

l2denom = dzero
l2Err = dzero
linfErr = dzero
do i=1, plane%particles%N
    ! compute tracer error, output to console
    ev = ThreeTracers([plane%particles%x(i), plane%particles%y(i)])
    err = abs(ev-plane%tracers(1)%comp1(i))
    if ( err > linfErr) linfErr = err
    l2denom = l2denom + ev*ev*plane%particles%weight(i)
    l2Err = l2Err + err*err*plane%particles%weight(i)    
enddo
l2Err = sqrt(l2Err/l2denom)

timeEnd = MPI_WTIME()

if (procRank==0) then
    write(logstring,'(3(A,G15.8))') "mesh size = ", plane%meshSize(), ', lInfErr = ', linfErr, ', l2err = ', l2Err
    call LogMessage(exeLog, TRACE_LOGGING_LEVEL, "", logString)
    write(logstring,'(A,G15.8,A)') ' elapsed time ', (timeEnd - timeStart) / 60.0_kreal, '  minutes.'
    call LogMessage(exeLog,TRACE_LOGGING_LEVEL, "PROGRAM COMPLETE: ", logString)
endif




call MPI_FINALIZE(mpiErrCode)

contains

function ThreeTracers(xy)
	real(kreal) :: ThreeTracers
	real(kreal), intent(in) :: xy(2)
	!
	real(kreal) :: r, cosArg
	real(kreal), dimension(2) :: gcent

    gcent = [0.25_kreal, 0.5_kreal]
! cosine hill
!	r = sqrt( (xy(1) - 0.25_kreal)*(xy(1) - 0.25_kreal) + (xy(2) - 0.5_kreal)*(xy(2)-0.5_kreal))
!	if ( r >= 0.15_kreal ) then
!		cosArg = 0.15_kreal
!	else
!		cosArg = r
!	endif
!	cosArg = cosArg / 0.15_kreal
!	ThreeTracers = 0.25_kreal * (1.0_kreal + cos(PI * cosArg))
!
! slotted cylinder
!	r = sqrt( (xy(1) - 0.5_kreal)*(xy(1)-0.5_kreal) + (xy(2) - 0.75_kreal)*(xy(2) - 0.75_kreal))
!	if ( r < 0.15_kreal) ThreeTracers = ThreeTracers + 1.0_kreal
!	if ( xy(1) < 0.525_kreal .AND. xy(1) > 0.475_kreal ) then
!		if ( xy(2) > 0.6_kreal .AND. xy(2) < 0.85_kreal ) then
!			ThreeTracers = ThreeTracers - 1.0_kreal
!		endif
!	endif
!
! cone
!	r = 0.15_kreal - sqrt( (xy(1) -0.5_kreal)*(xy(1) - 0.5_kreal) + (xy(2) - 0.25_kreal)*(xy(2) - 0.25_kreal))
!	if ( r < 0.0_kreal ) r = 0.0_kreal
!	ThreeTracers = ThreeTracers + r
    ThreeTracers = exp(-10.0_kreal*sum((xy-gcent)*(xy-gcent)))
end function

subroutine getInput(initNest, doOutput)
    integer(kint), intent(out) :: initNest
    logical(klog), intent(out) :: doOutput
    !
    integer(kint) :: narg
    character(len=128) :: arg
    integer(kint) :: iread
    narg = IARGC()
    if (narg == 0) then
        initNest = 0
    else
        call GETARG(1,arg)
        read(arg,*) initNest
        if (narg > 1) then
            call GETARG(2,arg)
            read(arg,*) iread
            doOutput = (iread>0)
        endif
    endif
end subroutine

function rigidRotation(x, y, t)
    real(kreal), dimension(2) :: rigidRotation
    real(kreal), intent(in) :: x, y, t
    rigidRotation(1) = -y *2.0_kreal * PI
    rigidRotation(2) = x * 2.0_kreal * PI
end function

end program
