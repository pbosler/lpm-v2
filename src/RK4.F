module RK4Module
!>@{
use NumberKindsModule
use UtilitiesModule
use LoggerModule
use ParticlesOOModule
use EdgesOOModule
use FieldOOModule
use PolyMeshOOModule
use MeshedParticlesModule
use MPIReplicatedDataModule

implicit none

include 'mpif.h'

private

public  Timestepper, TransportStepper

type Timestepper
    real(kreal), dimension(:), allocatable :: xIn !< input to each RK stage
    real(kreal), dimension(:), allocatable :: yIn !< input to each RK stage
    real(kreal), dimension(:), allocatable :: zIn !< input to each RK stage
    real(kreal), dimension(:), allocatable :: weightIn !< input to each RK stage
    real(kreal), dimension(:), allocatable :: rhoIn !< input to each RK stage
    real(kreal), dimension(:), allocatable :: areaIn !< input to each RK stage
    
    real(kreal), dimension(:), allocatable :: xStage1
    real(kreal), dimension(:), allocatable :: xStage2
    real(kreal), dimension(:), allocatable :: xStage3
    real(kreal), dimension(:), allocatable :: xStage4
    
    real(kreal), dimension(:), allocatable :: yStage1
    real(kreal), dimension(:), allocatable :: yStage2
    real(kreal), dimension(:), allocatable :: yStage3
    real(kreal), dimension(:), allocatable :: yStage4
    
    real(kreal), dimension(:), allocatable :: zStage1
    real(kreal), dimension(:), allocatable :: zStage2
    real(kreal), dimension(:), allocatable :: zStage3
    real(kreal), dimension(:), allocatable :: zStage4
    
    real(kreal), dimension(:), allocatable :: weightStage1
    real(kreal), dimension(:), allocatable :: weightStage2
    real(kreal), dimension(:), allocatable :: weightStage3
    real(kreal), dimension(:), allocatable :: weightStage4
    
    real(kreal), dimension(:), allocatable :: rhoStage1
    real(kreal), dimension(:), allocatable :: rhoStage2
    real(kreal), dimension(:), allocatable :: rhoStage3
    real(kreal), dimension(:), allocatable :: rhoStage4
    
    real(kreal), dimension(:), allocatable :: areaStage1
    real(kreal), dimension(:), allocatable :: areaStage2
    real(kreal), dimension(:), allocatable :: areaStage3
    real(kreal), dimension(:), allocatable :: areaStage4
    
    contains
        procedure :: init
        final :: deleteRK4
end type

type, extends(TimeStepper) :: TransportStepper
!    procedure(vectorFnOf2DSpaceAndTime), pointer, nopass :: velFn2d => null()
!    procedure(vectorFnOf3DSpaceAndTime), pointer, nopass :: velFn3d => null()
    contains
        procedure :: stepTransport2d
 !       procedure :: stepTransport3d
end type

!interface setVelocityFn
!    module procedure :: setVel2d
!    module procedure :: setVel3d
!end interface

logical(klog), save :: logInit = .FALSE.
type(Logger) :: log
character(len=28), save :: logKey = 'RK4Log'
integer(kint), parameter :: logLevel = DEBUG_LOGGING_LEVEL

contains

subroutine init(self, mp, mpiParticles, mpiFaces, isCompressible)
    class(Timestepper), intent(inout) :: self
    class(MeshedParticles), intent(in) :: mp
    class(MPIReplicatedData), intent(in) :: mpiParticles, mpiFaces
    logical(klog), intent(in) :: isCompressible
    !
    integer(kint) :: nParticles, nFaces
    
    nParticles = mpiParticles%messageSize(procRank)
    nFaces = mpiFaces%messageSize(procRank)
    
    allocate(self%xIn(nParticles))
    allocate(self%xStage1(nParticles))
    allocate(self%xStage2(nParticles))
    allocate(self%xStage3(nParticles))
    allocate(self%xStage4(nParticles))
    
    allocate(self%yIn(nParticles))
    allocate(self%yStage1(nParticles))
    allocate(self%yStage2(nParticles))
    allocate(self%yStage3(nParticles))
    allocate(self%yStage4(nParticles))
    
    if (mp%geomKind == SPHERE_GEOM) then
        allocate(self%zIn(nParticles))
        allocate(self%zStage1(nParticles))
        allocate(self%zStage2(nParticles))
        allocate(self%zStage3(nParticles))
        allocate(self%zStage4(nParticles))
    endif
    
    if (isCompressible) then
        allocate(self%weightIn(nParticles))
        allocate(self%weightStage1(nParticles))
        allocate(self%weightStage2(nParticles))
        allocate(self%weightStage3(nParticles))
        allocate(self%weightStage4(nParticles))
        allocate(self%rhoIn(nParticles))
        allocate(self%rhoStage1(nParticles))
        allocate(self%rhoStage2(nParticles))
        allocate(self%rhoStage3(nParticles))
        allocate(self%rhoStage4(nParticles))
        allocate(self%areaIn(nFaces))
        allocate(self%areaStage1(nFaces))
        allocate(self%areaStage2(nFaces))
        allocate(self%areaStage3(nFaces))
        allocate(self%areaStage4(nFaces))
    endif
end subroutine

subroutine deleteRK4(self)
    type(TimeStepper), intent(inout) :: self
    if (allocated(self%xIn)) then 
        deallocate(self%xIn)
        deallocate(self%xStage1)
        deallocate(self%xStage2)
        deallocate(self%xStage3)
        deallocate(self%xStage4)
        deallocate(self%yIn)
        deallocate(self%yStage1)
        deallocate(self%yStage2)
        deallocate(self%yStage3)
        deallocate(self%yStage4)        
    endif
    if (allocated(self%zIn)) then
        deallocate(self%zIn)
        deallocate(self%zStage1)
        deallocate(self%zStage2)
        deallocate(self%zStage3)
        deallocate(self%zStage4)
    endif
    if (allocated(self%weightIn)) then
        deallocate(self%weightIn)
        deallocate(self%weightStage1)
        deallocate(self%weightStage2)
        deallocate(self%weightStage3)
        deallocate(self%weightStage4)
        deallocate(self%rhoIn)
        deallocate(self%rhoStage1)
        deallocate(self%rhoStage2)
        deallocate(self%rhoStage3)
        deallocate(self%rhoStage4)
        deallocate(self%areaIn)
        deallocate(self%areaStage1)
        deallocate(self%areaStage2)
        deallocate(self%areaStage3)
        deallocate(self%areaStage4)
    endif
end subroutine

subroutine stepTransport2d(self, mp, mpiParticles, mpiFaces, dt, velFn, divFn)
    class(TransportStepper), intent(inout) :: self
    class(MeshedParticles), intent(inout) :: mp
    class(MPIReplicatedData), intent(in) :: mpiParticles, mpiFaces
    real(kreal), intent(in) :: dt
    procedure(vectorFnOf2dSpaceAndTime) :: velFn
    procedure(scalarFnOf2dSpaceAndTime), optional :: divFn
    !
    integer(kint) :: i, nP, poffset, foffset, j, mpiErrCode
    real(kreal), dimension(2) :: vec
    real(kreal), dimension(3) :: fcent
    real(kreal) :: tin, div
    logical(klog) :: compressible
    
    tin = mp%t
    np = mp%particles%n
    poffset = mpiParticles%startIndex(procRank)
    foffset = mpiFaces%startIndex(procRank)
    compressible = allocated(self%weightIn)
    !
    !RK stage 1
    !
    do i=0, mpiParticles%messageSize(procRank)-1
        vec = velFn(mp%particles%x(i+poffset), mp%particles%y(i+poffset), tin)
        self%xStage1(i+1) = dt*vec(1)
        self%yStage1(i+1) = dt*vec(2)
    enddo
    if (compressible) then
        do i=0, mpiParticles%messageSize(procRank)-1
            div = divFn(mp%particles%x(i+poffset), mp%particles%y(i+poffset), tin)
            self%weightStage1(i+1) = dt * div * mp%particles%weight(i+poffset)
            self%rhoStage1(i+1) = -dt * div * mp%density%comp1(i+poffset)
        enddo
        do i=0, mpiFaces%messageSize(procRank)-1
            if (.not. mp%faces%hasChildren(i+foffset)) then
                fcent = mp%faces%physCentroid(i+foffset, mp%particles)
                div = divFn(fcent(1), fcent(2), tin)
                self%areaStage1(i+1) = dt*div*mp%faces%area(i+foffset)
            endif
        enddo
    endif
    
    !
    !RK stage 2
    !
    tIn = mp%t + 0.5_kreal * dt
    do i=0, mpiParticles%messageSize(procRank)-1
        self%xIn(i+1) = mp%particles%x(i+poffset) + 0.5_kreal*self%xStage1(i+1)
        self%yIn(i+1) = mp%particles%y(i+poffset) + 0.5_kreal*self%yStage1(i+1)
    enddo
    do i=0, mpiParticles%messageSize(procRank)-1
        vec = velFn(self%xIn(i+1), self%yIn(i+1), tin)
        self%xStage2(i+1) = dt*vec(1)
        self%yStage2(i+1) = dt*vec(2)
    enddo
    
    if (compressible) then
        do i=0, mpiParticles%messageSize(procRank)-1
            self%weightIn(i+1) = mp%particles%weight(i+poffset) + 0.5_kreal*self%weightStage1(i+1)
            self%rhoIn(i+1) = mp%density%comp1(i+poffset) + 0.5_kreal*self%rhoStage1(i+1)
        enddo
        do i=0, mpiFaces%messageSize(procRank)-1
            if (.not. mp%faces%hasChildren(i+foffset)) &
            self%areaIn(i+1) = mp%faces%area(i+foffset) + 0.5_kreal*self%areaStage1(i+1)
        enddo
        do i=0, mpiParticles%messageSize(procRank)-1
            div = divFn(self%xIn(i+1), self%yIn(i+1), tin)
            self%weightStage2(i+1) = dt*div*self%weightIn(i+1)
            self%rhoStage2(i+1) = -dt*div*self%rhoIn(i+1)
        enddo
        if (mp%faceKind == TRI_PANEL) then
            do i=0, mpiFaces%messageSize(procRank)-1
                if (.not. mp%faces%hasChildren(i+foffset)) then
                    fcent = dzero
                    do j=1, 3
                        fcent = fcent + [self%xIn(mp%faces%vertices(j,i+foffset)), &
                                        self%yIn(mp%faces%vertices(j,i+foffset)), dzero]
                    enddo
                    fcent = fcent / 3.0_kreal
                    div = divFn(fcent(1), fcent(2), tin)
                    self%areaStage2(i+1) = dt*div*self%areaIn(i+1)
                endif
            enddo
        elseif (mp%faceKind == QUAD_PANEL) then
            do i=0, mpiFaces%messageSize(procRank)-1
                if (.not. mp%faces%hasChildren(i+foffset)) then
                    fcent = dzero
                    do j=1, 4
                        fcent = fcent + [self%xIn(mp%faces%vertices(j,i+foffset)), &
                                        self%yIn(mp%faces%vertices(j,i+foffset)),dzero]
                    enddo
                    fcent = fcent / 4.0_kreal
                    div = divFn(fcent(1), fcent(2), tin)
                    self%areaStage2(i+1) = dt*div*self%areaIn(i+1)
                endif
            enddo
        elseif (mp%faceKind == QUAD_CUBIC_PANEL) then
            do i=0, mpiFaces%messageSize(procRank)-1
                if (.not. mp%faces%hasChildren(i+foffset)) then
                    fcent = dzero
                    do j=1, 4
                        fcent = fcent + [self%xIn(mp%faces%vertices(mod(3*j+9,12)+1,i+foffset)), &
                                         self%yIn(mp%faces%vertices(mod(3*j+9,12)+1,i+foffset)),dzero]
                    enddo
                    fcent = fcent / 4.0_kreal
                    div = divFn(fcent(1), fcent(2), tin)
                    self%areaStage2(i+1) = dt*div*self%areaIn(i+1)
                endif
            enddo
        endif
    endif
    
    !
    !RK stage 3
    !    
    do i=0, mpiParticles%messageSize(procRank)-1
        self%xIn(i+1) = mp%particles%x(i+poffset) + 0.5_kreal * self%xStage2(i+1)
        self%yIn(i+1) = mp%particles%y(i+poffset) + 0.5_kreal * self%yStage2(i+1)
    enddo
    do i=0, mpiParticles%messageSize(procRank)-1
        vec = velFn(self%xIn(i+1), self%yIn(i+1), tin)
        self%xStage3(i+1) = dt*vec(1)
        self%yStage3(i+1) = dt*vec(2)
    enddo
    if (compressible) then
        do i=0, mpiParticles%messageSize(procRank)-1
            self%weightIn(i+1) = mp%particles%weight(i+poffset) + 0.5_kreal*self%weightStage2(i+1)
            self%rhoIn(i+1) = mp%density%comp1(i+poffset) + 0.5_kreal*self%rhoStage2(i+1)
        enddo
        do i=0, mpiFaces%messageSize(procRank)-1
            if (.not. mp%faces%hasChildren(i+foffset)) &
            self%areaIn(i+1) = mp%faces%area(i+foffset) + 0.5_kreal*self%areaStage2(i+1)
        enddo
        do i=0, mpiParticles%messageSize(procRank)-1
            div = divFn(self%xIn(i+1), self%yIn(i+1),tin)
            self%weightStage3(i+1) = dt*div*self%weightIn(i+1)
            self%rhoStage3(i+1) = -dt*div*self%rhoIn(i+1)
        enddo
        if (mp%faceKind == TRI_PANEL) then
            do i=0, mpiFaces%messageSize(procRank)-1
                if (.not. mp%faces%hasChildren(i+foffset)) then
                    fcent = dzero
                    do j=1, 3
                        fcent = fcent + [self%xIn(mp%faces%vertices(j,i+foffset)), &
                                        self%yIn(mp%faces%vertices(j,i+foffset)),dzero]
                    enddo
                    fcent = fcent / 3.0_kreal
                    div = divFn(fcent(1), fcent(2), tin)
                    self%areaStage3(i+1) = dt*div*self%areaIn(i+1)
                endif
            enddo
        elseif (mp%faceKind == QUAD_PANEL) then
            do i=0, mpiFaces%messageSize(procRank)-1
                if (.not. mp%faces%hasChildren(i+foffset)) then
                    fcent = dzero
                    do j=1, 4
                        fcent = fcent + [self%xIn(mp%faces%vertices(j,i+foffset)), &
                                        self%yIn(mp%faces%vertices(j,i+foffset)),dzero]
                    enddo
                    fcent = fcent / 4.0_kreal
                    div = divFn(fcent(1), fcent(2), tin)
                    self%areaStage3(i+1) = dt*div*self%areaIn(i+1)
                endif
            enddo
        elseif (mp%faceKind == QUAD_CUBIC_PANEL) then
            do i=0, mpiFaces%messageSize(procRank)-1
                if (.not. mp%faces%hasChildren(i+foffset)) then
                    fcent = dzero
                    do j=1, 4
                        fcent = fcent + [self%xIn(mp%faces%vertices(mod(3*j+9,12)+1,i+foffset)), &
                                         self%yIn(mp%faces%vertices(mod(3*j+9,12)+1,i+foffset)),dzero]
                    enddo
                    fcent = fcent / 4.0_kreal
                    div = divFn(fcent(1), fcent(2), tin)
                    self%areaStage3(i+1) = dt*div*self%areaIn(i+1)
                endif
            enddo
        endif
    endif
    !
    !RK stage 4
    !
    tin = mp%t + dt
    do i=0, mpiParticles%messageSize(procRank)-1
        self%xIn(i+1) = mp%particles%x(i+poffset) + self%xStage3(i+1)
        self%yIn(i+1) = mp%particles%y(i+poffset) + self%yStage3(i+1)
    enddo
    do i=0, mpiParticles%messageSize(procRank)-1
        vec = velFn(self%xIn(i+1), self%yIn(i+1), tIn)
        self%xStage4(i+1) = dt*vec(1)
        self%yStage4(i+1) = dt*vec(2)
    enddo
    if (compressible) then
        do i=0, mpiParticles%messageSize(procRank)-1
            self%weightIn(i+1) = mp%particles%weight(i+poffset) + self%weightStage3(i+1)
            self%rhoIn(i+1) = mp%particles%weight(i+poffset) + self%weightStage3(i+1)
        enddo
        do i=0, mpiParticles%messageSize(procRank)-1
            div = divFn(self%xIn(i+1), self%yIn(i+1), tIn)
            self%weightStage4(i+1) = dt*div*self%weightIn(i+1)
            self%rhoStage4(i+1) = -dt*div*self%rhoIn(i+1)
        enddo
        do i=0, mpiFaces%messageSize(procRank)-1
            if (.not. mp%faces%hasChildren(i+foffset)) &
                self%areaIn(i+1) = mp%faces%area(i+foffset) + self%areaStage3(i+1)
        enddo
        if (mp%faceKind == TRI_PANEL) then
            do i=0, mpiFaces%messageSize(procRank)-1
                if (.not. mp%faces%hasChildren(i+foffset)) then
                    fcent = dzero
                    do j=1, 3
                        fcent = fcent + [self%xIn(mp%faces%vertices(j,i+foffset)), &
                                        self%yIn(mp%faces%vertices(j,i+foffset)),dzero]
                    enddo
                    fcent = fcent / 3.0_kreal
                    div = divFn(fcent(1), fcent(2), tin)
                    self%areaStage4(i+1) = dt*div*self%areaIn(i+1)
                endif
            enddo
        elseif (mp%faceKind == QUAD_PANEL) then
            do i=0, mpiFaces%messageSize(procRank)-1
                if (.not. mp%faces%hasChildren(i+foffset)) then
                    fcent = dzero
                    do j=1, 4
                        fcent = fcent + [self%xIn(mp%faces%vertices(j,i+foffset)), &
                                        self%yIn(mp%faces%vertices(j,i+foffset)),dzero]
                    enddo
                    fcent = fcent / 4.0_kreal
                    div = divFn(fcent(1), fcent(2), tin)
                    self%areaStage4(i+1) = dt*div*self%areaIn(i+1)
                endif
            enddo
        elseif (mp%faceKind == QUAD_CUBIC_PANEL) then
            do i=0, mpiFaces%messageSize(procRank)-1
                if (.not. mp%faces%hasChildren(i+foffset)) then
                    fcent = dzero
                    do j=1, 4
                        fcent = fcent + [self%xIn(mp%faces%vertices(mod(3*j+9,12)+1,i+foffset)), &
                                         self%yIn(mp%faces%vertices(mod(3*j+9,12)+1,i+foffset)), dzero]
                    enddo
                    fcent = fcent / 4.0_kreal
                    div = divFn(fcent(1), fcent(2), tin)
                    self%areaStage4(i+1) = dt*div*self%areaIn(i+1)
                endif
            enddo
        endif
    endif
    !
    ! RK update
    !
    do i=0, mpiParticles%messageSize(procRank)-1
        mp%particles%x(i+poffset) = mp%particles%x(i+poffset) + self%xStage1(i+1)/6.0_kreal + &
            self%xStage2(i+1)/3.0_kreal + self%xStage3(i+1)/3.0_kreal + self%xStage4(i+1)/6.0_kreal
        mp%particles%y(i+poffset) = mp%particles%y(i+poffset) + self%yStage1(i+1)/6.0_kreal + &
            self%yStage2(i+1)/3.0_kreal + self%yStage3(i+1)/3.0_kreal + self%yStage4(i+1)/6.0_kreal
    enddo
    if (compressible) then
        do i=0, mpiParticles%messageSize(procRank)-1
            mp%particles%weight(i+poffset) = mp%particles%weight(i+poffset) + self%weightStage1(i+1)/6.0_kreal + &
                self%weightStage2(i+1)/3.0_kreal + self%weightStage3(i+1)/3.0_kreal + self%weightStage4(i+1)/6.0_kreal
            mp%density%comp1(i+poffset) = mp%density%comp1(i+poffset) + self%rhoStage1(i+1)/6.0_kreal + &
                self%rhoStage2(i+1)/3.0_kreal + self%rhoStage3(i+1)/3.0_kreal + self%rhoStage4(i+1)/6.0_kreal
        enddo
        do i=0, mpiFaces%messageSize(procRank)-1
            if (.not. mp%faces%hasChildren(i+foffset)) &
                mp%faces%area(i+foffset) = mp%faces%area(i+foffset) + self%areaStage1(i+1)/6.0_kreal + &
                self%areaStage2(i+1)/3.0_kreal + self%areaStage3(i+1)/3.0_kreal + self%areaStage4(i+1)/6.0_kreal
        enddo
    endif
    
    !
    !   broadcast
    !
    do i=0, numProcs-1
        call MPI_BCAST(mp%particles%x(mpiParticles%startIndex(i):mpiParticles%endIndex(i)), &
            mpiParticles%messageSize(i), MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)
        call MPI_BCAST(mp%particles%y(mpiParticles%startIndex(i):mpiParticles%endIndex(i)), &
            mpiParticles%messageSize(i), MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)
        if (compressible) then
            call MPI_BCAST(mp%particles%weight(mpiParticles%startIndex(i):mpiParticles%endIndex(i)), &
                mpiParticles%messageSize(i), MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)
            call MPI_BCAST(mp%density%comp1(mpiParticles%startIndex(i):mpiParticles%endIndex(i)), &
                mpiParticles%messageSize(i), MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)
            call MPI_BCAST(mp%faces%area(mpiFaces%startIndex(i):mpiFaces%endIndex(i)), &
                mpiFaces%messageSize(i), MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)
        endif
    enddo
    
    tIn = mp%t + dt
    do i=1, mp%particles%N
        vec = velFn(mp%particles%x(i), mp%particles%y(i), tIn)
        mp%velocity%comp1(i) = vec(1)
        mp%velocity%comp2(i) = vec(2)
    enddo
    mp%t = mp%t + dt
end subroutine

!subroutine stepTransport3d(self, mpiParticles, dt, velFn, divFn)
!    class(TransportStepper), intent(inout) :: self
!    class(MPIReplicatedData), intent(in) :: mpiParticles
!    real(kreal), intent(in) :: dt
!    procedure(vectorFnOf3dSpaceAndTime) :: velFn
!    procedure(scalarFnOf3dSpaceAndTime), optional :: divFn
!end subroutine

!subroutine setVel2d(self, velFn)
!    class(TransportStepper), intent(inout) :: self
!    procedure(vectorFnOf2DSpaceAndTime), intent(in) :: velFn
!    self%velFn2d => velFn
!    nullify(self%velFn3d)
!end subroutine
!
!subroutine setVel3d(self, velFn)
!    class(TransportStepper), intent(inout) :: self
!    procedure(vectorFnOf3DSpaceAndTime), intent(in) :: velFn
!    self%velFn3d => velFn
!    nullify(self%velFn2d)
!end subroutine

subroutine InitLogger(aLog,rank)
	type(Logger), intent(out) :: aLog
	integer(kint), intent(in) :: rank
	write(logKey,'(A,A,I0.3,A)') trim(logKey),'_',rank,' : '
	call New(aLog,logLevel)
	logInit = .TRUE.
end subroutine

!>@}
end module