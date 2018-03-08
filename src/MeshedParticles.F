module MeshedParticlesModule

use NumberKindsModule
use UtilitiesModule
use LoggerModule
use ParticlesOOModule
use EdgesOOModule
use FacesOOModule
use PolyMeshOOModule
use FieldOOModule
use MPIReplicatedDataModule

implicit none 

private

public MeshedParticles

type, extends(PolyMesh2d) :: MeshedParticles
    type(Field) :: density
    type(Field) :: velocity
    type(Field), dimension(:), allocatable :: tracers
    integer(kint) :: nTracers = 0
    contains
        procedure :: init => initTransport
        procedure :: addTracers
        procedure :: set2DVelocityWithFunction
        procedure :: set3DVelocityWithFunction
        procedure :: set2DDensityWithFunction
        procedure :: set3DDensityWithFunction
        procedure :: writeVTK => writeTransportVTK
!        procedure :: computeVelocity2d
!        procedure :: computeVelocity3d
        
        final :: deleteMP
end type

!interface
!    subroutine computeVelocity2d(u, v, x, y, mpiParticles)
!        import :: kreal
!        import :: MPIReplicatedData
!        implicit none
!        real(kreal), dimension(:), intent(out) :: u, v
!        real(kreal), dimension(:), intent(in) :: x, y
!        class(MPIReplicatedData), intent(in) :: mpiParticles
!    end subroutine
!end interface
!
!interface
!    subroutine computeVelocity2d(u, v, w, x, y, z, mpiParticles)
!        import :: kreal
!        import :: MPIReplicatedData
!        implicit none
!        real(kreal), dimension(:), intent(out) :: u, v, w
!        real(kreal), dimension(:), intent(in) :: x, y, z
!        class(MPIReplicatedData), intent(in) :: mpiParticles
!    end subroutine
!end interface

logical(klog), save :: logInit = .FALSE.
type(Logger) :: log
character(len=28), save :: logKey = 'MeshedParticlesLog'
integer(kint), parameter :: logLevel = DEBUG_LOGGING_LEVEL

contains

subroutine initTransport(self, mesh_type, initNest, maxNest, amrLimit, maxRadius)
    class(MeshedParticles), intent(inout) :: self
    character(len=*), intent(in) :: mesh_Type
    integer(kint), intent(in) :: initNest
    integer(kint), intent(in) :: maxNest
    integer(kint), intent(in) :: amrLimit
    real(kreal), intent(in) :: maxRadius
    
    
    if (.not. logInit) call InitLogger(log, procRank)
    call LogMessage(log, DEBUG_LOGGING_LEVEL, trim(logkey)//" ", "entering init MeshedParticles.")
    call self%polymesh2d%init(mesh_type, initNest, maxNest, amrLimit, maxRadius)
    call self%density%init(1, self%particles%N_Max, "density")
    if (self%geomKind == PLANAR_GEOM) then
        call self%velocity%init(2, self%particles%N_Max, "velocity")
    elseif (self%geomKind == SPHERE_GEOM) then
        call self%velocity%init(3, self%particles%N_Max, "velocity")
    endif
end subroutine

subroutine addTracers(self, nTracers, tracerDims)
    class(MeshedParticles), intent(inout) :: self
    integer(kint), intent(in) :: nTracers
    integer(kint), dimension(:), intent(in), optional :: tracerDims
    !
    character(len=10) :: trString
    integer(kint) :: i
    
    if (allocated(self%tracers)) then
        call LogMessage(log, ERROR_LOGGING_LEVEL, trim(logkey)//" addTracers error: ", "tracers already allocated.")
        return
    endif
    
    if (present(tracerDims)) then        
        if (size(tracerDims) == nTracers) then
            ! specify the dimensions of each tracer (vector or scalar)
            allocate(self%tracers(nTracers))
            do i=1, nTracers
                write(trString,'(A,I1)') 'tracer', i
                call self%tracers(i)%init(tracerDims(i), self%particles%N_Max, trString)
            enddo
        else
            call LogMessage(log, ERROR_LOGGING_LEVEL, trim(logKey)//" ",&
                 "tracers error, nTracers /= size(tracerDims)")
        endif
    else
        ! assume all tracers are scalar
        allocate(self%tracers(nTracers))
        do i=1,nTracers
            write(trString,'(A,I1)') 'tracer', i
            call self%tracers(i)%init(1, self%particles%N_Max, trString)
        enddo
    endif
    self%nTracers = nTracers
end subroutine

subroutine deleteMP(self)
    type(MeshedParticles), intent(inout) :: self
    if (allocated(self%tracers)) deallocate(self%tracers)
end subroutine

subroutine set2DVelocityWithFunction(self, velFn, t)
    class(MeshedParticles), intent(inout) :: self
    procedure(vectorFnOf2DSpaceAndTime) :: velFn
    real(kreal), intent(in) :: t
    !
    integer(kint) :: i
    
    if (.not. self%velocity%n == self%particles%n) self%velocity%n = self%particles%N
    
    do i=1, self%particles%n
        call self%velocity%replaceVector(i, velFn(self%particles%x(i), self%particles%y(i), t))
    enddo
end subroutine

subroutine set2DDensityWithFunction(self, rhoFn)
    class(MeshedParticles), intent(inout) :: self
    procedure(scalarFnOf2DSpace) :: rhoFn
    !
    integer(kint) :: i
    if (.not. self%density%n == self%particles%n) self%density%n = self%particles%n
    do i=1, self%particles%n
        call self%density%replaceScalar(i, rhoFn(self%particles%x(i), self%particles%y(i)))
    enddo
end subroutine

subroutine set3DDensityWithFunction(self, rhoFn)
    class(MeshedParticles), intent(inout) :: self
    procedure(scalarFnOf3DSpace) :: rhoFn
    !
    integer(kint) :: i
    if (.not. self%density%n == self%particles%n) self%density%n = self%particles%n
    do i=1, self%particles%n
        call self%density%replaceScalar(i, rhoFn(self%particles%x(i), self%particles%y(i), self%particles%z(i)))
    enddo
end subroutine

subroutine set3DVelocityWithFunction(self, velFn)
    class(MeshedParticles), intent(inout) :: self
    procedure(vectorFnOf3DSpace) :: velFn
    !
    integer(kint) :: i
    
    if (.not. self%velocity%n == self%particles%n) self%velocity%n = self%particles%N
    
    do i=1, self%particles%n
        call self%velocity%replaceVector(i, velFn(self%particles%x(i), self%particles%y(i), self%particles%z(i)))
    enddo
end subroutine

subroutine writeTransportVTK(self, fileunit)
    class(MeshedParticles), intent(in) :: self
    integer(kint), intent(in) :: fileunit
    !
    character(len=MAX_STRING_LENGTH) :: datastring
    integer(kint) :: i, sCount, vCount
    integer(kint), allocatable :: tracerIds(:)
    
    allocate(tracerIds(self%nTracers))
    tracerIds = 0
    if (self%geomKind == PLANAR_GEOM) then
        call self%writePlanarVTKStartWithHeightFromField(self%tracers(1), fileunit)
    else
        call self%writeVTKSerialStartXML(fileunit)
    endif
    sCount=0
    vCount=1
    do i=1,self%nTracers
        if (self%tracers(i)%nDim==1) then
            sCount = sCount + 1
            tracerIds(i) = i
        else
            tracerIds(i) = -i
            vCount = vCount + 1
        endif
    enddo
    write(fileunit,'(A,A)') '     <PointData Scalars="density" Vectors="velocity"','>'
    call self%velocity%writeVtkPointDataXML(fileunit)
    do i=1, self%ntracers
        call self%tracers(i)%writeVtkPointDataXML(fileunit)
    enddo
    write(fileunit,'(A)') '     </PointData>'
    call self%writeVTKSerialEndXML(fileunit)
    deallocate(tracerIds)
end subroutine



subroutine InitLogger(aLog,rank)
	type(Logger), intent(out) :: aLog
	integer(kint), intent(in) :: rank
	write(logKey,'(A,A,I0.3,A)') trim(logKey),'_',rank,' : '
	call New(aLog,logLevel)
	logInit = .TRUE.
end subroutine

end module

