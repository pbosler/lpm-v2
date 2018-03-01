module PolyMeshOOModule

use NumberKindsModule
use STDIntVectorModule
use OutputWriterModule
use LoggerModule
use ParticlesOOModule
use EdgesOOModule
use FacesOOModule
use PlaneGeomModule
use SphereGeomModule

implicit none
private

public PolyMesh2d
!public PlanarTriMesh, PlanarQuadMesh ! Linear Faces & edges
!public IcosTriSphereMesh, CubedSphereMesh ! Linear faces & edges
!public PlanarCubicQuadMesh ! Cubic faces & edges
!public CubedSphereCubicMesh ! Cubic faces & edges

type PolyMesh2d
    character(len=56) :: mesh_type
    class(Particles), pointer :: particles
    class(Edges), pointer :: edges
    class(Faces), pointer :: faces
    integer(kint) :: meshSeed = 0
    integer(kint) :: initNest = 0
    integer(kint) :: maxNest = 0
    integer(kint) :: amrLimit = 0
    real(kreal) :: t = dzero
    real(kreal) :: maxRadius = 1.0_kreal
    
    contains
        procedure :: init
        procedure :: createFromSeed
!        procedure :: copy
!        procedure :: refine
        
        final :: deleteMesh
end type

logical(klog), save :: logInit = .FALSE.
type(Logger) :: log
character(len=28), save :: logKey = 'PolyMesh2d'
integer(kint), parameter :: logLevel = TRACE_LOGGING_LEVEL
character(len=MAX_STRING_LENGTH) :: logstring

contains

subroutine init(self, mesh_type, initNest, maxNest, amrLimit, maxRadius)
    class(PolyMesh2d), intent(inout) :: self
    character(len=*), intent(in) ::  mesh_type
    integer(kint), intent(in) :: initNest, maxNest, amrLimit
    real(kreal), intent(in) :: maxRadius
    
    if ( .NOT. logInit) call InitLogger(log, procRank)
    
    self%mesh_type = trim(mesh_type)
    self%initNest = initNest
    self%maxNest = maxNest
    self%amrLimit = amrLimit
    self%maxRadius = maxRadius
    
    self%particles => null()
    self%edges => null()
    self%faces => null()
end subroutine

subroutine createFromSeed(self)
    class(PolyMesh2d), intent(inout) :: self
    !
    character(len=128) :: seedfilename
    
    if (associated(self%particles)) deallocate(self%particles)
    if (associated(self%edges)) deallocate(self%edges)
    if (associated(self%faces)) deallocate(self%faces)
    
    allocate( Particles :: self%particles)
    
    if (self%mesh_type == "planar_tri") then
         allocate( LinearEdges :: self%edges)
         allocate( TriLinearFaces :: self%faces)
    elseif (self%mesh_type == "planar_quad") then
        allocate(LinearEdges :: self%edges)
        allocate(QuadLinearFaces :: self%faces)
    elseif (self%mesh_type == "planar_cubic_quad") then
        allocate(CubicEdges :: self%edges)
        allocate(QuadCubicFaces :: self%faces)
    elseif (self%mesh_type == "icosTriSphere") then
        allocate(LinearEdges :: self%edges)
        allocate(TriLinearFaces :: self%faces)
    elseif (self%mesh_type == "cubedSphere") then
        allocate(LinearEdges :: self%edges)
        allocate(QuadLinearFaces :: self%faces)
    elseif (self%mesh_type == "cubedSphereCubic") then
        allocate(CubicEdges :: self%edges)
        allocate(QuadCubicFaces :: self%faces)
    endif
end subroutine

subroutine deleteMesh(self)
    type(PolyMesh2d), intent(inout) :: self
    deallocate(self%particles)
    nullify(self%particles)
    deallocate(self%edges)
    nullify(self%edges)
    deallocate(self%faces)
    nullify(self%faces)
end subroutine

!> @brief Initializes a logger for the Faces module
!> 
!> Output is controlled both by message priority and by MPI Rank
!> @param[out] aLog Target Logger object
!> @param[in] rank Rank of this processor
subroutine InitLogger(aLog,rank)
! Initialize a logger for this module and processor
	type(Logger), intent(out) :: aLog
	integer(kint), intent(in) :: rank
	write(logKey,'(A,A,I0.2,A)') trim(logKey),'_',rank,' : '
	if ( rank == 0 ) then
		call New(aLog,logLevel)
	else
		call New(aLog,ERROR_LOGGING_LEVEL)
	endif
	logInit = .TRUE.
end subroutine

end module