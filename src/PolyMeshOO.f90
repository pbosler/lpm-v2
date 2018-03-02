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

type PolyMesh2d
    character(len=56) :: mesh_type
    class(Particles), pointer :: particles
    class(Edges), pointer :: edges
    class(Faces), pointer :: faces
    integer(kint) :: faceKind = 0
    integer(kint) :: geomKind = 0
    integer(kint) :: initNest = 0
    integer(kint) :: maxNest = 0
    integer(kint) :: amrLimit = 0
    real(kreal) :: t = dzero
    real(kreal) :: maxRadius = 1.0_kreal
    
    contains
        procedure :: init
        procedure, private :: getSeed
        procedure, private :: readSeedFile
        procedure, private :: nVerticesInMesh
        procedure, private :: nEdgesInMesh
        procedure, private :: nFacesInMesh
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
    !
    integer(kint) :: i,j, nFacesOld, startIndex
    
    if ( .NOT. logInit) call InitLogger(log, procRank)
    
    self%mesh_type = trim(mesh_type)
    self%initNest = initNest
    self%maxNest = maxNest
    self%amrLimit = amrLimit
    self%maxRadius = maxRadius
    
    self%particles => null()
    self%edges => null()
    self%faces => null()
    
    call self%getSeed()
    
    startIndex = 1
    do i=1, initNest
        nFacesOld = self%faces%N
        do j=startIndex, nFacesOld
            call self%faces%divide(j, self%particles, self%edges)
        enddo
        startIndex = nFacesOld + 1
    enddo
end subroutine

pure function nVerticesInMesh(self, nestLevel)
    integer(kint) :: nVerticesInMesh
    class(PolyMesh2d), intent(in) :: self
    integer(kint), intent(in) :: nestLevel
    !
    integer(kint) :: i
    nVerticesInMesh = 0
    select case (self%mesh_type)
        case ("planar_tri")
            do i = 2**nestLevel +1, 2**(nestLevel+1)
                nVerticesInMesh = nVerticesInMesh + i
            enddo
            nVerticesInMesh = 2*nVerticesInMesh + 2**(nestLevel+1) + 1
        case ("planar_quad")
            nVerticesInMesh = 3
            do i=1,nestLevel
                nVerticesInMesh = nVerticesInMesh + 2**i
            enddo
        case ("planar_cubic_quad")
        case ("icosTriSphere")
            nVerticesInMesh = 2 + 10*4**nestLevel
        case ("cubedSphere")
            nVerticesInMesh = 2 + 6*4**nestLevel
    end select
end function

pure function nFacesInMesh(self, nestLevel)
    integer(kint) :: nFacesInMesh
    class(PolyMesh2d),intent(in) :: self
    integer(kint), intent(in) :: nestLevel
    !
    integer(kint) :: i
    
    nFacesInMesh = 0
    select case (self%mesh_type)
        case ("planar_tri")
            nFacesInMesh = 6*4**nestLevel
        case ("planar_quad")
            nFacesInMesh = 4*4**nestLevel
        case ("planar_cubic_quad")
            nFacesInMesh = 4*4**nestLevel
        case ("icosTriSphere")
            nFacesInMesh = 20*4**nestLevel
        case ("cubedSphere")
            nFacesInMesh = 6*4**nestLevel
    end select
end function

pure function nEdgesInMesh(self, nVerts, nFaces)
    integer(kint) :: nEdgesInMesh
    class(PolyMesh2d), intent(in) :: self
    integer(kint), intent(in) :: nVerts, nFaces
    
    nEdgesInMesh = 0
    select case (self%geomKind)
        case (PLANAR_GEOM)
            nEdgesInMesh = nFaces + nVerts -1
        case (SPHERE_GEOM)
            nEdgesInMesh = nFaces + nVerts -2
    end select
end function

subroutine getSeed(self)
    class(PolyMesh2d), intent(inout) :: self
    !
    character(len=128) :: seedfilename
    integer(kint) :: i, nSeedParticles, nSeedEdges, nSeedFaces, nSeedVerts, nv, nc
    integer(kint), allocatable :: seedEdgeOrigs(:), seedEdgeDests(:), seedEdgeLefts(:), seedEdgeRights(:)
    integer(kint), allocatable :: seedEdgeInts(:,:), seedFaceEdges(:,:), seedFaceVerts(:,:), seedFaceCenters(:,:)
    real(kreal), allocatable :: seedXYZ(:,:)
    integer(kint) :: nMaxParticles, nMaxEdges, nMaxFaces
    
    if (associated(self%particles)) deallocate(self%particles)
    if (associated(self%edges)) deallocate(self%edges)
    if (associated(self%faces)) deallocate(self%faces)
    
    allocate( Particles :: self%particles)
    
    if (self%mesh_type == "planar_tri") then
         allocate( LinearEdges :: self%edges)
         allocate( TriLinearFaces :: self%faces)
         
        seedfilename = "triHexSeed.dat"
         
        nSeedParticles = 13
        nSeedEdges = 12
        nSeedFaces = 6
        nSeedVerts = 7
        
        self%faceKind = TRI_PANEL
        self%geomKind = PLANAR_GEOM
        nv = 3
        nc=1
    elseif (self%mesh_type == "planar_quad") then
        allocate(LinearEdges :: self%edges)
        allocate(QuadLinearFaces :: self%faces)
        
        seedfilename = "quadRectSeed.dat"
        nSeedParticles = 13
        nSeedEdges = 12
        nSeedFaces = 4
        nv = 4
        nc=1
        self%faceKind = QUAD_PANEL
        self%geomKind = PLANAR_GEOM
    elseif (self%mesh_type == "planar_cubic_quad") then
        allocate(CubicEdges :: self%edges)
        allocate(QuadCubicFaces :: self%faces)
        seedfilename = "quadCubicSeed.dat"
        nSeedParticles = 49
        nSeedEdges = 12
        nSeedFaces = 4
        nv=12
        nc=4
        self%faceKind= QUAD_CUBIC_PANEL
        self%geomKind = PLANAR_GEOM
    elseif (self%mesh_type == "icosTriSphere") then
        allocate(LinearEdges :: self%edges)
        allocate(TriLinearFaces :: self%faces)
        
        seedfilename = "icosTriSeed.dat"
        nSeedParticles = 32
        nSeedEdges = 30
        nSeedFaces = 20
        nv=3
        nc=1
        self%faceKind = TRI_PANEL
        self%geomKind = SPHERE_GEOM
    elseif (self%mesh_type == "cubedSphere") then
        allocate(LinearEdges :: self%edges)
        allocate(QuadLinearFaces :: self%faces)
        
        seedfilename = "cubedSphereSeed.dat"
        nSeedParticles = 14
        nSeedEdges = 12
        nSeedFaces = 6
        nv = 4
        nc=1
        self%faceKind = QUAD_PANEL
        self%geomKind = SPHERE_GEOM
    elseif (self%mesh_type == "cubedSphereCubic") then
        allocate(CubicEdges :: self%edges)
        allocate(QuadCubicFaces :: self%faces)
        
        self%faceKind = QUAD_CUBIC_PANEL
        self%geomKind = SPHERE_GEOM
        nv = 12
        nc = 4
    endif
    
    nMaxParticles = self%nVerticesInMesh(self%maxNest) + self%nFacesInMesh(self%maxNest)
    nMaxFaces = 0
    nMaxEdges = 0
    do i=0, self%maxNest
        nMaxFaces = nMaxFaces + self%nFacesInMesh(i)
        nMaxEdges = nMaxEdges + self%nEdgesInMesh(self%nVerticesInMesh(i), self%nFacesInMesh(i))
    enddo
    
    call self%particles%init(nMaxParticles, self%geomKind)
    call self%edges%init(nMaxEdges)
    call self%faces%init(self%faceKind, nMaxFaces)
    
    allocate(seedXYZ(3,nSeedParticles))
    allocate(seedEdgeOrigs(nSeedEdges))
    allocate(seedEdgeDests(nSeedEdges))
    allocate(seedEdgeLefts(nSeedEdges))
    allocate(seedEdgeRights(nSeedEdges))
    allocate(seedEdgeInts(2,nSeedEdges))
    
    if (self%faceKind == TRI_PANEL) then
        allocate(seedFaceVerts(3,nSeedFaces))
        allocate(seedFaceEdges(3,nSeedFaces))
        allocate(seedFaceCenters(1,nSeedFaces))
    elseif (self%faceKind == QUAD_PANEL) then
        allocate(seedFaceVerts(4,nSeedFaces))
        allocate(seedFaceEdges(4,nSeedFaces))
        allocate(seedFaceCenters(1,nSeedFaces))
    elseif (self%faceKind == QUAD_CUBIC_PANEL) then
        allocate(seedFaceVerts(12,nSeedFaces))
        allocate(seedFaceEdges(4,nSeedFaces))
        allocate(seedFaceCenters(4,nSeedFaces))
    endif
    
    call self%readSeedFile(seedfilename, nSeedParticles, nSeedEdges, nSeedFaces, nv, nc, &
        seedXYZ, seedEdgeOrigs, seedEdgeDests, seedEdgeLefts, seedEdgeRights, seedEdgeInts, &
        seedFaceVerts, seedFaceEdges, seedFaceCenters)
    
    seedXYZ = self%maxRadius * seedXYZ
    
    ! initialize mesh
    do i=1, nSeedParticles
        call self%particles%insert(seedXYZ(:,i), seedXYZ(:,i))
    enddo
    if ( self%particles%N /= nSeedParticles ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL, logkey//" initMeshFromSeed ERROR : "," particles%N.")
	endif
	
	do i=1, nSeedEdges
	    call self%edges%insert(seedEdgeOrigs(i), seedEdgeDests(i), seedEdgeLefts(i), seedEdgeRights(i), seedEdgeInts(:,i))
	enddo
	if ( self%edges%N /= nSeedEdges ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL, logkey//" initMeshFromSeed ERROR : "," edges%N.")
	endif
	
	do i=1, nSeedFaces
	    call self%faces%insert(seedFaceCenters(:,i), seedFaceVerts(:,i), seedFaceEdges(:,i))
	enddo
end subroutine

subroutine readSeedFile(self, seedfilename, nSeedParticles, nSeedEdges, nSeedFaces, nv, nc, &
    seedXYZ, seedEdgeOrigs, seedEdgeDests, seedEdgeLefts, seedEdgeRights, seedEdgeInts, &
    seedFaceVerts, seedFaceEdges, seedFaceCenters)
    class(PolyMesh2d), intent(in) :: self
    character(len=*), intent(in) :: seedfilename
    integer(kint), intent(in) :: nSeedParticles, nSeedEdges, nSeedFaces, nv, nc
    real(kreal), dimension(3,nSeedParticles), intent(out) :: seedXYZ
    integer(kint), dimension(nSeedEdges), intent(out) :: seedEdgeOrigs, seedEdgeDests, seedEdgeLefts, seedEdgeRights
    integer(kint), dimension(2,nSeedEdges), intent(out) :: seedEdgeInts
    integer(kint), dimension(nv, nSeedFaces), intent(out) :: seedFaceVerts, seedFaceEdges
    integer(kint), dimension(nc, nSeedFaces), intent(out) :: seedFaceCenters
    !
    character(len=128) :: linestring
    integer(kint) :: i, readStat
    
    open(unit=READ_UNIT, file=seedfilename, status='OLD', action='READ', iostat=readStat)
    if (readStat /= 0) then
        call LogMessage(log, ERROR_LOGGING_LEVEL, trim(logkey)//" io error, cannot read file ", seedfilename)
        return
    endif
    
    !
    !   read particle starting positions
    !
    seedXYZ = dzero
    read(READ_UNIT, *) linestring ! x    y    z
    if (self%geomKind == PLANAR_GEOM) then
        do i=1, nSeedParticles
            read(READ_UNIT,*) seedXYZ(1:2,i)
        enddo
    elseif(self%geomKind == SPHERE_GEOM) then
        do i=1, nSeedParticles
            read(READ_UNIT,*) seedXYZ(:,i)
        enddo
    endif
    
    !
    !   read root edges
    !
    read(READ_UNIT,*) linestring ! edgeO edgeD edgeLeft edgeRight edgeInts
    seedEdgeOrigs = 0
    seedEdgeDests = 0
    seedEdgeLefts = 0
    seedEdgeRights = 0
    seedEdgeInts = 0
    if (self%faceKind /= QUAD_CUBIC_PANEL) then
        do i=1, nSeedEdges
            read(READ_UNIT,*) seedEdgeOrigs(i), seedEdgeDests(i), seedEdgeLefts(i), seedEdgeRights(i)        
        enddo
    else
        do i=1, nSeedEdges
            read(READ_UNIT,*) seedEdgeOrigs(i), seedEdgeDests(i), seedEdgeLefts(i), seedEdgeRights(i), seedEdgeInts(:,i) 
        enddo
    endif
    !
    !   change for fortran base 1 indexing
    !
    seedEdgeOrigs = seedEdgeOrigs + 1
    seedEdgeDests = seedEdgeDests + 1
    seedEdgeLefts = seedEdgeLefts + 1
    seedEdgeRights = seedEdgeRights + 1
    seedEdgeInts = seedEdgeInts + 1
    
    !
    !   read root faces
    !
    read(READ_UNIT,*) linestring ! faceverts
    seedFaceCenters = 0
    seedFaceEdges = 0
    seedFaceVerts = 0
    do i=1, nSeedFaces
        read(READ_UNIT,*) seedFaceVerts(:,i)
    enddo
    read(READ_UNIT,*) linestring    ! faceedges
    do i=1, nSeedFaces
        read(READ_UNIT,*) seedFaceEdges(:,i)
    enddo
    
    if (self%faceKind == QUAD_CUBIC_PANEL) then
        read(READ_UNIT, *) linestring ! facecenters
        do i=1, nSeedFaces
            read(READ_UNIT,*) seedFaceCenters(:,i)
        enddo
    endif
    seedFaceCenters = seedFaceCenters + 1
    seedFaceEdges = seedFaceEdges + 1
    seedFaceVerts = seedFaceVerts +1
    
    close(READ_UNIT)
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