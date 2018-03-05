module PolyMeshOOModule

use NumberKindsModule
use UtilitiesModule
use STDIntVectorModule
use OutputWriterModule
use LoggerModule
use ParticlesOOModule
use EdgesOOModule
use FacesOOModule
use PlaneGeomModule
use SphereGeomModule
use FieldOOModule

implicit none
private

public PolyMesh2d

type PolyMesh2d
    integer(kint) :: meshSeed
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
        procedure :: logStats
        procedure, private :: getSeed
        procedure, private :: readSeedFile
        procedure, private :: nVerticesInMesh
        procedure, private :: nEdgesInMesh
        procedure, private :: nFacesInMesh
        procedure :: writeMatlab
        procedure :: writeVTKSerialStartXML
        procedure :: writeVTKSerialEndXML
        procedure :: integrateScalar
!        procedure :: copy
!        procedure :: refine

        final :: deleteMesh
end type

logical(klog), save :: logInit = .FALSE.
type(Logger) :: log
character(len=28), save :: logKey = 'PolyMesh2d'
integer(kint), parameter :: logLevel = DEBUG_LOGGING_LEVEL
character(len=MAX_STRING_LENGTH) :: logstring

contains

function meshKindFromString(str)
    integer(kint) :: meshKindFromString
    character(len=*), intent(in) :: str

    meshKindFromString = 0
    if (trim(str) == "planar_tri") then
        meshKindFromString = TRI_HEX_SEED
    elseif (trim(str) == "planar_quad") then
        meshKindFromString = QUAD_RECT_SEED
    elseif (trim(str) == "planar_cubic_quad") then
        meshKindFromString = CUBIC_PLANE_SEED
    elseif (trim(str) == "cubed_sphere") then
        meshKindFromString = CUBED_SPHERE_SEED
    elseif (trim(str) == "icos_tri_sphere") then
        meshKindFromString = ICOS_TRI_SPHERE_SEED
    else
        call LogMessage(log, ERROR_LOGGING_LEVEL, "invalid mesh_type string: ", str)
    endif
end function

subroutine init(self, mesh_type, initNest, maxNest, amrLimit, maxRadius)
    class(PolyMesh2d), intent(inout) :: self
    character(len=*), intent(in) ::  mesh_type
    integer(kint), intent(in) :: initNest, maxNest, amrLimit
    real(kreal), intent(in) :: maxRadius
    !
    integer(kint) :: i,j, nFacesOld, startIndex

    if ( .NOT. logInit) call InitLogger(log, procRank)

    self%meshSeed = meshKindFromString(mesh_type)
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

pure function integrateScalar(self, scalar)
    real(kreal) :: integrateScalar
    class(PolyMesh2d), intent(in) ::  self
    class(Field), intent(in) :: scalar
    !
    integer(kint) :: i, j, pindex
    
    integrateScalar = dzero
    if (self%faceKind == TRI_PANEL .or. self%faceKind == QUAD_PANEL) then
        do i=1, self%faces%N
            if (.not. self%faces%hasChildren(i)) then
                pindex = self%faces%centerParticles(1,i)
                integrateScalar = integrateScalar + self%particles%weight(pindex) * scalar%comp1(pindex)
            endif
        enddo
    elseif (self%faceKind == QUAD_CUBIC_PANEL) then
        do i=1, self%faces%N
            if (.not. self%faces%hasChildren(i)) then
                do j=1,12
                    integrateScalar = integrateScalar + self%particles%weight(self%faces%vertices(j,i)) * &
                        scalar%comp1(self%faces%vertices(j,i))
                enddo
                do j=1,4
                    integrateScalar = integrateScalar + self%particles%weight(self%faces%centerParticles(j,i)) * &
                        scalar%comp1(self%faces%centerParticles(j,i))
                enddo
            endif
        enddo
    endif
end function

subroutine logStats(self, aLog)
    class(PolyMesh2d), intent(in) :: self
    type(Logger), intent(inout) :: aLog
    !
    call LogMessage(aLog, TRACE_LOGGING_LEVEL, logKey, " PolyMesh2d stats:")
    call StartSection(aLog)
    call LogMessage(aLog, TRACE_LOGGING_LEVEL, "meshSeed = ", self%meshSeed)
    call LogMessage(aLog, TRACE_LOGGING_LEVEL, "initNest = ", self%initNest)
    call LogMessage(aLog, TRACE_LOGGING_LEVEL, "maxNest = ", self%maxNest)
    call LogMessage(alog, TRACE_LOGGING_LEVEL, "amrLimit = ", self%amrLimit)
    call LogMessage(alog, TRACE_LOGGING_LEVEL, "maxRadius = ", self%maxRadius)
    call self%particles%logStats(aLog)
    call self%edges%logStats(alog)
    call self%faces%logStats(aLog)
    call EndSection(aLog)
end subroutine

subroutine writeMatlab(self, fileunit)
    class(PolyMesh2d), intent(in) :: self
    integer(kint), intent(in) :: fileunit

    call self%particles%writeMatlab(fileunit)
    call self%edges%writeMatlab(fileunit)
    call self%faces%writeMatlab(fileunit)
end subroutine

subroutine writeVTKSerialStartXML(self, fileunit)
    class(PolyMesh2d), intent(in) :: self
    integer(kint), intent(in) :: fileunit
    !
    integer(kint) :: i, j
    
    write(fileunit,'(A)') '<?xml version="1.0"?>'
    write(fileunit,'(A)') '<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian">'
    write(fileunit,'(A)') ' <PolyData>'
    write(fileunit,'(A,I8,A)',advance='no') '   <Piece NumberOfPoints="', self%particles%N, &
        '" NumberOfVerts="0" NumberOfLines="0" '
    if (self%faceKind == TRI_PANEL) then
        write(fileunit,'(A,I8,A)') 'NumberOfStrips="0" NumberOfPolys="', 3*self%faces%N_Active, '">'
    elseif(self%faceKind == QUAD_PANEL) then
        write(fileunit,'(A,I8,A)') 'NumberOfStrips="0" NumberOfPolys="', 4*self%faces%N_Active, '">'
    elseif(self%faceKind == QUAD_CUBIC_PANEL) then
        write(fileunit,'(A,I8,A)') 'NumberOfStrips="0" NumberOfPolys="', 9*self%faces%N_Active, '">'
    endif
    
    write(fileunit,'(A)') '    <Points>'
    if (self%geomKind == PLANAR_GEOM) then
        write(fileunit,'(A)') '      <DataArray type="Float32" NumberOfComponents="3" format="ascii">'
        do i=1, self%particles%N
            !write(fileunit,'(3(F24.9,A))',advance='no') self%particles%x(i), ' ', self%particles%y(i), ' ', dzero, ' '
            write(fileunit,'(3(F24.9,A))') self%particles%x(i), ' ', self%particles%y(i), ' ', dzero, ' '
        enddo
    elseif (self%geomKind == SPHERE_GEOM) then
        write(fileunit,'(A)') '      <DataArra type="Float32" NumberOfComponents="3" format="ascii">'
        do i=1, self%particles%N
            !write(fileunit,'(3(F24.9,A))',advance='no') self%particles%x(i), ' ', self%particles%y(i), ' ', &
            !    self%particles%z(i), ' '
            write(fileunit,'(3(F24.9,A))') self%particles%x(i), ' ', self%particles%y(i), ' ', &
                self%particles%z(i), ' '
        enddo
    endif
    !write(fileunit,'(A)') ''
    write(fileunit,'(A)') '      </DataArray>'
    write(fileunit,'(A)') '    </Points>'
    write(fileunit,'(A)') '    <Polys>'
    write(fileunit,'(A)') '      <DataArray type="Int32" Name="connectivity" format="ascii">'
    if (self%faceKind == TRI_PANEL) then
        do i=1, self%faces%N
            if (.not. self%faces%hasChildren(i)) then
                do j=1,3
                    write(fileunit,'(3(I20,A))') self%faces%vertices(j,i)-1, ' ', &
                        self%faces%vertices(mod(j,3)+1,i)-1, ' ', self%faces%centerParticles(1,i)-1, ' '
                enddo
            endif
        enddo
        !write(fileunit,'(A)') ''
        write(fileunit,'(A)') '      </DataArray>'
        write(fileunit,'(A)') '      <DataArray type="Int32" Name="offsets" format="ascii">'
        do i = 1, 3*self%faces%N_Active
            write(fileunit,'(I20)') 3*i
        enddo
    elseif (self%faceKind == QUAD_PANEL) then
        do i=1, self%faces%N
            if (.not. self%faces%hasChildren(i)) then
                do j=1,4
                    write(fileunit,'(4(I20,A))') self%faces%vertices(j,i)-1, ' ', &
                        self%faces%vertices(mod(j,4)+1,i)-1, ' ', self%faces%centerParticles(1,i)-1, ' '
                enddo
            endif
        enddo
        !write(fileunit,'(A)') ''
        write(fileunit,'(A)') '      </DataArray>'
        write(fileunit,'(A)') '      <DataArray type="Int32" Name="offsets" format="ascii">'
        do i = 1, 4*self%faces%N_Active
            write(fileunit,'(I20)') 3*i
        enddo
    elseif (self%faceKind == QUAD_CUBIC_PANEL) then
        do i=1,self%faces%N 
            if (.not. self%faces%hasChildren(i)) then
                write(fileunit,'(4(I20,A))') self%faces%vertices(1,i)-1, ' ', self%faces%vertices(2,i)-1, ' ', &
                    self%faces%centerParticles(1,i)-1, ' ', self%faces%vertices(12,i)-1, ' '
                write(fileunit,'(4(I20,A))') self%faces%vertices(2,i)-1, ' ', self%faces%vertices(3,i)-1, ' ', &
                    self%faces%centerParticles(2,i)-1, ' ', self%faces%centerParticles(1,i)-1, ' '
                write(fileunit,'(4(I20,A))') self%faces%vertices(3,i)-1, ' ', self%faces%vertices(4,i)-1, ' ', &
                    self%faces%vertices(5,i)-1, ' ', self%faces%centerParticles(2,i)-1, ' '
                write(fileunit,'(4(I20,A))') self%faces%centerParticles(2,i)-1, ' ', &
                    self%faces%vertices(5,i)-1, ' ', self%faces%vertices(6,i)-1, ' ', self%faces%centerParticles(3,i)-1, ' '
                write(fileunit,'(4(I20,A))') self%faces%centerParticles(3,i)-1, ' ', &
                    self%faces%vertices(6,i)-1, ' ', self%faces%vertices(7,i)-1, ' ', self%faces%vertices(8,i)-1, ' '
                write(fileunit,'(4(I20,A))') self%faces%centerParticles(4,i)-1, ' ', &
                    self%faces%centerParticles(3,i)-1, ' ', self%faces%vertices(8,i)-1, ' ', self%faces%vertices(9,i)-1, ' '
                write(fileunit,'(4(I20,A))') self%faces%vertices(11,i)-1, ' ', &
                    self%faces%centerParticles(4,i)-1, ' ', self%faces%vertices(9,i)-1, ' ', self%faces%vertices(10,i)-1, ' '
                write(fileunit,'(4(I20,A))') self%faces%vertices(12,i)-1, ' ', &
                    self%faces%centerParticles(1,i)-1, ' ', self%faces%centerParticles(4,i)-1, ' ', self%faces%vertices(11,i)-1, ' '
                do j=1,4
                    write(fileunit,'(I20,A)') self%faces%centerParticles(j,i)-1, ' '
                enddo
            endif
        enddo
        !write(fileunit,'(A)') ''
        write(fileunit,'(A)') '      </DataArray>'
        write(fileunit,'(A)') '      <DataArray type="Int32" Name="offsets" format="ascii">'
        do i = 1, 9*self%faces%N_Active
            write(fileunit,'(I20,A)') 4*i
        enddo
    endif
    write(fileunit,'(A)') ''
    write(fileunit,'(A)') '      </DataArray>'
    write(fileunit,'(A)') '    </Polys>'
end subroutine

subroutine writeVTKSerialEndXML(self, fileunit)
    class(PolyMesh2d), intent(in) :: self
    integer(kint), intent(in) :: fileunit
    write(fileunit,'(A)') '    </Piece>'
    write(fileunit,'(A)') '  </PolyData>'
    write(fileunit,'(A)') '</VTKFile>'
end subroutine

pure function nVerticesInMesh(self, nestLevel)
    integer(kint) :: nVerticesInMesh
    class(PolyMesh2d), intent(in) :: self
    integer(kint), intent(in) :: nestLevel
    !
    integer(kint) :: i
    nVerticesInMesh = 0
    select case (self%meshSeed)
        case (TRI_HEX_SEED)
            do i = 2**nestLevel +1, 2**(nestLevel+1)
                nVerticesInMesh = nVerticesInMesh + i
            enddo
            nVerticesInMesh = 2*nVerticesInMesh + 2**(nestLevel+1) + 1
        case (QUAD_RECT_SEED)
            nVerticesInMesh = 3
            do i=1,nestLevel
                nVerticesInMesh = nVerticesInMesh + 2**i
            enddo
            nVerticesInMesh = nVerticesInMesh * nVerticesInMesh
        case (CUBIC_PLANE_SEED)
            nVerticesInMesh = 49
            do i=1,nestLevel
                nVerticesInMesh = (2*sqrt(real(nVerticesInMesh,kreal))-1)**2
            enddo
        case (ICOS_TRI_SPHERE_SEED)
            nVerticesInMesh = 2 + 10*4**nestLevel
        case (CUBED_SPHERE_SEED)
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
    select case (self%meshSeed)
        case (TRI_HEX_SEED)
            nFacesInMesh = 6*4**nestLevel
        case (QUAD_RECT_SEED)
            nFacesInMesh = 4*4**nestLevel
        case (CUBIC_PLANE_SEED)
            nFacesInMesh = 4*4**nestLevel
        case (ICOS_TRI_SPHERE_SEED)
            nFacesInMesh = 20*4**nestLevel
        case (CUBED_SPHERE_SEED)
            nFacesInMesh = 6*4**nestLevel
    end select
end function

pure function nEdgesInMesh(self, nVerts, nFaces, nestLevel)
    integer(kint) :: nEdgesInMesh
    class(PolyMesh2d), intent(in) :: self
    integer(kint), intent(in) :: nVerts, nFaces
    integer(kint), intent(in), optional :: nestLevel
    !
    integer(kint) :: nv, nf, i

    nEdgesInMesh = 0
    select case (self%meshSeed)
        case (QUAD_RECT_SEED)
            nEdgesInMesh = nFaces + nVerts -1
        case (TRI_HEX_SEED)
            nEdgesInMesh = nFaces + nVerts -1
        case (CUBED_SPHERE_SEED)
            nEdgesInMesh = nFaces + nVerts -2
        case (ICOS_TRI_SPHERE_SEED)
            nEdgesInMesh = nFaces + nVerts -2
        case (CUBIC_PLANE_SEED)
            nf = 4*4**nestLevel
            nv = 3
            do i=1,nestLevel
                nv = nv + 2**i
            enddo
            nv = nv*nv
            nEdgesInMesh = nf + nv - 1
    end select
end function

subroutine getSeed(self)
    class(PolyMesh2d), intent(inout) :: self
    !
    character(len=128) :: seedfilename
    integer(kint) :: i, j, nSeedParticles, nSeedEdges, nSeedFaces, nSeedVerts, nv, nc
    integer(kint), allocatable :: seedEdgeOrigs(:), seedEdgeDests(:), seedEdgeLefts(:), seedEdgeRights(:)
    integer(kint), allocatable :: seedEdgeInts(:,:), seedFaceEdges(:,:), seedFaceVerts(:,:), seedFaceCenters(:,:)
    real(kreal), allocatable :: seedXYZ(:,:)
    real(kreal) :: jac(16)
    integer(kint) :: nMaxParticles, nMaxEdges, nMaxFaces

!    call LogMessage(log, DEBUG_LOGGING_LEVEL, trim(logkey)//" entering : ", "getSeed")

    if (associated(self%particles)) deallocate(self%particles)
    if (associated(self%edges)) deallocate(self%edges)
    if (associated(self%faces)) deallocate(self%faces)

    allocate( Particles :: self%particles)

    if (self%meshSeed == TRI_HEX_SEED) then
         allocate( LinearEdges :: self%edges)
         allocate( TriLinearFaces :: self%faces)

        seedfilename = "triHex.namelist"

        nSeedParticles = 13
        nSeedEdges = 12
        nSeedFaces = 6
        nSeedVerts = 7

        self%faceKind = TRI_PANEL
        self%geomKind = PLANAR_GEOM
        nv = 3
        nc=1
    elseif (self%meshSeed == QUAD_RECT_SEED) then
        allocate(LinearEdges :: self%edges)
        allocate(QuadLinearFaces :: self%faces)

        seedfilename = "quadRect.namelist"
        nSeedParticles = 13
        nSeedEdges = 12
        nSeedFaces = 4
        nv = 4
        nc=1
        self%faceKind = QUAD_PANEL
        self%geomKind = PLANAR_GEOM
    elseif (self%meshSeed == CUBIC_PLANE_SEED) then
        allocate(CubicEdges :: self%edges)
        allocate(QuadCubicFaces :: self%faces)
        seedfilename = "quadCubic.namelist"
        nSeedParticles = 49
        nSeedEdges = 12
        nSeedFaces = 4
        nv=12
        nc=4
        self%faceKind= QUAD_CUBIC_PANEL
        self%geomKind = PLANAR_GEOM
    elseif (self%meshSeed == ICOS_TRI_SPHERE_SEED) then
        allocate(LinearEdges :: self%edges)
        allocate(TriLinearFaces :: self%faces)

        seedfilename = "icosTri.namelist"
        nSeedParticles = 32
        nSeedEdges = 30
        nSeedFaces = 20
        nv=3
        nc=1
        self%faceKind = TRI_PANEL
        self%geomKind = SPHERE_GEOM
    elseif (self%meshSeed == CUBED_SPHERE_SEED) then
        allocate(LinearEdges :: self%edges)
        allocate(QuadLinearFaces :: self%faces)

        seedfilename = "cubedSphere.namelist"
        nSeedParticles = 14
        nSeedEdges = 12
        nSeedFaces = 6
        nv = 4
        nc=1
        self%faceKind = QUAD_PANEL
        self%geomKind = SPHERE_GEOM
!    elseif (self%mesh_type == "cubedSphereCubic") then
!        allocate(CubicEdges :: self%edges)
!        allocate(QuadCubicFaces :: self%faces)
!
!        self%faceKind = QUAD_CUBIC_PANEL
!        self%geomKind = SPHERE_GEOM
!        nv = 12
!        nc = 4
    endif

    if (self%meshSeed == CUBIC_PLANE_SEED) then
        nMaxParticles = self%nVerticesInMesh(self%maxNest)
    else
        nMaxParticles = self%nVerticesInMesh(self%maxNest) + self%nFacesInMesh(self%maxNest)
    endif
    nMaxFaces = 0
    nMaxEdges = 0
    do i=0, self%maxNest
        nMaxFaces = nMaxFaces + self%nFacesInMesh(i)
        nMaxEdges = nMaxEdges + self%nEdgesInMesh(self%nVerticesInMesh(i), self%nFacesInMesh(i),i)
    enddo


!    print *, "nMaxParticles = ", nMaxParticles, ", nMaxFaces = ", nMaxFaces, ", nMaxEdges = ", nMaxEdges

    call self%particles%init(nMaxParticles, self%geomKind)
    call self%edges%init(nMaxEdges)
    call self%faces%init(self%faceKind, nMaxFaces)

!    call LogMessage(log,DEBUG_LOGGING_LEVEL,trim(logKey)//" return from ", "initial allocations.")
!    call self%logStats(log)

    allocate(seedXYZ(3,nSeedParticles))
    allocate(seedEdgeOrigs(nSeedEdges))
    allocate(seedEdgeDests(nSeedEdges))
    allocate(seedEdgeLefts(nSeedEdges))
    allocate(seedEdgeRights(nSeedEdges))
    if (self%faceKind == TRI_PANEL) then
        allocate(seedFaceVerts(3,nSeedFaces))
        allocate(seedFaceEdges(3,nSeedFaces))
        allocate(seedFaceCenters(1,nSeedFaces))
        allocate(seedEdgeInts(1,1))
    elseif (self%faceKind == QUAD_PANEL) then
        allocate(seedFaceVerts(4,nSeedFaces))
        allocate(seedFaceEdges(4,nSeedFaces))
        allocate(seedFaceCenters(1,nSeedFaces))
        allocate(seedEdgeInts(1,1))
    elseif (self%faceKind == QUAD_CUBIC_PANEL) then
        allocate(seedFaceVerts(12,nSeedFaces))
        allocate(seedFaceEdges(4,nSeedFaces))
        allocate(seedFaceCenters(4,nSeedFaces))
        allocate(seedEdgeInts(2,nSeedEdges))
    endif

    call self%readSeedFile(seedfilename, seedXYZ, seedEdgeOrigs, seedEdgeDests, seedEdgeLefts, seedEdgeRights, &
        seedEdgeInts, seedFaceVerts, seedFaceEdges, seedFaceCenters)

    seedXYZ = self%maxRadius * seedXYZ

!    call LogMessage(log, DEBUG_LOGGING_LEVEL, trim(logkey)//" returned from ", "reedSeedFile")

    ! initialize mesh
    do i=1, nSeedParticles
        call self%particles%insert(seedXYZ(:,i), seedXYZ(:,i))
    enddo
!    call LogMessage(log, DEBUG_LOGGING_LEVEL, trim(logkey)//" returned from ", "insert particles")

    if ( self%particles%N /= nSeedParticles ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL, logkey//" initMeshFromSeed ERROR : "," particles%N.")
	endif

	do i=1, nSeedEdges
	    call self%edges%insert(seedEdgeOrigs(i), seedEdgeDests(i), seedEdgeLefts(i), seedEdgeRights(i), seedEdgeInts(:,i))
	enddo
!	call LogMessage(log, DEBUG_LOGGING_LEVEL, trim(logkey)//" returned from ", "insert edges")

	if ( self%edges%N /= nSeedEdges ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL, logkey//" initMeshFromSeed ERROR : "," edges%N.")
	endif

	do i=1, nSeedFaces
	    call self%faces%insert(seedFaceCenters(:,i), seedFaceVerts(:,i), seedFaceEdges(:,i))
	enddo
!    call LogMessage(log, DEBUG_LOGGING_LEVEL, trim(logkey)//" returned from ", "insert faces")
!    print *, self%faces%vertices(:,1:4)
!    call LogMessage(log, DEBUG_LOGGING_LEVEL, trim(logkey)//" returned from ", "insert faces2")
	!
	!   set face areas and particle weights
	!
	self%faces%N_Active = nSeedFaces
	self%faces%area = dzero
    do i=1, nSeedFaces
        self%faces%area(i) = self%faces%setArea(i, self%particles)
    enddo

!    call LogMessage(log, DEBUG_LOGGING_LEVEL, trim(logkey)//" returned from ", " set face area.")

    if (self%faceKind == QUAD_PANEL .or. self%faceKind == TRI_PANEL) then
        do i=1, nSeedFaces
            self%particles%weight(self%faces%centerParticles(1,i)) = self%faces%area(i)
        enddo
    elseif (self%faceKind == QUAD_CUBIC_PANEL) then
        do i=1, nSeedFaces
            seedXYZ(:,1) = self%particles%physCoord(self%faces%vertices(1,i))
            seedXYZ(:,2) = self%particles%physCoord(self%faces%vertices(4,i))
            seedXYZ(:,3) = self%particles%physCoord(self%faces%vertices(7,i))
            seedXYZ(:,4) = self%particles%physCoord(self%faces%vertices(10,i))
            do j=1,12
                jac(j) = bilinearPlaneJacobian(seedXYZ(:,1:4), quad16_vertex_qp(1,j), quad16_vertex_qp(2,j))
            enddo
            do j=1,4
                jac(12+j) = bilinearPlaneJacobian(seedXYZ(:,1:4), quad16_center_qp(1,j), quad16_center_qp(2,j))
            enddo
!            do j=1,16
!                print *, "jac(",j,") = ", jac(j)
!            enddo
            do j=1,12
                self%particles%weight(self%faces%vertices(j,i)) = quad16_vertex_qw(j) * jac(j)
            enddo
            do j=1,4
                self%particles%weight(self%faces%centerParticles(j,i)) = quad16_center_qw(j)*jac(12+j)
            enddo
        enddo
    endif
end subroutine

subroutine readSeedFile(self, seedfilename, xyz, origs, dests, lefts, rights, ints, faceverts, faceedges, facecenters)
    class(PolyMesh2d), intent(in) :: self
    character(len=*), intent(in) :: seedfilename
    real(kreal), dimension(:,:), intent(inout) :: xyz
    integer(kint), dimension(:), intent(inout) :: origs, dests, lefts, rights
    integer(kint), dimension(:,:), intent(inout) :: ints
    integer(kint), dimension(:,:), intent(inout) :: faceverts, faceedges, facecenters
    !
    integer(kint) :: readStat
    integer(kint) :: nParticles, nEdges, nFaces, nVertsPerFace, nCentersPerFace, nIntsPerEdge
    namelist /seed/ xyz, origs, dests, lefts, rights, ints, faceverts, faceedges, facecenters

!    call LogMessage(log, DEBUG_LOGGING_LEVEL,trim(logkey)//" entering : ", "readSeedFile")

!    print *, "nParticles = ", size(xyz,1)
!    print *, "size(origs) = ", size(origs)
!    print *, "size(dests) = ", size(dests)
!    print *, "size(faceVerts,1)", size(faceverts,1)
!    print *, "size(faceverts,2)", size(faceverts,2)
!    print *, "size(ints)", size(ints,1), size(ints,2)

    xyz = dzero
    origs = 0
    dests = 0
    lefts = 0
    rights = 0
    ints = 0
    faceverts = 0
    faceedges = 0
    facecenters = 0
    open(unit=READ_UNIT, file=seedfilename, action='read', status='old', iostat=readStat)
    if (readStat /= 0) then
        call LogMessage(log, ERROR_LOGGING_LEVEL, "readSeedFile error : cannot open file ", seedfilename)
        return
    endif
    read(READ_UNIT, nml=seed, iostat=readstat)
!    print *, "readstat ", readstat
!    print *, faceVerts
    origs = origs+1
    dests = dests+1
    lefts = lefts+1
    rights= rights+1
    ints = ints+1
    faceverts = faceverts+1
    faceedges = faceedges+1
    facecenters = facecenters+1
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
