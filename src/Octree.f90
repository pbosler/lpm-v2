module OctreeModule

use NumberKindsModule
use UtilitiesModule
use SphereGeomModule
use LoggerModule
use Box3dModule
use ParticlesOOModule
use FacesOOModule

#include "LpmConfig.h"

implicit none
private

public :: Octree, TreeNode

type :: TreeNode
    class(Box3d) :: box
    integer(kint) :: level
    class(TreeNode), pointer :: parent => null()
    class(TreeNode), pointer, dimension(8) ::  kids => null()
    integer(kint) :: nKids = 0
    
    integer(kint), dimension(:), allocatable :: indices
    integer(kint) :: nInds = 0
    
    contains
        procedure :: init
        procedure :: hasKids
        procedure :: isLeaf
        procedure :: shrinkBoxParticles
        procedure :: shrinkBoxFaces
        generic :: shrinkBox => shrinkBoxParticles, shrinkBoxFaces
        final :: deleteNode
end type

type :: Octree
    class(TreeNode) :: root
    class(Particles), pointer :: particles => null()
    class(Faces), pointer :: faces => null()
    
    integer(kint) :: nnodes
    
    
    contains
        procedure :: initParticles
        procedure :: initFaces
        generic :: init => initParticles, initFaces
        procedure :: buildTreeMaxDepth
        procedure :: buildTreeMaxCoordsPerNode
        final :: deleteTree
end type

contains

subroutine init(node, box, parent, inds)
    class(TreeNode), intent(inout) :: node
    class(Box3d), intent(in) :: box
    class(TreeNode), pointer, intent(in) :: parent
    integer(kint), dimension(:), intent(in) :: inds
    integer(kint) :: i
    
    call node%box%init(box%xmin, box%xmax, box%ymin, box%ymax, box%zmin, box%zmax)
    if (associated(parent)) then
        node%parent => parent
        node%level = parent%level + 1
    else
        nullify(node%parent)
        node%level = 0
    endif
    
    do i=1,8
        nullify(node%kids(i))
    enddo
    
    allocate(node%indices(size(inds)))
    node%indices = inds
end subroutine

recursive subroutine deleteNode(node)
    type(TreeNode), intent(inout) :: node
    integer(kint) :: i
    deallocate(node%indices)
    do i = 1, node%nKids
        call deleteNode(node%kids(i))
        deallocate(node%kids(i))
    enddo
end subroutine

pure function hasKids(node)
    logical(klog) :: hasKids
    class(TreeNode), intent(in) :: node
    hasKids = (node%nKids > 0)
end function

pure function isLeaf(node)
    logical(klog) :: isLeaf
    class(TreeNode), intent(in) :: node
    isLeaf = (node%nKids == 0)
end function

subroutine shrinkBoxFaces(node, aFaces, aParticles)
    class(TreeNode), intent(inout) :: node
    class(Faces), intent(in) :: aFaces
    class(Particles), intent(inout) :: aParticles
    integer(kint) :: i
    real(kreal) :: xmin, xmax, ymin, ymax, zmin, zmax
    real(kreal), dimension(3) :: cntd
    cntd = aFaces%physCentroid(node%indices(1), aParticles)
    xmin = cntd(1)
    xmax = cntd(1)
    ymin = cntd(2)
    ymax = cntd(2)
    zmin = cntd(3)
    zmax = cntd(3)
    do i = 2, node%nInds
        cntd = aFaces%physCentroid(node%indices(i), aParticles)
        if (cntd(1) < xmin) xmin = cntd(1)
        if (cntd(1) > xmax) xmax = cntd(1)
        if (cntd(2) < ymin) ymin = cntd(2)
        if (cntd(2) > ymax) ymax = cntd(2)
        if (cntd(3) < zmin) zmin = cntd(3)
        if (cntd(3) > zmax) zmax = cntd(3)
    enddo
    call node%box%init(xmin, xmax, ymin, ymax, zmin, zmax)
end subroutine

subroutine shrinkBoxParticles(node, aParticles)
    class(TreeNode), intent(inout)
    class(Particles), intent(in) :: aParticles
    integer(kint) :: i
    real(kreal) :: xmin, xmax, ymin, ymax, zmin, zmax
    real(kreal), dimension(3) :: pos
    pos = aParticles%physCoord(node%indices(1))
    xmin = pos(1)
    xmax = pos(1)
    ymin = pos(2)
    ymax = pos(2)
    zmin = pos(3)
    zmax = pos(3)
    do i = 2, node%nInds
        pos = aParticles%physCoord(node%indices(i))
        if (pos(1) < xmin) xmin = pos(1)
        if (pos(1) > xmax) xmax = pos(1)
        if (pos(2) < ymin) ymin = pos(2)
        if (pos(2) > ymax) ymax = pos(2)
        if (pos(3) < zmin) zmin = pos(3)
        if (pos(3) > zmax) zmax = pos(3)
    enddo
    call node%box%init(xmin, xmax, ymin, ymax, zmin, zmax)
end subroutine

subroutine initFaces(tree, radius, aParticles, aFaces)
    class(Octree), intent(inout) :: tree
    real(kreal), intent(in) :: radius
    class(Particles), intent(in) :: aParticles
    class(Faces), intent(in) :: aFaces
    !
    class(Box3d) :: rbox
    integer(kint), dimension(:), allocatable :: rinds
    integer(kint) :: i, ins
    
    tree%particles => aParticles
    tree%faces => aFaces
    
    allocate(rinds(aFaces%N_Active))
    ins = 1
    do i = 1, aFaces%N
        if (.not. aFaces%hasChildren(i)) then
            rinds(ins) = i
            ins = ins+1
        endif
    enddo
    
    call rbox%init(radius)
    call tree%root%init(rbox, null(), rinds)
    tree%nnodes = 1
    
    deallocate(rinds)
end subroutine

subroutine initParticles(tree, radius, aParticles)
    class(Octree), intent(inout) :: tree
    real(kreal), intent(in) :: radius
    class(Particles), intent(in) :: aParticles
    !
    class(Box3d) :: rbox
    integer(kint), dimension(:), allocatable :: rinds
    integer(kint) :: i
    
    tree%particles => aParticles
    nullify(tree%faces)
    
    call rbox%init(radius)
    allocate(rinds(aParticles%N))
    do i=1, aParticles%n
        rinds(i) = i
    enddo
    
    call tree%root%init(rbox, null(), rinds)
    tree%nnodes = 1
    
    deallocate(rinds)
end subroutine

subroutine deleteTree(tree)
    type(Octree), intent(inout) :: tree
    call deleteNode(tree%root)
    nullify(tree%particles)
    nullify(tree%faces)
end subroutine

recursive subroutine buildTreeMaxDepth(tree, node, maxLevels)
    class(Octree), intent(inout) :: tree
    class(TreeNode), intent(inout) :: node
    integer(kint), intent(in) :: maxLevels
    !
    class(Box3d), dimension(8) :: kidboxes
    integer(kint), allocatable, dimension(:,:) :: kidIndices
    integer(kint), dimension(8) :: nKidInds
    integer(kint) :: i, j, ins
    
end subroutine

recursive subroutine buildTreeMaxCoordsPerNode(tree, node, maxN)
    class(Octree), intent(inout) :: tree
    class(TreeNode), intent(inout) :: node
    integer(kint), intent(in) :: maxN
    !
    class(Box3d), dimension(8) :: kidboxes
    integer(kint), allocatable, dimension(:,:) :: kidIndices
    integer(kint), dimension(8) :: nKidInds
    integer(kint) :: i, j, ins
    real(kreal), dimension(3) :: pos
    
    if (node%nInds <= maxN) then
        return
    else
        kidboxes = node%box%bisectAll()
        allocate(kidIndices(node%nInds, 8))
        kidIndices = 0
        nKidInds = 0
        if (associated(tree%faces)) then
            do i = 1, node%nInds
                pos = tree%faces%physCentroid(node%indices(i), tree%particles)
                do j = 1,8
                    if (kidboxes(j)%containsPoint(pos)) then
                        nKidInds(j) = nKidInds(j) + 1
                        kidIndices(nKidInds(j), j) = node%indices(i)
                    endif
                enddo
            enddo
        else
            do i = 1, node%nInds
                pos = tree%particles%physCoord(node%indices(i))
                do j = 1,8
                    if (kidboxes(j)%containsPoint(pos)) then
                        nKidInds(j) = nKidInds(j) + 1
                        kidIndices(nKidInds(j), j) = node%indices(i)
                    endif
                enddo
            enddo
        endif
        ins = 1
        do j = 1,8
            if (nKidInds(j) > 0) then
                allocate(node%kids(ins))
                call node%kids(ins)%init(kidboxes(j), node, kidIndices(1:nKidInds(j), j))
                node%nKids = node%nKids + 1
                ins = ins+1
            endif
        enddo
        deallocate(kidIndices)
    endif
end subroutine

end module