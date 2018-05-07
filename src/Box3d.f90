module Box3dModule

use NumberKindsModule
use UtilitiesModule
use SphereGeomModule

implicit none
private

public Box3d

type :: Box3d
   real(kreal) :: xmin, xmax
   real(kreal) :: ymin, ymax
   real(kreal) :: zmin, zmax
   real(kreal) :: min_radius, max_radius
   real(kreal), dimension(3) :: centroid
   real(kreal) :: volume
   real(kreal) :: aspectRatio
   
   contains
     procedure :: init
     procedure :: containsPoint
     procedure :: bisectAll
     !procedure :: bisectAlongDims
     !procedure :: logStats
     procedure :: minEdgeLength
     procedure :: maxEdgeLength
     procedure, private :: calc_min_radius
     procedure, private :: calc_max_radius
     procedure :: corners
     procedure :: faceCentroids
end type

contains

pure subroutine init(self, xmin, xmax, ymin, ymax, zmin, zmax)
  class(Box3d), intent(inout) :: self
  real(kreal), intent(in) :: xmin, xmax, ymin, ymax, zmin, zmax
  
  self%xmin = xmin
  self%xmax = xmax
  self%ymin = ymin
  self%ymax = ymax
  self%zmin = zmin
  self%zmax = zmax

  self%volume = (xmax-xmin)*(ymax-ymin)*(zmax-zmin)
  self%centroid = 0.5_kreal * [xmin+xmax, ymin+ymax, zmin+zmax]
  self%aspectRatio = self%maxEdgeLength() / self%minEdgeLength()
  self%min_radius = self%calc_min_radius()
  self%max_radius = self%calc_max_radius()
end subroutine

pure function minEdgeLength(self)
  real(kreal) :: minEdgeLength
  class(Box3d), intent(in) :: self
  real(kreal), dimension(3) :: lens
  lens(1) = self%xmax - self%xmin
  lens(2) = self%ymax - self%ymin
  lens(3) = self%zmax - self%zmin
  minEdgeLength = minval(lens)
end function

pure function maxEdgeLength(self)
  real(kreal) :: maxEdgeLength
  class(Box3d), intent(in) :: self
  real(kreal), dimension(3) :: lens
  lens(1) = self%xmax-self%ymin
  lens(2) = self%ymax-self%ymin
  lens(3) = self%zmax-self%zmin
  maxEdgeLength = maxval(lens)
end function 

pure function calc_min_radius(self)
  class(Box3d), intent(in) :: self
  real(kreal) :: calc_min_radius
  real(kreal), dimension(3, 6) :: face_locs
  integer(kint) :: i
  real(kreal) :: test_dist
  
  face_locs = self%faceCentroids()
  calc_min_radius = huge(1.0_kreal)
  do i=1, 6
    test_dist = ChordDistance(self%centroid, face_locs(:,i))
    if (test_dist < calc_min_radius) then
        calc_min_radius = test_dist
    endif
  enddo
end function

pure function calc_max_radius(self)
    real(kreal) :: calc_max_radius
    class(Box3d), intent(in) :: self
    real(kreal), dimension(3,8) :: corners
    integer(kint) :: i
    real(kreal) :: test_dist
    
    corners = self%corners()
    calc_max_radius = 0.0_kreal
    do i=1,8
        test_dist = ChordDistance(self%centroid, corners(:,i))
        if (test_dist > calc_max_radius) then
            calc_max_radius = test_dist
        endif
    enddo
end function

pure function faceCentroids(self)
    real(kreal), dimension(3,6) :: faceCentroids
    class(Box3d), intent(in) :: self
    real(kreal) :: midx, midy, midz
    midx = 0.5_kreal*(self%xmin + self%xmax)
    midy = 0.5_kreal*(self%ymin + self%ymax)
    midz = 0.5_kreal*(self%zmin + self%zmax)
    faceCentroids(:,1) = [midx, self%ymin, midz]
    faceCentroids(:,2) = [self%xmax, midy, midz]
    faceCentroids(:,3) = [midx, self%ymax, midz]
    faceCentroids(:,4) = [self%xmin, midy, midz]
    faceCentroids(:,5) = [midx, midy, self%zmin]
    faceCentroids(:,6) = [midx, midy, self%zmax]
end function 

pure function corners(self)
    real(kreal), dimension(3,8) :: corners
    class(Box3d), intent(in) :: self
    corners(:,1) = [self%xmin, self%ymin, self%zmin]
    corners(:,2) = [self%xmin, self%ymax, self%zmin]
    corners(:,3) = [self%xmin, self%ymin, self%zmax]
    corners(:,4) = [self%xmin, self%ymax, self%zmax]
    corners(:,5) = [self%xmax, self%ymin, self%zmin]
    corners(:,6) = [self%xmax, self%ymax, self%zmin]
    corners(:,7) = [self%xmax, self%ymin, self%zmax]
    corners(:,8) = [self%xmax, self%ymax, self%zmax]
end function

pure function containsPoint(self, queryPt)
    logical(klog) :: containsPoint
    class(Box3d), intent(in) :: self
    real(kreal), dimension(3), intent(in) :: queryPt
    containsPoint = (self%xmin <= queryPt(1) .and. queryPt(1) <= self%xmax) .and. &
                    (self%ymin <= queryPt(2) .and. queryPt(2) <= self%ymax) .and. &
                    (self%zmin <= queryPt(3) .and. queryPt(3) <= self%zmax) 
end function 

pure function bisectAll(self)
    type(Box3d), dimension(8) :: bisectAll
    class(Box3d), intent(in) :: self
    real(kreal) :: midx, midy, midz
    midx = 0.5_kreal*(self%xmin + self%xmax)
    midy = 0.5_kreal*(self%ymin + self%ymax)
    midz = 0.5_kreal*(self%zmin + self%zmax)
    
    call bisectAll(1)%init(self%xmin, midx, self%ymin, midy, self%zmin, midz)
    call bisectAll(2)%init(midx, self%xmax, self%ymin, midy, self%zmin, midz)
    call bisectAll(3)%init(self%xmin, midx, midy, self%ymax, self%zmin, midz)
    call bisectAll(4)%init(midx, self%xmax, midy, self%ymax, self%zmin, midz)
    call bisectAll(5)%init(self%xmin, midx, self%ymin, midy, midz, self%zmax)
    call bisectAll(6)%init(midx, self%xmax, self%ymin, midy, midz, self%zmax)
    call bisectAll(7)%init(self%xmin, midx, midy, self%ymax, midz, self%zmax)
    call bisectAll(8)%init(midx, self%xmax, midy, self%ymax, midz, self%zmax)
end function

end module