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
     procedure :: bisectAlongDims
     procedure :: logStats
     procedure :: minEdgeLength
     procedure :: maxEdgeLength
     procedure, private :: calc_min_radius
     procedure, private :: calc_max_radius
     procedure :: corners
     procedure :: faceCentroids
end type

contains

subroutine init(self, xmin, xmax, ymin, ymax, zmin, zmax)
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
  lens(1) = self.xmax - self.xmin
  lens(2) = self.ymax - self.ymin
  lens(3) = self.zmax - self.zmin
  minEdgeLength = minval(lens)
end function

pure function maxEdgeLength(self)
  real(kreal) :: maxEdgeLength
  class(Box3d), intent(in) :: self
  real(kreal), dimension(3) :: lens
  lens(1) = xmax-ymin
  lens(2) = ymax-ymin
  lens(3) = zmax-zmin
  maxEdgeLength = maxval(lens)
end function 

pure function calc_min_radius(self)
  class(Box3d), intent(in) :: self
  real(kreal) :: calc_min_radius
  real(kreal), dimension(3, 8) :: face_locs
  
  face_locs = self%faceCentroids()
  
end function




