module Box3dModule

use NumberKindsModule
use UtilitiesModule
use SphereGeomModule
use LoggerModule

#include "LpmConfig.h"

implicit none
private

public Box3d
public BOX_PAD

type :: Box3d
   real(kreal) :: xmin, xmax
   real(kreal) :: ymin, ymax
   real(kreal) :: zmin, zmax
   real(kreal) :: min_radius, max_radius
   real(kreal), dimension(3) :: centroid
   real(kreal) :: volume
   real(kreal) :: aspectRatio

   contains
     procedure :: initSpecific
     procedure :: initPadded
     generic :: init => initSpecific, initPadded
     procedure :: containsPoint
     procedure :: bisectAll
     !procedure :: bisectAlongDims
     procedure :: logStats
     procedure :: minEdgeLength
     procedure :: maxEdgeLength
     procedure, private :: calc_min_radius
     procedure, private :: calc_max_radius
     procedure :: corners
     procedure :: faceCentroids
end type

real(kreal), parameter :: BOX_PAD = 0.00001

contains

pure subroutine initPadded(self, radius)
    class(Box3d), intent(inout) ::  self
    real(kreal), intent(in) :: radius

    self%xmin = -radius - BOX_PAD
    self%xmax = radius + BOX_PAD
    self%ymin = -radius - BOX_PAD
    self%ymax = radius + BOX_PAD
    self%zmin = -radius - BOX_PAD
    self%zmax = radius + BOX_PAD

    self%volume = (2.0_kreal * (radius + BOX_PAD))**3
    self%centroid = 0.0_kreal
    self%aspectRatio = 1.0_kreal
    self%min_radius = self%calc_min_radius()
    self%max_radius = self%calc_max_radius()
end subroutine

pure subroutine initSpecific(self, xmin, xmax, ymin, ymax, zmin, zmax)
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

subroutine logStats(self, aLog)
    class(Box3d), intent(in) :: self
    type(Logger), intent(inout) :: aLog
    character(len=256) :: infostring

    call StartSection(aLog)
    write(infostring,'(6(A,G10.4),A)') 'x in [', self%xmin, ', ', self%xmax, &
        ']; y in [', self%ymin, ', ', self%ymax, ']; z in [', self%zmin, ', ', self%zmax, ']'
    call LogMessage(aLog, TRACE_LOGGING_LEVEL, "box info: ", trim(infostring))
    call LogMessage(aLog, TRACE_LOGGING_LEVEL, "volume = ", self%volume)
    call LogMessage(aLog, TRACE_LOGGING_LEVEL, "aspect ratio = ", self%aspectRatio)
    call LogMessage(aLog, TRACE_LOGGING_LEVEL, "min_radius = ", self%min_radius)
    call LogMessage(aLog, TRACE_LOGGING_LEVEL, "max_radius = ", self%max_radius)
    call LogMessage(aLog, TRACE_LOGGING_LEVEL, "centroid = ", self%centroid)
    call EndSection(aLog)
end subroutine

end module
