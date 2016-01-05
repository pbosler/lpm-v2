module PlaneGeomModule
!------------------------------------------------------------------------------
! Lagrangian Particle / Panel Method
!------------------------------------------------------------------------------
!
!> @author
!> Peter Bosler, Department of Mathematics, University of Michigan
!
!> @defgroup PlaneGeom Planar Geometry
!> defines functions for performing geometric calculations in the plane
!
!
! DESCRIPTION:
!> @file
!> defines functions for performing geometric calculations in the plane
!
!------------------------------------------------------------------------------

use NumberKindsModule

implicit none

public

contains
!----------------
! Basic geometry : length, area, centers of mass, etc.
!----------------


!------------------------------------------------------------------------------
!> @author Peter Bosler
!> @brief Euclidean distance between two points in the plane
!
!> @ingroup PlaneGeom
!
!> @param[in] xyA double precision real size(2); vector in the plane
!> @param[in] xyB double precision real size(2); vector in the plane
!> @return Distance double precision real, scalar distance between xyA and xyB
!------------------------------------------------------------------------------
function Distance( xyA, xyB)
	real(kreal) :: Distance
	real(kreal), intent(in) :: xyA(2), xyB(2)
	Distance = sqrt( sum( (xyB - xyA)*(xyB - xyA) ) )
end function

!------------------------------------------------------------------------------
!> @author Peter Bosler
!> @brief Find midpoint of line segment connecting two points in the plane
!
!> @ingroup PlaneGeom
!
!> @param[in] xyA double precision real size(2); vector in the plane
!> @param[in] xyB double precision real size(2); vector in the plane
!> @return Midpoint double precision real size(2)
!------------------------------------------------------------------------------
function Midpoint( xyA, xyB )
	real(kreal) :: Midpoint(2)
	real(kreal), intent(in) :: xyA(2), xyB(2)
	Midpoint = 0.5_kreal*(xyA + xyB)
end function

!------------------------------------------------------------------------------
!> @author Peter Bosler
!> @brief Cetroid of a triangle created by three points in the plane
!
!> @ingroup PlaneGeom
!
!> @param[in] xyA double precision real size(2); vector in the plane
!> @param[in] xyB double precision real size(2); vector in the plane
!> @param[in] xyC double precision real size(2); vector in the plane
!> @return TriCentroid double precision real size(2)
!------------------------------------------------------------------------------
function TriCentroid( xyA, xyB, xyC )
	real(kreal) :: TriCentroid(2)
	real(kreal), intent(in) :: xyA(2), xyB(2), xyC(2)
	TriCentroid = (xyA + xyB + xyC)/3.0_kreal
end function

!------------------------------------------------------------------------------
!> @author Peter Bosler
!> @brief Cetroid of a quadrilateral created by four points in the plane
!
!> @ingroup PlaneGeom
!
!> @param[in] xyA double precision real size(2); vector in the plane
!> @param[in] xyB double precision real size(2); vector in the plane
!> @param[in] xyC double precision real size(2); vector in the plane
!> @param[in] xyD double precision real size(2); vector in the plane
!> @return QuadCentroid double precision real size(2)
!------------------------------------------------------------------------------
function QuadCentroid( xyA, xyB, xyC, xyD )
	real(kreal) :: QuadCentroid(2)
	real(kreal), intent(in) :: xyA(2), xyB(2), xyC(2), xyD(2)
	QuadCentroid = 0.25_kreal*(xyA + xyB + xyC + xyD)
end function

!------------------------------------------------------------------------------
!> @author Peter Bosler
!> @brief area of a planar triangle whose vertices are the input vectors
!
!> @ingroup PlaneGeom
!
!> @param[in] xyA double precision real size(2); vector in the plane
!> @param[in] xyB double precision real size(2); vector in the plane
!> @param[in] xyC double precision real size(2); vector in the plane
!> @return TriArea double precision real; scalar area of triangle 
!------------------------------------------------------------------------------
function TriArea( xA, xB, xC )
	real(kreal) :: TriArea
	real(kreal), intent(in) :: xA(2), xB(2), xC(2)
	TriArea = 0.5_kreal*abs( -xB(1)*xA(2) + xC(1)*xA(2) + xA(1)*xB(2) - xC(1)*xB(2) - xA(1)*xC(2) + xB(1)*xC(2))
end function


!----------------
! Misc. functions
!----------------
!

end module


