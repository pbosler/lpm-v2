module SphereGeomModule
!------------------------------------------------------------------------------
! Lagrangian Particle / Panel Method - Spherical Model
!------------------------------------------------------------------------------
!
!> @author
!> Peter Bosler, Department of Mathematics, University of Michigan
!
!> @defgroup SphereGeom Spherical Geometry
!> defines functions for performing geometric calculations on the surface of a sphere.
!
!
! DESCRIPTION:
!> @file
!> defines functions for performing geometric calculations on the surface of a sphere.
!
!------------------------------------------------------------------------------

use NumberKindsModule

implicit none

public

interface SphereTriArea
	module procedure SphereTriAreaVector
	module procedure SphereTriAreaComponents
end interface

interface SphereDistance
	module procedure SphereDistanceVector
	module procedure SphereDistanceComponents
end interface

interface SphereArcLength
	module procedure SphereArcLengthVector
	module procedure SphereArcLengthComponents
end interface

interface Latitude
	module procedure LatitudeVector
	module procedure LatitudeComponents
end interface

interface Longitude
	module procedure LongitudeVector
	module procedure LongitudeComponents
	module procedure LongitudeComponents2
end interface

contains

!----------------
! Basic geometry : length, area, coordinates, etc.
!----------------

!------------------------------------------------------------------------------
!> @author Peter Bosler
!> @brief Euclidean distance between two points on the sphere
!
!> @ingroup SphereGeom
!
!> @param[in] xyzA double precision real size(3); vector on the sphere
!> @param[in] xyzB double precision real size(3); vector on the sphere
!> @return Distance double precision real, straight-line distance between xyzA and xyzB
pure function ChordDistance(xyzA, xyzB)
	! Outputs the Euclidean distance between two points in R3.
	! Units of length
	real(kreal), intent(in) :: xyzA(3), xyzB(3)
	real(kreal) :: ChordDistance
	ChordDistance = sqrt( sum( (xyzB - xyzA)*(xyzB-xyzA) ) )
end function


!------------------------------------------------------------------------------
!> @author Peter Bosler
!> @brief Great-circle distance between two points in the plane
!
!> @ingroup SphereGeom
!
!> @param[in] xyzA double precision real size(3); vector on the sphere
!> @param[in] xyzB double precision real size(3); vector on the sphere
!> @return Distance double precision real, great-circle distance between xyzA and xyzB
!------------------------------------------------------------------------------
pure function SphereDistanceVector( xyzA, xyzB )
! Finds the great circle distance between two points (xyzA and xyzB) on the sphere
! Units of length
	! Calling parameters
	real(KREAL), dimension(3), intent(in) :: xyzA, xyzB
	real(KREAL) :: SphereDistanceVector
	! Local variables
	real(KREAL), dimension(3) :: crossProd
	real(KREAL) :: dotProd , crossNorm

	crossProd = [xyzA(2)*xyzB(3)-xyzB(2)*xyzA(3),xyzB(1)*xyzA(3)-xyzA(1)*xyzB(3),&
    xyzA(1)*xyzB(2)-xyzB(1)*xyzA(2) ]

    crossNorm = sqrt(sum(crossProd*crossProd))

    dotProd = xyzA(1)*xyzB(1)+xyzA(2)*xyzB(2)+xyzA(3)*xyzB(3)

    SphereDistanceVector = atan2(crossNorm,dotProd)*SphereRadius
end function

!------------------------------------------------------------------------------
!> @author Peter Bosler
!> @brief Great-circle distance between two points in the plane
!
!> @ingroup SphereGeom
!
!> @param[in] xA double precision real  
!> @param[in] yA double precision real 
!> @param[in] zA
!> @param[in] xB
!> @param[in] yB
!> @param[in] zB
!> @return Distance double precision real, great-circle distance between xyzA and xyzB
!------------------------------------------------------------------------------
pure function SphereDistanceComponents(xA, yA, zA, xB, yB, zB)
! Finds the great circle distance between two points (xyzA and xyzB) on the sphere
! Units of length
	real(kreal), intent(in) :: xA, yA, zA
	real(kreal), intent(in) :: xB, yB, zB
	real(kreal) :: SphereDistanceComponents
	real(kreal) :: cp1, cp2, cp3, cpNorm, dp

	cp1 = yA*zB - yB*zA
	cp2 = xB*zA - xA*zB
	cp3 = xA*yB - xB*yA

	cpNorm = sqrt( cp1*cp1 + cp2*cp2 + cp3*cp3)

	dp = xA*xB + yA*yB + zA*zB

	SphereDistanceComponents = atan2(cpNorm,dp)*SphereRadius
end function

!------------------------------------------------------------------------------
!> @author Peter Bosler
!> @brief central angle between two points in the plane
!
!> @ingroup SphereGeom
!
!> @param[in] xyzA double precision real size(3); vector on the sphere
!> @param[in] xyzB double precision real size(3); vector on the sphere
!> @return Distance double precision real, angular separation between xyzA and xyzB
!------------------------------------------------------------------------------
pure function SphereArcLengthVector(xyzA, xyzB)
! returns the arc length between two vectors on the surface of a sphere.
! dimensionless (angle)
	! Calling parameters
	real(KREAL), dimension(3), intent(in) :: xyzA, xyzB
	real(KREAL) :: SphereArcLengthVector
	! Local variables
	real(KREAL), dimension(3) :: crossProd
	real(KREAL) :: dotProd , crossNorm

	crossProd = [xyzA(2)*xyzB(3)-xyzB(2)*xyzA(3),xyzB(1)*xyzA(3)-xyzA(1)*xyzB(3),&
    xyzA(1)*xyzB(2)-xyzB(1)*xyzA(2) ]

    crossNorm = sqrt(sum(crossProd*crossProd))

    dotProd = xyzA(1)*xyzB(1)+xyzA(2)*xyzB(2)+xyzA(3)*xyzB(3)

    SphereArcLengthVector = atan2(crossNorm,dotProd)
end function

!------------------------------------------------------------------------------
!> @author Peter Bosler
!> @brief central angle between two points in the plane
!
!> @ingroup SphereGeom
!
!> @param[in] xA double precision real  
!> @param[in] yA double precision real 
!> @param[in] zA
!> @param[in] xB
!> @param[in] yB
!> @param[in] zB
!> @return Distance double precision real, angular separation between xyzA and xyzB
!------------------------------------------------------------------------------
pure function SphereArcLengthComponents(xA, yA, zA, xB, yB, zB)
! dimensionless (angle)
	real(kreal), intent(in) :: xA, yA, zA
	real(kreal), intent(in) :: xB, yB, zB
	real(kreal) :: SphereArcLengthComponents
	real(kreal) :: cp1, cp2, cp3, cpNorm, dp

	cp1 = yA*zB - yB*zA
	cp2 = xB*zA - xA*zB
	cp3 = xA*yB - xB*yA

	cpNorm = sqrt( cp1*cp1 + cp2*cp2 + cp3*cp3)

	dp = xA*xB + yA*yB + zA*zB

	SphereArcLengthComponents = atan2(cpNorm,dp)
end function

!------------------------------------------------------------------------------
!> @author Peter Bosler
!> @brief central angle between two points in the plane
!
!> @ingroup SphereGeom
!
!> @param[in] xyz double precision real size(3); vector on the sphere
!> @return northward-pointing unit vector at position xyz
!------------------------------------------------------------------------------
pure function LatUnitVector(xyz)
!	Returns the latitudinal unit vector at the point xyz
	real(kreal) :: LatUnitVector(3)
	real(kreal), intent(in) :: xyz(3)

	if ( xyz(3) == 0.0_kreal) then
		LatUnitVector = [0.0_kreal,0.0_kreal,0.0_kreal]
	else
		LatUnitVector(1) = -xyz(1)*xyz(3)
		LatUnitVector(2) = -xyz(2)*xyz(3)
		LatUnitVector(3) = 1.0_kreal - xyz(3)*xyz(3)

		LatUnitVector = LatUnitVector / sqrt(LatUnitVector(3))
	endif
end function

!------------------------------------------------------------------------------
!> @author Peter Bosler
!> @brief central angle between two points in the plane
!
!> @ingroup SphereGeom
!
!> @param[in] xyz double precision real size(3); vector on the sphere
!> @return eastward-pointing unit vector at position xyz
!-----------------------------------------------------------------------------
pure function LonUnitVector(xyz)
!	Returns the Longitudinal unit vector at the point xyz
	real(kreal) :: LonUnitVector(3)
	real(kreal), intent(in) :: xyz(3)

	if ( xyz(3) == 0.0_kreal) then
		LonUnitVector = [0.0_kreal,0.0_kreal,0.0_kreal]
	else
		LonUnitVector(1) = -xyz(2)
		LonUnitVector(2) = -xyz(1)
		LonUnitVector(3) = 0.0_kreal

		LonUnitVector = LonUnitVector / sqrt( 1.0_kreal - xyz(3)*xyz(3))
	endif
end function

!------------------------------------------------------------------------------
!> @author Peter Bosler
!> @brief midpoint of a great-circle segment connecting two points on the sphere
!
!> @ingroup SphereGeom
!
!> @param[in] xyzA double precision real size(3); vector on the sphere
!> @param[in] xyzB double precision real size(3); vector on the sphere
!> @return Midpoint vector
!------------------------------------------------------------------------------
pure function SphereMidpoint(xyzA, xyzB)
! Finds the midpoint of two points on the sphere by finding the midpoint of the chord
! connecting the two points, then projecting the chord midpoint to the sphere.
	! Calling parameters
	real(KREAL), dimension(3), intent(in) :: xyzA, xyzB
	real(KREAL), dimension(3) :: SphereMidpoint

	SphereMidpoint = (xyzA + xyzB)/2.0_KREAL
	SphereMidpoint = SphereMidpoint/sqrt(sum(sphereMidpoint*sphereMidpoint))*SphereRadius
end function

!------------------------------------------------------------------------------
!> @author Peter Bosler
!> @brief Finds the centroid of a spherical triangle defined by three vertices
!
!> @ingroup SphereGeom
!
!> @param[in] xyzA double precision real size(3); vector on the sphere
!> @param[in] xyzB double precision real size(3); vector on the sphere
!> @param[in] xyzC double precision real size(3); vector on the sphere
!> @return centroid vector
!------------------------------------------------------------------------------
pure function SphereTriCenter(xyzA, xyzB, xyzC)
! Finds the midpoint of three points on the sphere by find their average position in Cartesian
! coordinates, then projecting that average onto the sphere.
	! Calling parameters
	real(KREAL), dimension(3), intent(in) :: xyzA, xyzB, xyzC
	real(KREAL), dimension(3) :: SphereTriCenter

	SphereTriCenter = (xyzA + xyzB + xyzC)/3.0_KREAL
	SphereTriCenter = SphereTriCenter/sqrt(sum(sphereTriCenter*sphereTriCenter))*SphereRadius
end function

!------------------------------------------------------------------------------
!> @author Peter Bosler
!> @brief Finds the centroid of a spherical quadrilateral defined by four vertices
!
!> @ingroup SphereGeom
!
!> @param[in] xyzA double precision real size(3); vector on the sphere
!> @param[in] xyzB double precision real size(3); vector on the sphere
!> @param[in] xyzC double precision real size(3); vector on the sphere
!> @param[in] xyzD double precision real size(3); vector on the sphere
!> @return centroid vector
!------------------------------------------------------------------------------
pure function SphereQuadCenter(xyzA, xyzB, xyzC, xyzD)
! Finds the midpoint of four points on the sphere by finding their average position in
! Cartesian coordinates, then projecting that average onto the sphere.
	! Calling parameters
	real(KREAL), dimension(3), intent(in) :: xyzA, xyzB, xyzC, xyzD
	real(KREAL), dimension(3) :: sphereQuadCenter

	SphereQuadCenter = (xyzA + xyzB + xyzC + xyzD)/4.0_KREAL
	SphereQuadCenter = SphereQuadCenter/sqrt(sum(sphereQuadCenter*sphereQuadCenter))*SphereRadius
end function

!------------------------------------------------------------------------------
!> @author Peter Bosler
!> @brief Finds the area of a spherical triangle defined by three vertices
!
!> @ingroup SphereGeom
!
!> @param[in] xyzA double precision real size(3); vector on the sphere
!> @param[in] xyzB double precision real size(3); vector on the sphere
!> @param[in] xyzC double precision real size(3); vector on the sphere
!> @return scalar area
!------------------------------------------------------------------------------
pure function SphereTriAreaVector(xyzA, xyzB, xyzC)
! Calculates the area of a spherical triangle on the unit sphere
!   NOTE : This function requires function sphereDistance.
	! Calling parameters
	real(KREAL), dimension(3), intent(in) :: xyzA, xyzB, xyzC
	real(KREAL) :: SphereTriAreaVector
	! Local variables
	real(KREAL) :: side1, side2, side3, halfPerimeter, zz

	side1 = SphereArcLength(xyzA,xyzB)
	side2 = SphereArcLength(xyzB,xyzC)
	side3 = SphereArcLength(xyzC,xyzA)

	halfPerimeter = (side1 + side2 + side3)/2.0_KREAL

	zz = tan(halfPerimeter/2.0_KREAL)*tan( (halfPerimeter-side1)/2.0_KREAL )*&
		tan( (halfPerimeter - side2)/2.0_KREAL )*tan( (halfPerimeter - side3)/2.0_KREAL )

	SphereTriAreaVector = 4.0_KREAL * atan2(sqrt(zz),1.0_KREAL)*SphereRadius*SphereRadius
end function

!------------------------------------------------------------------------------
!> @author Peter Bosler
!> @brief Finds the area of a spherical triangle defined by three vertices
!
!> @ingroup SphereGeom
!
!> @param[in] xA double precision real  
!> @param[in] yA double precision real 
!> @param[in] zA
!> @param[in] xB
!> @param[in] yB
!> @param[in] zB
!> @param[in] xC
!> @param[in] yC
!> @param[in] zC
!> @return scalar area
!------------------------------------------------------------------------------
pure function SphereTriAreaComponents(xa,ya,za, xb,yb,zb, xc,yc,zc)
	real(kreal), intent(in) :: xa,ya,za
	real(kreal), intent(in) :: xb,yb,zb
	real(kreal), intent(in) :: xc,yc,zc
	real(kreal) :: SphereTriAreaComponents
	real(kreal) :: s1,s2,s3, halfPerim, zz

	s1 = SphereArcLengthComponents(xa,ya,za,xb,yb,zb)
	s2 = SphereArcLengthComponents(xb,yb,zb,xc,yc,zc)
	s3 = SphereArcLengthComponents(xc,yc,zc,xa,ya,za)

	halfPerim = (s1+s2+s3)/2.0_kreal
	zz = tan( halfPerim/2.0_kreal)*tan( (halfPerim-s1)/2.0_kreal) * &
		 tan( (halfPerim-s2)/2.0_kreal)*tan( (halfPerim-s3)/2.0_kreal)

	SphereTriAreaComponents = 4.0_kreal*atan(sqrt(zz))*SphereRadius*SphereRadius
end function


!------------------------------------------------------------------------------
!> @author Peter Bosler
!> @brief Finds the area of a planar triangle defined by three vertices
!
!> @ingroup SphereGeom
!
!> @param[in] xyzA double precision real size(3); vector on the sphere
!> @param[in] xyzB double precision real size(3); vector on the sphere
!> @param[in] xyzC double precision real size(3); vector on the sphere
!> @return scalar area
!------------------------------------------------------------------------------
pure function PlaneTriArea(xyzA,xyzB,xyzC)
	! Outputs the area of a planar triangle in R3 with vertices xyzA,B,C.
	real(kreal) :: PlaneTriArea
	real(kreal), intent(in) :: xyzA(3), xyzB(3), xyzC(3)
	PlaneTriArea = crossMagnitude( xyzB-xyzA, xyzC-xyzA)/2.0_kreal
end function

!------------------------------------------------------------------------------
!> @author Peter Bosler
!> @brief Finds the longitude of a point on the sphere
!
!> @ingroup SphereGeom
!
!> @param[in] xyz double precision real size(3); vector on the sphere in Cartesian coordinates
!> @return Longitude of xyz
!------------------------------------------------------------------------------
pure function longitudeVector(xyz)
	! Outputs the longitude of a point on the sphere.
	real(KREAL) :: longitudeVector
	real(KREAL), intent(in) :: xyz(3)
	longitudeVector = atan4(xyz(2),xyz(1))
end function

!------------------------------------------------------------------------------
!> @author Peter Bosler
!> @brief Finds the longitude of a point on the sphere
!
!> @ingroup SphereGeom
!
!> @param[in] xyz double precision real size(3); vector on the sphere in Cartesian coordinates
!> @return Longitude of xyz
!------------------------------------------------------------------------------
pure function LongitudeComponents(x,y,z)
	real(kreal), intent(in) :: x, y, z
	real(kreal) :: LongitudeComponents
	LongitudeComponents = atan4(y,x)
end function

!------------------------------------------------------------------------------
!> @author Peter Bosler
!> @brief Finds the longitude of a point on the sphere
!
!> @ingroup SphereGeom
!
!> @param[in] xyz double precision real size(3); vector on the sphere in Cartesian coordinates
!> @return Longitude of xyz
!------------------------------------------------------------------------------
pure function LongitudeComponents2(x,y)
	real(kreal), intent(in) :: x, y
	real(kreal) :: LongitudeComponents2
	LongitudeCOmponents2 = atan4(y,x)
end function

!------------------------------------------------------------------------------
!> @author Peter Bosler
!> @brief Finds the latitude of a point on the sphere
!
!> @ingroup SphereGeom
!
!> @param[in] xyz double precision real size(3); vector on the sphere in Cartesian coordinates
!> @return Latitude of xyz
!------------------------------------------------------------------------------
pure function latitudeVector(xyz)
	! Outputs the latitude of a point on the unit sphere.
	real(KREAL) :: latitudeVector
	real(KREAL), intent(in) :: xyz(3)
	latitudeVector = atan2(xyz(3),sqrt(xyz(1)**2 + xyz(2)**2))
end function

!------------------------------------------------------------------------------
!> @author Peter Bosler
!> @brief Finds the latitude of a point on the sphere
!
!> @ingroup SphereGeom
!
!> @param[in] xyz double precision real size(3); vector on the sphere in Cartesian coordinates
!> @return Latitude of xyz
!------------------------------------------------------------------------------
pure function LatitudeComponents(x,y,z)
	real(kreal), intent(in) :: x, y, z
	real(kreal) :: LatitudeComponents
	LatitudeComponents = atan2(z,sqrt(x*x+y*y))
end function

!------------------------------------------------------------------------------
!> @author Peter Bosler
!> @brief A 4-quadrant inverse tangent function with range 0 to 2*pi instead of -pi to pi.
!
!> @ingroup SphereGeom
!
!> @param[in] y double precision real
!> @param[in] x double precision real
!> @return atan(y/x) in range 0 to 2*pi
!------------------------------------------------------------------------------
pure function atan4(y,x)
	!This function computes the inverse tangent (like atan2) but outputs angles in the range
	! 0 to 2 pi (rather than -pi to pi).
	! Adapted from John Burkhardt: http://people.sc.fsu.edu/~jburkhardt/m_src/halton/atan4.m
	real(kreal), intent(in) :: y,x
	real(kreal) :: atan4
	real(kreal) :: absY, absX, theta
	if ( x == 0.0_kreal) then
		if ( y > 0.0_kreal) then
			atan4 = PI/2.0_kreal
		elseif ( y < 0.0_kreal) then
			atan4 = 3.0_kreal*PI/2.0_kreal
		elseif ( y == 0.0_kreal) then
			atan4 = 0.0_kreal
		endif
	elseif ( y == 0.0_kreal) then
		if ( x > 0.0_kreal) then
			atan4 = 0.0_kreal
		elseif ( x < 0.0_kreal) then
			atan4 = PI
		endif
	else
		absY = abs(y)
		absX = abs(x)
		theta = atan2(absY,absX)
		if ( (x>0.0_kreal) .and. (y>0.0_kreal)) then
			atan4 = theta
		elseif ( (x < 0.0_kreal) .and. (y > 0.0_kreal)) then
			atan4 = pi -theta
		elseif ( (x < 0.0_kreal) .and. (y < 0.0_kreal)) then
			atan4 = pi+ theta
		elseif ( (x > 0.0_kreal) .and. (y<0.0_kreal)) then
			atan4 = 2.0_kreal*PI - theta
		endif
	endif
end function

pure function SphereProjection( xyz )
	real(kreal) :: SphereProjection(3,3)
	real(kreal), intent(in) :: xyz(3)
	SphereProjection(1,1) = 1.0_kreal - xyz(1) * xyz(1)
	SphereProjection(2,1) = - xyz(2) * xyz(1)
	SphereProjection(3,1) = - xyz(3) * xyz(1)
	SphereProjection(1,2) = - xyz(1) * xyz(2)
	SphereProjection(2,2) = 1.0_kreal - xyz(2) * xyz(2)
	SphereProjection(3,2) = - xyz(3) * xyz(2)
	SphereProjection(1,3) = - xyz(1) * xyz(3)
	SphereProjection(2,3) = - xyz(2) * xyz(3)
	SphereProjection(3,3) = 1.0_kreal - xyz(3) * xyz(3)
end function

!------------------------------------------------------------------------------
!> @author Peter Bosler
!> @brief Magnitude of the vector xyzA X xyzB
!
!> @ingroup SphereGeom
!
!> @param[in] xyzA double precision real size(3)
!> @param[in] xyzB double precision real size(3)
!> @return |xyzA x xyzB|
!------------------------------------------------------------------------------
pure function crossMagnitude(xyzA,xyzB)
	! Computes the magnitude of xyzA cross xyzB
	real(kreal) :: crossMagnitude
	real(kreal), intent(in) :: xyzA(3), xyzB(3)
	crossMagnitude = sqrt( (xyzA(2)*xyzB(3) - xyzA(3)*xyzB(2))*(xyzA(2)*xyzB(3) - xyzB(2)*xyzA(3)) + &
	   					   (xyzA(3)*xyzB(1) - xyzA(1)*xyzB(3))*(xyzA(3)*xyzB(1) - xyzA(1)*xyzB(3)) + &
	   					   (xyzA(1)*xyzB(2) - xyzA(2)*xyzB(1))*(xyzA(1)*xyzB(2) - xyzA(2)*xyzB(1)))
end function

!------------------------------------------------------------------------------
!> @author Peter Bosler
!> @brief Computes a vector cross product 
!
!> @ingroup SphereGeom
!
!> @param[in] xyzA double precision real size(3)
!> @param[in] xyzB double precision real size(3)
!> @return xyzA x xyzB double precision real size(3)
!------------------------------------------------------------------------------
pure function crossProduct(xyzA,xyzB)
	! Computes the cross product vector xyzA cross xyzB
	real(kreal) :: crossProduct(3)
	real(kreal), intent(in) :: xyzA(3), xyzB(3)
	crossProduct(1) = xyzA(2)*xyzB(3) - xyzA(3)*xyzB(2)
	crossProduct(2) = xyzA(3)*xyzB(1) - xyzA(1)*xyzB(3)
	crossProduct(3) = xyzA(1)*xyzB(2) - xyzA(2)*xyzB(1)
end function

!------------------------------------------------------------------------------
!> @author Peter Bosler
!> @brief Computes a vector triple product
!
!> @ingroup SphereGeom
!
!> @param[in] xA double precision real size(3)
!> @param[in] xB double precision real size(3)
!> @param[in] xC double precision real size(3)
!> @return \f$ x_A \cdot (x_b \times x_c) $\f double precision real size(3)
!------------------------------------------------------------------------------
pure function Determinant(xA,xB,xC)
	! Computes the vector determinant (triple product) of xA, xB, and xC
	! This result will be positive if xA lies to the left of the directed arc
	! xB to xC, and negative if xA lies to the right of the arc xB->xC.
	real(kreal) :: Determinant
	real(kreal), intent(in) :: xA(3), xB(3), xC(3)
	real(kreal) :: cross(3)
	cross = CrossProduct(xB,xC)
	Determinant = sum(xA*cross)
end function

!------------------------------------------------------------------------------
!> @author Peter Bosler
!> @brief Generates a random point on the surface of a sphere
!
!> @ingroup SphereGeom
!
!> @return random vector satisfying x^2 + y^2 + z^2 = R^2
!------------------------------------------------------------------------------
function RandomSpherePoint()
	! Outputs a random point on the unit sphere.
	! Note : the seed is not initialized, so each sequence of
	!        random points will be the same.
	real(kreal) :: RandomSpherePoint(3)
	real(kreal) :: x, y, z
	!call Random_Seed()
	call Random_Number(x)
	call Random_Number(y)
	call Random_Number(z)
	x = -1.0_kreal + 2.0_kreal*x
	y = -1.0_kreal + 2.0_kreal*y
	z = -1.0_kreal + 2.0_kreal*z
	RandomSpherePoint(1) = x
	RandomSpherePoint(2) = y
	RandomSpherePoint(3) = z
	RandomSpherePoint = RandomSpherePoint/sqrt(sum(RandomSpherePoint*RandomSpherePoint))*SphereRadius
end function

subroutine NormalizePositionVectors(xyz)
	real(kreal), intent(inout) :: xyz(:,:)
	integer(kint) :: j, n
!	real(kreal) :: norm
	if ( size(xyz,1) /= 3) stop 'ERROR : position vector array shape error.'
	n = size(xyz,2)
	do j=1,n
		xyz(:,j) = xyz(:,j)/(sqrt(sum(xyz(:,j)*xyz(:,j))))*SphereRadius
	enddo
end subroutine


!function NorthPoleRotationMatrix(xyz)
!! returns a rotation matrix R so that R*x = (0,0,1)^T
!	real(kreal) :: NorthPoleRotationMatrix(3,3)
!	real(kreal), intent(in) :: xyz(3)
!	real(kreal) :: cy, sy, cx, sx
!
!	cy = sqrt(xyz(2)*xyz(2) + xyz(3)*xyz(3))
!	sy = xyz(1)
!
!	if (abs(cy) > ZERO_TOL) then
!		cx = xyz(3)/cy
!		sx = xyz(2)/cy
!	else
!		cx = 1.0_kreal
!		sx = 0.0_kreal
!	endif
!
!	NorthPoleRotationMatrix(1,1) = cy
!	NorthPoleRotationMatrix(2,1) = 0.0_kreal
!	NorthPoleRotationMatrix(3,1) = sy
!
!	NorthPoleRotationMatrix(1,2) = -sx*sy
!	NorthPoleRotationMatrix(2,2) = cx
!	NorthPoleRotationMatrix(3,2) = cy*sx
!
!	NorthPoleRotationMatrix(1,3) = -cx*sy
!	NorthPoleRotationMatrix(2,3) = -sx
!	NorthPoleRotationMatrix(3,3) = cx*cy
!
!end function

end module
