module SphereTracersModule

use NumberKindsModule
use SphereGeomModule

implicit none

public 

contains

function SlottedCylinderTracer( x, y, z )
	real(kreal) :: SlottedCylinderTracer
	real(kreal), intent(in) :: x, y, z 
	real(kreal), parameter :: xx1 = -0.866025403784439_kreal, &	
							  yy1 = 0.5_kreal, zz1 = 0.0_kreal
	real(kreal), parameter :: xx2 = -0.866025403784439_kreal, &
							  yy2 = -0.5_kreal, zz2 = 0.0_kreal
	real(kreal), parameter :: lat1 = 0.0_kreal, long1 = 5.0_kreal*PI/6.0_kreal
	real(kreal), parameter :: lat2 = 0.0_kreal, long2 = 7.0_kreal*PI/6.0_kreal
	real(kreal), parameter :: RR = 0.5_kreal, b = 0.1_kreal, c = 1.0_kreal
	real(kreal) :: r1, r2, lat, lon
	
	lat = Latitude(x, y, z)
	lon = Longitude(x, y, z)
	r1 = SphereArcLength([x,y,z], [xx1,yy1,zz1])
	r2 = SphereArcLength([x,y,z], [xx2,yy2,zz2])
	
	SlottedCylinderTracer = b
	
	if ( r1 <= RR ) then
		if ( abs(lon - long1) >= RR / 6.0_kreal ) then
			SlottedCylinderTracer = c
		elseif ( lat - lat1 < -5.0_kreal * RR / 12.0_kreal ) then
			SlottedCylinderTracer = c
		endif
	endif
	
	if ( r2 <= RR ) then
		if ( abs(lon - long2) >= RR / 6.0_kreal ) then
			SlottedCylinderTracer = c
		elseif ( lat - lat2 > 5.0_kreal * RR / 12.0_kreal ) then
			SlottedCylinderTracer = c
		endif
	endif
end function

function CosineBellsTracer( x, y, z )
	real(kreal) CosineBellsTracer
	real(kreal), intent(in) :: x, y, z
	real(kreal), parameter :: xx1 = -0.866025403784439_kreal, &
							  yy1 = 0.5_kreal, zz1 = 0.0_kreal
	real(kreal), parameter :: xx2 = -0.866025403784439_kreal,&
		  					  yy2 = -0.5_kreal, zz2 = 0.0_kreal
	real(kreal), parameter :: hmax = 1.0_kreal, RR = 0.5_kreal, b = 0.1_kreal, c = 0.9_kreal
	real(kreal) :: r1, r2, h1, h2
	r1 = SphereArcLength([x,y,z], [xx1,yy1,zz1])
	r2 = SphereArcLength([x,y,z], [xx2,yy2,zz2])
	h1 = 0.5_kreal * hmax * (1.0_kreal + cos(PI*r1/RR))
	h2 = 0.5_kreal * hmax * (1.0_kreal + cos(PI*r2/RR))
	if ( r1 < RR ) then
		CosineBellsTracer = b + c * h1
	elseif ( r2 < RR) then	
		CosineBellsTracer = b + c * h2
	else 
		CosineBellsTracer = b
	endif
end function

function MovingVorticesTracer(x, y, z, t) 
	real(kreal) :: MovingVorticesTracer
	real(kreal), intent(in) :: x, y, z, t
	real(kreal) :: lat, lon, rho, omg, lamPrime, rhoDenom
	real(kreal), parameter :: u0 = 2.0_kreal * PI / 12.0_kreal
	real(kreal), parameter :: Omega = u0
	
	lon = Longitude(x,y,z)
	lat = Latitude(x,y,z)
	
	lamPrime = atan4( - cos(lon - Omega * t), tan(lat) )
	
	rho = 3.0_kreal * sqrt( 1.0_kreal - cos(lat) * cos(lat) * sin( lon - Omega * t) * sin( lon - Omega * t))
	rhoDenom = rho / ( rho * rho + ZERO_TOL * ZERO_TOL )
	
	omg = 1.5_kreal * sqrt(3.0_kreal) * u0 * tanh(rho) / cosh(rho) / cosh(rho)
	
	MovingVorticesTracer = 1.0_kreal - tanh( 0.2_kreal * rho * sin( lamPrime - omg * t ))	
end function 

function GaussianHillsTracer( x, y, z )
	real(kreal) :: GaussianHillsTracer
	real(kreal), intent(in) :: x, y, z
	real(kreal), parameter :: beta = 5.0_kreal
	real(kreal), parameter :: center1x = cos(5.0_kreal * PI / 6.0_kreal), &
							  center1y = sin(5.0_kreal * PI / 6.0_kreal ), &
							  center1z = 0.0_kreal
	real(kreal), parameter :: center2x = cos(7.0_kreal * PI / 6.0_kreal), &
						   	  center2y = sin(7.0_kreal * PI / 6.0_kreal ), &
						   	  center2z = 0.0_kreal
	real(kreal), parameter :: hmax = 0.95_kreal
	real(kreal) :: h1, h2
	
	h1 = hmax * exp( -beta * ( (x - center1x) * (x - center1x) + (y - center1y) * (y - center1y) + &
							   (z - center1z) * (z - center1z) ) )
	h2 = hmax * exp( -beta * ( (x - center2x) * (x - center2x) + (y - center2y) * (y - center2y) + &
							   (z - center2z) * (z - center2z) ) )
	GaussianHillsTracer = h1 + h2
end function

function BlockMVector(x, y, z)
	real(kreal) :: BlockMVector
	real(kreal), intent(in) :: x, y, z
	real(kreal), parameter :: unit = PI/6.0_kreal
	real(kreal) :: lat, lon, xx, yy

	lat = Latitude(x,y,z)
	lon = atan2(y,x)
	xx = lon/unit
	yy = lat/unit
	BlockMVector = ScaledBlockM(xx, yy)
end function

function ScaledBlockM(x,y)
	real(kreal) :: ScaledBlockM
	real(kreal), intent(in) :: x, y
	real(kreal), parameter :: m = 1.6625_kreal, &
							  vert1 = 0.38421_kreal, &
							  horiz1 = 0.315789_kreal, &
							  horiz2 = 0.8421_kreal,&
							  horiz3 = 0.4736842,&
							  top = 0.7_kreal, &
							  right = 1.0_kreal, &
							  bottom = -0.7_kreal, &
							  left = -1.0_kreal, &
							  one = 1.0_kreal, &
							  zero = 0.0_kreal

	ScaledBlockM = 0.1_kreal
	! Left Foot
	if ( (( y < -vert1) .AND. ( y > bottom )) .AND. ( (x < - horiz1) .AND. (x > left))) ScaledBlockM = one
	! Right Foot
	if ( ( (y < -vert1) .AND. ( y > bottom )) .AND. ( (x > horiz1) .AND. (x < right))) ScaledBlockM = one
	! Left Leg
	if ( ( (y > bottom) .AND. (y < top) ) .AND. ( (x > -horiz2) .AND. (x < -horiz3 ))) ScaledBlockM = one
	! Right Leg
	If ( ( (y > bottom) .AND. (y<top)) .AND. ( (x>horiz3) .AND. (x<horiz2))) ScaledBlockM = one
	! Left top
	If ( ( (y > vert1) .AND. (y<top)) .AND. ( (x<-horiz1) .AND. (x> left))) ScaledBlockM = one
	! Right top
	If ( ( (y > vert1) .AND. (y<top)) .AND. ( (x<right) .AND. (x>horiz1))) ScaledBlockM = one
	! Left triangle
	If ( (x <= zero) .AND. ( y < top) ) then
		if (( y > -m*x - top) .AND. ( y < -m*(x + horiz1) + top)) ScaledBlockM = one
	endif
	If ( (x >= zero) .AND. ( y<top)) then
		if ((y > m*x - top) .AND. (y<m*(x-horiz1)+top)) ScaledBlockM = one
	endif
end function

pure function InitLatTracer( x0, y0, z0 )
	real(kreal) :: InitLatTracer
	real(kreal), intent(in) :: x0, y0, z0
	InitLatTracer = Latitude( x0, y0, z0 )
end function

pure function RH54Vorticity(x0, y0, z0 )
	real(kreal) :: RH54Vorticity
	real(kreal), intent(in) :: x0, y0, z0
	!
	real(kreal) :: lat, lon
	
	lat = Latitude(x0, y0, z0)
	lon = Longitude(x0, y0, z0)
	
	RH54Vorticity = 30.0_kreal * cos( 4.0_kreal * lon ) * cos( lat ) **4 * sin(lat)
end function


end module
