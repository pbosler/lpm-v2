module SphereTransportVelocitiesModule

use NumberKindsModule
use SphereGeomModule

implicit none

public 

contains

function MovingVorticesVelocity( x, y, z, t )
	real(kreal), dimension(3) :: MovingVorticesVelocity
	real(kreal), intent(in) :: x, y, z, t
	real(kreal), parameter :: u0 = 2.0_kreal * PI / 12.0_kreal
	real(kreal) :: lat, lon, rho, rhoDenom, omg, u, v
	
	lat = Latitude(x, y, z)
	lon = Longitude(x,y,z)
	
	rho = 3.0_kreal * sqrt(1.0_kreal - cos(lat) * cos(lat) * sin( lon - u0 * t) * sin(lon - u0*t ))
	rhoDenom = rho / ( rho * rho + ZERO_TOL * ZERO_TOL)
	omg = 1.5_kreal * sqrt(3.0_kreal) * u0 * tanh(rho) * rhoDenom / cosh(rho) / cosh(rho)
	
	u = omg * sin( lon - u0 * t) * sin(lat) + u0 * cos(lat)
	v = omg * cos( lon - u0 * t)
	
	MovingVorticesVelocity(1) = -u * sin(lon) - v * sin(lat) * cos(lon)
	MovingVorticesVelocity(2) =  u * cos(lon) - v * sin(lat) * sin(lon)
	MovingVorticesVelocity(3) =  v * cos(lat)
end function


function LauritzenEtalDeformationalVelocity( x, y, z, t )
	real(kreal), dimension(3) :: LauritzenEtalDeformationalVelocity
	real(kreal), intent(in) :: x, y, z, t
	real(kreal), parameter :: RR = 1.0_kreal
	real(kreal), parameter :: TT = 5.0_kreal
	real(kreal) :: u, v, lat, lon
	
	lat = Latitude(x, y, z)
	lon = Longitude(x, y, z)
	u = 10.0_kreal * RR / TT * sin(lon - 2.0_kreal * PI * t / TT) * sin(lon - 2.0_kreal * PI * t / TT ) * &
		sin(2.0_kreal * lat ) * cos( PI * t / TT) + 2.0_kreal * PI * RR / TT * cos(lat)
	v = 10.0_kreal * RR / TT * sin(2.0_kreal * (lon - 2.0_kreal * PI * t / TT)) * cos(lat) * cos(PI * t/TT )
	
	LauritzenEtalDeformationalVelocity(1) = -u * sin(lon) - v * sin(lat) * cos(lon)
	LauritzenEtalDeformationalVelocity(2) =  u * cos(lon) - v * sin(lat) * sin(lon)
	LauritzenEtalDeformationalVelocity(3) =  v * cos(lat)
end function

function LauritzenEtalDivergentFlowVelocity( x, y, z, t )
	real(kreal), dimension(3) :: LauritzenEtalDivergentFlowVelocity
	real(kreal), intent(in) :: x, y, z, t
	real(kreal), parameter :: RR = 1.0_kreal
	real(kreal), parameter :: TT = 5.0_kreal
	real(kreal), parameter :: relPeriod = RR / TT
	real(kreal) :: u, v, lat, lon, lamPrime
	
	lat = Latitude(x, y, z)
	lon = Longitude(x, y, z)
	lamPrime = lon - 2.0_kreal * PI * t / TT
	
	u = - 5.0_kreal * relPeriod * sin( 0.5_kreal * lamPrime) * sin( 0.5_kreal * lamPrime ) * sin( 2.0_kreal * lat ) * &
		cos(lat) * cos(lat) * cos( PI * t / TT ) + 2.0_kreal * PI * relPeriod * cos(lat)
	v = 2.5_kreal * relPeriod * sin(lamPrime) * cos(lat) * cos(lat) * cos(lat) * cos( PI * t / TT )
	
	LauritzenEtalDivergentFlowVelocity(1) = -u * sin(lon) - v * sin(lat) * cos(lon)
	LauritzenEtalDivergentFlowVelocity(2) =  u * cos(lon) - v * sin(lat) * sin(lon)
	LauritzenEtalDivergentFlowVelocity(3) =  v * cos(lat)
end function

function LauritzenEtalDivergentFlowDivergence( x, y, z, t )
	real(kreal) :: LauritzenEtalDivergentFlowDivergence
	real(kreal), intent(in) :: x, y, z, t
	real(kreal), parameter :: TT = 5.0_kreal
	real(kreal) :: u, v, lat, lon, lamPrime
	
	lat = Latitude(x, y, z)
	lon = Longitude(x, y, z)
	lamPrime = lon - 2.0_kreal * PI * t / TT
	
	LauritzenEtalDivergentFlowDivergence = - 15.0_kreal * cos(lat) * cos(lat) * cos( PI * t / TT ) * sin(lat) * &
		sin(lamPrime) / TT
end function 

function RH4Velocity( x, y, z, t )
	real(kreal), dimension(3) :: RH4Velocity
	real(kreal), intent(in) :: x, y, z, t
	real(kreal) :: u, v, lat, lon
	real(kreal), parameter :: u0 = 0.0_kreal
	real(kreal), parameter :: amp = 1.0_kreal
	
	lat = Latitude(x, y, z)
	lon = Longitude(x, y, z)
	
	u = cos(lat) * ( u0 - 0.5_kreal * amp * cos( 4.0_kreal * lon ) * cos(lat) * cos(lat) * &
		(-3.0_kreal + 5.0_kreal * cos(2.0_kreal * lat )) )
	v = -4.0_kreal * amp * cos(lat) * cos(lat) * cos(lat) * sin(4.0_kreal * lon ) * sin(lat)
	
	RH4Velocity(1) = -u * sin(lon) - v * sin(lat) * cos(lon)
	RH4Velocity(2) =  u * cos(lon) - v * sin(lat) * sin(lon)
	RH4Velocity(3) =  v * cos(lat)
end function

end module