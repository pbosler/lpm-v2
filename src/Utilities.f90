module UtilitiesModule

use NumberKindsModule

implicit none

public

real(kreal), parameter :: dzero = 0.0_kreal !< timesaver
real(kreal), parameter :: one = 1.0_kreal !< timesaver

logical(klog), save :: testPass = .TRUE. !< global test variable


!> @brief Provides a generic function interface
interface
    real(8) function scalarFnOf2DSpace( x, y)
        real(8), intent(in) :: x, y
    end function
end interface

!> @brief Provides a generic function interface
interface
    real(8) function scalarFnOf2DSpaceAndTime( x, y, t)
        real(8), intent(in) :: x, y, t
    end function
end interface

!> @brief Provides a generic function interface
interface
    real(8) function scalarFnOf3DSpace( x, y, z)
        real(8), intent(in) :: x, y, z
    end function
end interface

!> @brief Provides a generic function interface
interface
    function scalarFnOf3DSpaceAndTime( x, y, z, t )
        real(8) :: scalarFnOf3DSpaceAndTime
        real(8), intent(in) :: x, y, z, t
    end function
end interface

!> @brief Provides a generic function interface
interface
    function vectorFnOf2DSpace( x, y )
        real(8), dimension(2) :: vectorFnOf2DSpace
        real(8), intent(in) :: x, y
    end function
end interface

!> @brief Provides a generic function interface
interface
    function vectorFnOf2DSpaceAndTime( x, y, t )
        real(8), dimension(2) :: vectorFnOf2DSpaceAndTime
        real(8), intent(in) :: x, y, t
    end function
end interface

!> @brief Provides a generic function interface
interface
    function vectorFnOf3DSpace( x, y, z)
        real(8), dimension(3) :: vectorFnOf3DSpace
        real(8), intent(in) :: x, y, z
    end function
end interface

!> @brief Provides a generic function interface
interface
    function vectorFnOf3DSpaceAndTime( x, y, z, t )
        real(8), dimension(3) :: vectorFnOf3DSpaceAndTime
        real(8), intent(in) :: x, y, z, t
    end function
end interface

contains

function geomString(geomKind)
    character(len=56) :: geomString
    integer(kint), intent(in) :: geomKind
    if (geomKind == PLANAR_GEOM) then
        geomString = "planar_geometry"
    elseif(geomKind == SPHERE_GEOM) then
        geomString = "sphere_geometry"
    else
        geomString = "invalid_geometry_kind"
    endif
end function

function faceString(faceKind)
    character(len=56) :: faceString
    integer(kint), intent(in) :: faceKind
    if (faceKind == QUAD_PANEL) then
        faceString = "quad_panel"
    elseif(faceKind == TRI_PANEL) then
        faceString = "tri_panel"
    elseif(faceKind == QUAD_CUBIC_PANEL) then
        faceString = "quad_cubic_panel"
    else
        faceString = "invalid faceKind"
    endif
end function

pure function square(x)
    real(kreal) :: square
    real(kreal), intent(in) :: x
    square = x*x
end function

pure function cube(x)
    real(kreal) :: cube
    real(kreal), intent(in) :: x
    cube = x*x*x
end function

end module
