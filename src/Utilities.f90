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
