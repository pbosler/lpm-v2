module UtilitiesModule

use NumberKindsModule

implicit none

public

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
