module UtilitiesModule

use NumberKindsModule

implicit none

public

real(kreal), parameter :: oor5 = 1.0_kreal / sqrt(5.0_kreal)

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

pure function bilinearMap(vertXyz, s1, s2)
    real(kreal), dimension(3) :: bilinearMap
    real(kreal), dimension(3,4), intent(in) :: vertXyz
    real(kreal), intent(in) :: s1, s2

    bilinearMap = 0.25_kreal * ( (1.0_kreal-s1)*(1.0_kreal+s2)*vertXyz(:,1) + (1.0_kreal-s1)*(1.0_kreal-s2)*vertXyz(:,2) + &
        (1.0_kreal+s1)*(1.0_kreal-s2)*vertXyz(:,3) + (1.0_kreal+s1)*(1.0_kreal+s2)*vertXyz(:,4))
end function

pure function bilinearPlaneJacobian(vertXyz, s1, s2)
    real(kreal) :: bilinearPlaneJacobian
    real(kreal), dimension(3,4), intent(in) :: vertXyz
    real(kreal), intent(in) :: s1, s2
    !
    real(kreal) :: a, b, c, d

    a = 0.25_kreal * (-(1.0_kreal+s2)*vertXyz(1,1) -(1.0_kreal-s2)*vertXyz(1,2) +&
         (1.0_kreal-s2)*vertXyz(1,3) + (1.0_kreal+s2)*vertXyz(1,4))

    b = 0.25_kreal * (-(1.0_kreal+s2)*vertXyz(2,1) -(1.0_kreal-s2)*vertXyz(2,2) +&
         (1.0_kreal-s2)*vertXyz(2,3) + (1.0_kreal+s2)*vertXyz(2,4))

    c = 0.25_kreal * ((1.0_kreal-s1)*vertXyz(1,1) - (1.0_kreal-s1)*vertXyz(1,2) - &
         (1.0_kreal+s1)*vertXyz(1,3) + (1.0_kreal+s2)*vertXyz(1,4))

    d = 0.25_kreal * ((1.0_kreal-s1)*vertXyz(2,1) - (1.0_kreal-s1)*vertXyz(2,2) - &
         (1.0_kreal+s1)*vertXyz(2,3) + (1.0_kreal+s2)*vertXyz(2,4))
    bilinearPlaneJacobian = abs(a*d - b*c)
end function

end module
