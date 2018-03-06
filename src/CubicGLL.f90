module CubicGLLModule

use NumberKindsModule
use UtilitiesModule

implicit none 

public

real(kreal), parameter :: r5 = sqrt(5.0_kreal)
real(kreal), parameter :: oor5 = 1.0_kreal / r5
real(kreal), dimension(4), parameter :: gll_cubic_qp = (/ -1.0_kreal, -oor5, oor5, 1.0_kreal /)
real(kreal), dimension(4), parameter :: gll_cubic_qw = (/ 1.0_kreal, 5.0_kreal, 5.0_kreal, 1.0_kreal /) / 6.0_kreal
real(kreal), dimension(2,12), parameter :: quad16_vertex_qp = reshape([ gll_cubic_qp(1), gll_cubic_qp(4), &
                                                                gll_cubic_qp(1), gll_cubic_qp(3), &
                                                                gll_cubic_qp(1), gll_cubic_qp(2), &
                                                                gll_cubic_qp(1), gll_cubic_qp(1), &
                                                                gll_cubic_qp(2), gll_cubic_qp(1), &
                                                                gll_cubic_qp(3), gll_cubic_qp(1), &
                                                                gll_cubic_qp(4), gll_cubic_qp(1), &
                                                                gll_cubic_qp(4), gll_cubic_qp(2), &
                                                                gll_cubic_qp(4), gll_cubic_qp(3), &
                                                                gll_cubic_qp(4), gll_cubic_qp(4), &
                                                                gll_cubic_qp(3), gll_cubic_qp(4), &
                                                                gll_cubic_qp(2), gll_cubic_qp(4)], shape(quad16_vertex_qp))
real(kreal), dimension(2,4), parameter :: quad16_center_qp = reshape([ gll_cubic_qp(2), gll_cubic_qp(3), &
                                                                gll_cubic_qp(2), gll_cubic_qp(2), &
                                                                gll_cubic_qp(3), gll_cubic_qp(2), &
                                                                gll_cubic_qp(3), gll_cubic_qp(3)], shape(quad16_center_qp))

real(kreal), dimension(12), parameter :: quad16_vertex_qw =(/gll_cubic_qw(1)*gll_cubic_qw(4), &
                                              gll_cubic_qw(1)*gll_cubic_qw(3), &
                                              gll_cubic_qw(1)*gll_cubic_qw(2), &
                                              gll_cubic_qw(1)*gll_cubic_qw(1), &
                                              gll_cubic_qw(2)*gll_cubic_qw(1), &
                                              gll_cubic_qw(3)*gll_cubic_qw(1), &
                                              gll_cubic_qw(4)*gll_cubic_qw(1), &
                                              gll_cubic_qw(4)*gll_cubic_qw(2), &
                                              gll_cubic_qw(4)*gll_cubic_qw(3), &
                                              gll_cubic_qw(4)*gll_cubic_qw(4), &
                                              gll_cubic_qw(3)*gll_cubic_qw(4), &
                                              gll_cubic_qw(2)*gll_cubic_qw(4) /)
real(kreal), dimension(4), parameter :: quad16_center_qw =(/gll_cubic_qw(2)*gll_cubic_qw(3), &
                                              gll_cubic_qw(2)*gll_cubic_qw(2), &
                                              gll_cubic_qw(3)*gll_cubic_qw(2), &
                                              gll_cubic_qw(3)*gll_cubic_qw(3) /)

contains
                                            
pure function bilinearMap(vertXyz, s1, s2)
    real(kreal), dimension(3) :: bilinearMap
    real(kreal), dimension(3,4), intent(in) :: vertXyz
    real(kreal), intent(in) :: s1, s2

    bilinearMap = 0.25_kreal * ( (1.0_kreal-s1)*(1.0_kreal+s2)*vertXyz(:,1) + (1.0_kreal-s1)*(1.0_kreal-s2)*vertXyz(:,2) + &
        (1.0_kreal+s1)*(1.0_kreal-s2)*vertXyz(:,3) + (1.0_kreal+s1)*(1.0_kreal+s2)*vertXyz(:,4))
end function

function sphereQuadMap(vertXyz, s1, s2)
    real(kreal), dimension(3) :: sphereQuadMap
    real(kreal), dimension(3,4), intent(in) :: vertXyz
    real(kreal), intent(in) :: s1, s2
    !
    real(kreal) :: norm
    sphereQuadMap = bilinearMap(vertXyz, s1, s2)
    norm = sqrt(sum(vertXyz*vertXyz))
    sphereQuadMap = sphereQuadMap / norm
end function

pure function quad16InteriorPts(vertXYZ)
    real(kreal), dimension(3,4) :: quad16InteriorPts
    real(kreal), dimension(3,4), intent(in) :: vertXyz
    quad16InteriorPts(:,1) = bilinearMap(vertXyz, gll_cubic_qp(2), gll_cubic_qp(3))
    quad16InteriorPts(:,2) = bilinearMap(vertXyz, gll_cubic_qp(2), gll_cubic_qp(2))
    quad16InteriorPts(:,3) = bilinearMap(vertXyz, gll_cubic_qp(3), gll_cubic_qp(2))
    quad16InteriorPts(:,4) = bilinearMap(vertXyz, gll_cubic_qp(3), gll_cubic_qp(3))
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
         (1.0_kreal+s1)*vertXyz(2,3) + (1.0_kreal+s1)*vertXyz(2,4))
    bilinearPlaneJacobian = abs(a*d - b*c)
end function

elemental function gll_cubic_basis1(x)
    real(kreal) :: gll_cubic_basis1
    real(kreal), intent(in) :: x
    gll_cubic_basis1 = (-1.0_kreal + x + 5.0_kreal * x**2 - 5.0_kreal * x**3) / 8.0_kreal
end function

elemental function gll_cubic_basis2(x)
    real(kreal) :: gll_cubic_basis2
    real(kreal), intent(in) :: x
    gll_cubic_basis2 = 5.0_kreal * (x-1.0_kreal)*(x+1.0_kreal)* (r5 *x-1.0_kreal) / 8.0_kreal
end function

elemental function gll_cubic_basis3(x)
    real(kreal) :: gll_cubic_basis3
    real(kreal), intent(in) :: x 
    gll_cubic_basis3 = -r5 * (x-1.0_kreal)*(x+1.0_kreal)*(r5 + 5.0_kreal * x) / 8.0_kreal
end function

elemental function gll_cubic_basis4(x)
    real(kreal) :: gll_cubic_basis4
    real(kreal), intent(in) :: x
    gll_cubic_basis4 = (-1.0_kreal - x + 5.0_kreal *x**2 + 5.0_kreal*x**3)/8.0_kreal
end function

end module