module CubicGLLModule

use NumberKindsModule
use UtilitiesModule

implicit none

public

integer(kint), parameter :: MAX_ITER = 20 !< maximum iterations for Newton's method (used in inverseSphereMap)

real(kreal), parameter :: r5 = sqrt(5.0_kreal) !< timesaver
real(kreal), parameter :: oor5 = 1.0_kreal / r5 !< timesaver

real(kreal), dimension(4), parameter :: gll_cubic_qp = (/ -1.0_kreal, -oor5, oor5, 1.0_kreal /) !< quadrature pts for cubic GLL
real(kreal), dimension(4), parameter :: gll_cubic_qw = (/ 1.0_kreal, 5.0_kreal, 5.0_kreal, 1.0_kreal /) / 6.0_kreal !< quadrature weights for cubic GLL

!> quadrature points re-indexed to match layout of quad16 faces, edges, and particles
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

!> quadrature points re-indexed to match layout of quad16 faces, edges, and particles
real(kreal), dimension(2,4), parameter :: quad16_center_qp = reshape([ gll_cubic_qp(2), gll_cubic_qp(3), &
                                                                gll_cubic_qp(2), gll_cubic_qp(2), &
                                                                gll_cubic_qp(3), gll_cubic_qp(2), &
                                                                gll_cubic_qp(3), gll_cubic_qp(3)], shape(quad16_center_qp))
!> quadrature weights re-indexed to match layout of quad16 faces, edges, and particles
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
!> quadrature weights re-indexed to match layout of quad16 faces, edges, and particles
real(kreal), dimension(4), parameter :: quad16_center_qw =(/gll_cubic_qw(2)*gll_cubic_qw(3), &
                                              gll_cubic_qw(2)*gll_cubic_qw(2), &
                                              gll_cubic_qw(3)*gll_cubic_qw(2), &
                                              gll_cubic_qw(3)*gll_cubic_qw(3) /)

contains

!> Computes the matrix associated with bilinearMap
pure function bilinearTransformMatrix(vertXyz)
    real(kreal), dimension(3,4) :: bilinearTransformMatrix
    real(kreal), dimension(3,4), intent(in) :: vertXyz
    real(kreal), dimension(4,4), parameter :: P = reshape([one, -one, one, -one, & !col 1
                                                  -one, one, dzero, dzero, & ! col2
                                                  -one, dzero, dzero, one, & ! col3
                                                  one, dzero, dzero, dzero], shape(P)) ! col4
    bilinearTransformMatrix = matmul(vertXyz, P)
end function

!> @brief Finds the point X inside a bilinear quadrilateral,
!> Given the vertices of a quadrilateral in R^d (d=2,3), and reference coordinates (r1,r2) in [-1,1]x[-1,1],
!>
!> @param vertXyz Vertex coordinates of quadrilateral, in CCW order
!> @param s1, s2 reference element coordinates in [-1,1]x[-1,1]
pure function bilinearMap(vertXyz, s1, s2)
    real(kreal), dimension(3) :: bilinearMap
    real(kreal), dimension(3,4), intent(in) :: vertXyz
    real(kreal), intent(in) :: s1, s2
    !
    real(kreal) :: a, b
    real(kreal), dimension(3,4) :: tMatrix
    integer(kint) :: i

    bilinearMap = 0.25_kreal * ( (1.0_kreal-s1)*(1.0_kreal+s2)*vertXyz(:,1) + &
                                 (1.0_kreal-s1)*(1.0_kreal-s2)*vertXyz(:,2) + &
                                (1.0_kreal+s1)*(1.0_kreal-s2)*vertXyz(:,3) + &
                                (1.0_kreal+s1)*(1.0_kreal+s2)*vertXyz(:,4))

end function

!> Inverts bilinearMap to find reference coordinates (r1,r2) in [-1,1]^2
!> Folows C. Hua, 1990, Finite Elem. Anal. Design 7:159--166.
function invertBilinear(vertXyz, xyz)
    real(kreal), dimension(2) :: invertBilinear
    real(kreal), dimension(3,4), intent(in) :: vertXyz
    real(kreal), dimension(3), intent(in) :: xyz
    !
    real(kreal), dimension(3) :: aa, bb, cc, dd
    integer(kint), dimension(2) :: ij
    real(kreal) :: xi, eta, a_b, a_c, a_d, b_c, b_d, c_d
    integer(kint) :: i, j

    aa = -vertXyz(:,1) + vertXyz(:,2) - vertXyz(:,3) + vertXyz(:,4)
    bb = -vertXyz(:,1) - vertXyz(:,2) + vertXyz(:,3) + vertXyz(:,4)
    cc = vertXyz(:,1) - vertXyz(:,2) - vertXyz(:,3) + vertXyz(:,4)

    dd = 4.0_kreal * xyz - sum(vertXyz,2)

    ij = pickBilinearIJ(vertXyz)
    i = ij(1)
    j = ij(2)

    xi = 999.0d20
    eta = 999.0d20

    a_b = twoByTwoDeterminant([aa(i),aa(j)],[bb(i),bb(j)])
    a_c = twoByTwoDeterminant([aa(i),aa(j)],[cc(i),cc(j)])
    a_d = twoByTwoDeterminant([aa(i),aa(j)],[dd(i),dd(j)])
    b_c = twoByTwoDeterminant([bb(i),bb(j)],[cc(i),cc(j)])
    b_d = twoByTwoDeterminant([bb(i),bb(j)],[dd(i),dd(j)])
    c_d = twoByTwoDeterminant([cc(i),cc(j)],[dd(i),dd(j)])

    if ( abs(aa(i)) < ZERO_TOL) then !I
        if (abs(aa(j)) < ZERO_TOL) then ! rectangle
            !I.A
            xi = -c_d / (aa(i)*dd(j) + b_c)
            eta = b_d /(aa(j)*dd(i) + b_c)
        elseif (abs(cc(i)) < ZERO_TOL) then
            !I.B(a)
            xi = dd(i) / bb(i)
            eta = (bb(i)*dd(j)-bb(j)*dd(i))/(aa(j)*dd(i)+bb(i)*cc(j))
        else
            !I.B(b)
            xi = quadraticRoot(aa(j)*bb(i), cc(j)*bb(i)-aa(j)*dd(i)-bb(j)*cc(i), dd(j)*cc(i)-cc(j)*dd(i))
            eta = (dd(i)-bb(i)*xi)/cc(i)
        endif
    else !II
        if (abs(aa(j)) > ZERO_TOL) then !II.A
            if (abs(a_b) > ZERO_TOL) then ! II.A(a)
                if (abs(a_c) > ZERO_TOL) then
                    !II.A.(a).1
                    xi = quadraticRoot(a_b, (-b_c-a_d), -c_d)
                    eta = (a_d-a_b*xi)/(a_c)
                else
                    !II.A.(a).2
                    xi = a_d/a_b
                    eta = -aa(i)*b_d/(cc(i)*a_b + aa(i)*a_d)
                endif
            else ! II.A(b)
                xi = -aa(i)*c_d/(bb(i)*a_c + aa(i)*a_d)
                eta = a_d/a_c
            endif
        else !II.B
            if (abs(bb(j)) < ZERO_TOL) then
                !II.B(a)
                xi = -c_d / (aa(i)*dd(j)+bb(i)*cc(j))
                eta = dd(j) / cc(j)
            else
                !II.B(b)
                xi = quadraticRoot(aa(i)*bb(j), (cc(i)*bb(j) - aa(i)*dd(j) - bb(i)*cc(j)), (dd(i)*cc(j)-cc(i)*dd(j)))
                eta = (dd(j)-bb(j)*xi)/cc(j)
            endif
        endif
    endif
    invertBilinear = [xi, eta]
end function

pure function twoByTwoDeterminant(c1, c2)
    real(kreal) :: twoByTwoDeterminant
    real(kreal), dimension(2), intent(in) :: c1, c2
    twoByTwoDeterminant = c1(1)*c2(2) - c1(2)*c2(1)
end function

function quadraticRoot(a, b, c)
    real(kreal) :: quadraticRoot
    real(kreal), intent(in) :: a, b, c
    !
    real(kreal) :: r1, r2, disc, aa, rdisc

    quadraticRoot = -999.d20

    aa = 2.0_kreal * a
    disc = b*b - 4.0_kreal * a * c
    if (abs(disc) < ZERO_TOL) then
        disc = dzero
    elseif (disc < -ZERO_TOL) then
        print *, "quadratic root error: complex root."
    endif
    rdisc = sqrt(disc)
    r1 = (-b + rdisc)/aa
    r2 = (-b - rdisc)/aa
    if (-one <= r1 .and. r1 <= one) quadraticRoot = r1
    if (-one <= r2 .and. r2 <= one) quadraticRoot = r2
end function

!> @brief Find the indices of the dimensions with the largest extent
!> for a bilinear quadrilateral in R^3; helps with stability in invertBilinear.
pure function pickBilinearIJ(vertXyz)
    integer(kint), dimension(2) :: pickBilinearIJ
    real(kreal), dimension(3,4), intent(in) :: vertXyz
    !
    real(kreal) :: dx, dy, dz

    dx = maxval(vertXyz(1,:)) - minval(vertXyz(1,:))
    dy = maxval(vertXyz(2,:)) - minval(vertXyz(2,:))
    dz = maxval(vertXyz(3,:)) - minval(vertXyz(3,:))

    ! start with x, y
    pickBilinearIJ(1) = 1
    pickBilinearIJ(2) = 2
    if (dy < dx .and. dy < dz) then ! skip y
        pickBilinearIJ(2) = 3
    elseif( dx < dy .and. dx < dz) then !skip x
        pickBilinearIJ(1) = 2
        pickBilinearIJ(2) = 3
    endif
end function

!> @brief Maps the reference element RQ = [-1,1]^2 to a spherical quadrilateral
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

!> @brief Computes the Jacobian matrix associated with sphereQuadMap
function sphereMapJacobian(xyz, vertXyz, s1, s2)
    real(kreal), dimension(3,3) :: sphereMapJacobian
    real(kreal), dimension(3), intent(in) :: xyz
    real(kreal), dimension(3,4), intent(in) :: vertXyz
    real(kreal), intent(in) :: s1, s2
    sphereMapJacobian = dzero
end function

!> @brief Computes the residual of sphereQuadMap; used for Newton iterations in invertSphereQuadMap.
pure function sphereMapResidual(xyz, vertXyz, s1, s2)
    real(kreal), dimension(3) :: sphereMapResidual
    real(kreal), dimension(3), intent(in) :: xyz
    real(kreal), dimension(3,4), intent(in) :: vertXyz
    real(kreal), intent(in) :: s1, s2
    !
    integer(kint) :: i
    real(kreal) :: norm
    real(kreal), dimension(3) :: blm

    blm = bilinearMap(vertXyz, s1, s2)
    norm = sqrt(sum(blm*blm))
    do i=1,3
        sphereMapResidual(i) = blm(i)/norm - xyz(i)
    enddo
end function

function invertSphereMap(xyz, vertXyz)
    real(kreal), dimension(2) :: invertSphereMap
    real(kreal), dimension(3) :: xyz
    real(kreal), dimension(3,4) :: vertXyz
    invertSphereMap = dzero
end function



!> @brief  Locates the interior points of a quadrilateral panel in R^2
pure function quad16InteriorPts(vertXYZ)
    real(kreal), dimension(3,4) :: quad16InteriorPts
    real(kreal), dimension(3,4), intent(in) :: vertXyz
    quad16InteriorPts(:,1) = bilinearMap(vertXyz, gll_cubic_qp(2), gll_cubic_qp(3))
    quad16InteriorPts(:,2) = bilinearMap(vertXyz, gll_cubic_qp(2), gll_cubic_qp(2))
    quad16InteriorPts(:,3) = bilinearMap(vertXyz, gll_cubic_qp(3), gll_cubic_qp(2))
    quad16InteriorPts(:,4) = bilinearMap(vertXyz, gll_cubic_qp(3), gll_cubic_qp(3))
end function

!> @brief Computes the Jacobian of the bilinearMap for planar panels.
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

!> @brief Lagrange interpolating polynomial, GLL basis function, cardinal function on 1d reference element [-1,1]
elemental function gll_cubic_basis1(x)
    real(kreal) :: gll_cubic_basis1
    real(kreal), intent(in) :: x
    gll_cubic_basis1 = (-1.0_kreal + x + 5.0_kreal * x**2 - 5.0_kreal * x**3) / 8.0_kreal
end function

!> @brief Lagrange interpolating polynomial, GLL basis function, cardinal function on 1d reference element [-1,1]
elemental function gll_cubic_basis2(x)
    real(kreal) :: gll_cubic_basis2
    real(kreal), intent(in) :: x
    gll_cubic_basis2 = 5.0_kreal * (x-1.0_kreal)*(x+1.0_kreal)* (r5 *x-1.0_kreal) / 8.0_kreal
end function

!> @brief Lagrange interpolating polynomial, GLL basis function, cardinal function on 1d reference element [-1,1]
elemental function gll_cubic_basis3(x)
    real(kreal) :: gll_cubic_basis3
    real(kreal), intent(in) :: x
    gll_cubic_basis3 = -r5 * (x-1.0_kreal)*(x+1.0_kreal)*(r5 + 5.0_kreal * x) / 8.0_kreal
end function

!> @brief Lagrange interpolating polynomial, GLL basis function, cardinal function on 1d reference element [-1,1]
elemental function gll_cubic_basis4(x)
    real(kreal) :: gll_cubic_basis4
    real(kreal), intent(in) :: x
    gll_cubic_basis4 = (-1.0_kreal - x + 5.0_kreal *x**2 + 5.0_kreal*x**3)/8.0_kreal
end function

end module
