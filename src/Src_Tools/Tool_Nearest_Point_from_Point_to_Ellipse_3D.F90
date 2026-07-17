!-----------------------------------------------------------
! Brief: Find the closest point on a 3D ellipse to a given query point
!
! Parameters:
!   Input:  In_Point - Query point in 3D
!   Input:  Center   - Ellipse center
!   Input:  a, b     - Semi-axes of the ellipse
!   Input:  vec_U, vec_V - In-plane unit basis vectors of the ellipse
!   Output: Out_Point - Closest point on the ellipse to the query
!
! Notes:   Minimizes the squared distance with a 1D bracket search using
!          the open-source fmin routine over the parameter t in [0, 2 pi].
!-----------------------------------------------------------

subroutine Tool_Nearest_Point_from_Point_to_Ellipse_3D( In_Point,Center,a,b,vec_U,vec_V,Out_Point)
    ! Calculate the point on the ellipse that is closest to the vertex In_Point.
    ! 2022-05-07.

    use Global_Float_Type
    use fmin_module

    implicit none
    real(kind=FT),intent(in)  :: In_Point(3),Center(3),a,b
    real(kind=FT),intent(in)  :: vec_U(3),vec_V(3)
    real(kind=FT),intent(out) :: Out_Point(3)
    real(kind=FT) :: ax, bx, xmin
    real(kind=FT),parameter :: tol = 1.0D-8
    ax = ZR
    bx = TWO*pi

    xmin = fmin(func,ax,bx,tol)

Out_Point(1:3) = Center(1:3) + a*cos(xmin)*vec_U(1:3)+b*sin(xmin)*vec_V(1:3)


    contains

    function func(x) result(f)
    use Global_Float_Type
    implicit none

    real(kind=FT),intent(in) :: x
    real(kind=FT)            :: f

f = (In_Point(1)-Center(1)-a*cos(x)*vec_U(1)-b*sin(x)*vec_V(1))**2 + &
(In_Point(2)-Center(2)-a*cos(x)*vec_U(2)-b*sin(x)*vec_V(2))**2 + &
(In_Point(3)-Center(3)-a*cos(x)*vec_U(3)-b*sin(x)*vec_V(3))**2

    end function func


    end SUBROUTINE Tool_Nearest_Point_from_Point_to_Ellipse_3D
