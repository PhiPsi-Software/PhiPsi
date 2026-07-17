!-----------------------------------------------------------
! Brief: Computes the area of a 3D triangle from its three vertex coordinates.
!
! Parameters:
!   Input:  A, B, C - 3D vertex coordinates
!   Output: area    - triangle area
!
! Notes:   Builds the two edge vectors AB and AC and returns half of the
!          magnitude of their cross product.
!-----------------------------------------------------------

    subroutine Tool_Area_Tri_3D(A,B,C,area)

    ! Given the coordinates of the vertices, calculate the area of a triangle in three-dimensional
    ! space.
    ! See my notes V3-P139 for details.

    use Global_Float_Type
    IMPLICIT NONE
    real(kind=FT),intent(in)::A(3),B(3),C(3)
    real(kind=FT),intent(out)::area
    real(kind=FT) Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz
    real(kind=FT) x1,x2,x3,y1,y2,y3

    Ax = A(1)
    Ay = A(2)
    Az = A(3)
    Bx = B(1)
    By = B(2)
    Bz = B(3)
    Cx = C(1)
    Cy = C(2)
    Cz = C(3)

    x1 = Bx - Ax;x2 = By - Ay;x3 = Bz - Az

    y1 = Cx - Ax;y2 = Cy - Ay;y3 = Cz - Az

area = HLF*sqrt((x2*y3-x3*y2)**2 + (x3*y1-x1*y3)**2 + (x1*y2-x2*y1)**2)

    RETURN
    END subroutine Tool_Area_Tri_3D
