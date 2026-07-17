!-----------------------------------------------------------
! Brief: Tests whether three 3D points are collinear.
!
! Parameters:
!   Input:  p1(3),p2(3),p3(3) - the three points to test
!   Output: res - .true. if the three points lie on a common line
!
! Notes:   Computes the cross product (p2-p1) x (p3-p1); if all components
!          are below 1e-6, the points are considered collinear.
!-----------------------------------------------------------

SUBROUTINE Tool_Yes_Three_Points_Collinear(p1, p2, p3, res)
! Determine whether three points are collinear. 2024-05-03.
use Global_Float_Type
IMPLICIT NONE
REAL(kind=FT), DIMENSION(3), INTENT(IN) :: p1, p2, p3
LOGICAL :: res

REAL(kind=FT) :: cross_product(3)


res = .False.

call Vector_Cross_Product_3(p2 - p1,p3 - p1,cross_product)   

res = ALL(ABS(cross_product) < 1.0D-6)
END SUBROUTINE  Tool_Yes_Three_Points_Collinear