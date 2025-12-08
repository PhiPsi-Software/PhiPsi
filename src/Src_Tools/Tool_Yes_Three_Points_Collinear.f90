 
SUBROUTINE Tool_Yes_Three_Points_Collinear(p1, p2, p3, res)
use Global_Float_Type
IMPLICIT NONE
REAL(kind=FT), DIMENSION(3), INTENT(IN) :: p1, p2, p3
LOGICAL :: res

REAL(kind=FT) :: cross_product(3)


res = .False.

call Vector_Cross_Product_3(p2 - p1,p3 - p1,cross_product)   

res = ALL(ABS(cross_product) < 1.0D-6)
END SUBROUTINE  Tool_Yes_Three_Points_Collinear