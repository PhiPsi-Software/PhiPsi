!-----------------------------------------------------------
! Brief: Solves for the coefficients (a,b,c) of a parabola passing through three 2D points.
!
! Parameters:
!   Input:  Point_1(2),Point_2(2),Point_3(2) - three collinear-x points on the parabola
!   Output: a,b,c - coefficients of y = a*x^2 + b*x + c
!
! Notes:   Builds the 3x3 Vandermonde system K*D = F, inverts K via
!          Matrix_Inverse, and returns the solution vector.
!-----------------------------------------------------------

subroutine Tool_abc_of_Parabola(Point_1,Point_2,Point_3,a,b,c)
! Calculate the coefficients a, b, c of the parabola
! The input variables are the x and y values of three points arranged sequentially on a parabolic curve.
! Calculate the equation of a parabola: ax**2 + bx + c = 0
use Global_Float_Type
implicit none
real(kind=FT),intent(in):: Point_1(2),Point_2(2),Point_3(2)
real(kind=FT),intent(out):: a,b,c

real(kind=FT) D(3),K(3,3),inv_K(3,3),F(3),x1,x2,x3,y1,y2,y3

x1 = Point_1(1); y1 = Point_1(2)
x2 = Point_2(1); y2 = Point_2(2)
x3 = Point_3(1); y3 = Point_3(2)

K(1,1:3) = [x1**2, x1, ONE]
K(2,1:3) = [x2**2, x2, ONE]
K(3,1:3) = [x3**2, x3, ONE]


F(1:3)   = [y1,y2,y3]

call Matrix_Inverse(K,inv_K,3)

D = MATMUL(inv_K,F)

a = D(1)
b = D(2)
c = D(3)


return
end SUBROUTINE Tool_abc_of_Parabola
