 
SUBROUTINE Tool_Fit_Least_Square_Straight_Line(x, y, n, k, b)
use Global_Float_Type
implicit none
integer, intent(in) :: n
real(kind=FT), intent(in) :: x(n), y(n)
real(kind=FT), intent(out) :: k, b
real(kind=FT) :: sum_x, sum_y, sum_xy, sum_x2
integer :: i

sum_x = ZR
sum_y = ZR
sum_xy = ZR
sum_x2 = ZR

do i = 1, n
    sum_x = sum_x + x(i)
    sum_y = sum_y + y(i)
    sum_xy = sum_xy + x(i) * y(i)
    sum_x2 = sum_x2 + x(i) * x(i)
end do

k = (n * sum_xy - sum_x * sum_y) / (n * sum_x2 - sum_x * sum_x)
b = (sum_y - k * sum_x) / n

RETURN
END SUBROUTINE Tool_Fit_Least_Square_Straight_Line
