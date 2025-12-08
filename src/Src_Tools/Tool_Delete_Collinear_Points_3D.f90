 
subroutine Tool_Delete_Collinear_Points_3D(input_points, input_point_num, output_points, output_point_num)
use Global_Float_Type
implicit none
integer, intent(in) :: input_point_num
integer :: output_point_num
real(kind=FT), dimension(input_point_num, 1:3), intent(in) :: input_points
real(kind=FT), dimension(input_point_num, 1:3), intent(out) :: output_points
integer :: i, c_count
logical :: Yes_Three_Points_Collinear
    
real(kind=FT) P1(3),P2(3),P3(3)

output_point_num = 1
c_count = 1


output_points(1, :) = input_points(1, :)

do i = 2, input_point_num-1


    P1(1:3) = input_points(i-1,1:3)
    P2(1:3) = input_points(i,1:3)
    P3(1:3) = input_points(i+1,1:3)
    
    call Tool_Yes_Three_Points_Collinear(p1, p2, p3, Yes_Three_Points_Collinear)
    
    
    if(Yes_Three_Points_Collinear .eqv. .False.) then
        c_count = c_count + 1
        output_points(c_count, :) = input_points(i, :)
    endif
    
end do

output_points(c_count+1, :) = input_points(input_point_num, :)

output_point_num = c_count + 1
end subroutine Tool_Delete_Collinear_Points_3D
