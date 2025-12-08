 
subroutine Tool_Yes_Point_in_3D_Plane_Polygon(num_Point,Polygon,Point,&
                                              Yes_Point_in)



use Global_Float_Type 
implicit none
integer,intent(in) :: num_Point
real(kind=FT),intent(in) :: Point(3),Polygon(num_Point,3)
logical,intent(out) :: Yes_Point_in
real(kind=FT) Area_Polygon,Sum_Area,c_Area
real(kind=FT) Point_1(3),Point_2(3),Point_3(3)
integer num_lines,i_Lines

Yes_Point_in = .False.

call Tool_Area_3D_Plane_Polygon(num_Point,Polygon(1:num_Point,1:3),Area_Polygon)

Sum_Area  = ZR
num_lines = num_Point
Point_1  = Point
do i_Lines = 1,num_lines
    Point_2 = Polygon(i_Lines,1:3)
    if(i_Lines<num_lines) then
        Point_3 = Polygon(i_Lines+1,1:3)
    else
        Point_3 = Polygon(1,1:3)
    endif
    call Tool_Area_Tri_3D(Point_1,Point_2,Point_3,c_Area)
    Sum_Area = Sum_Area + c_Area
enddo


if(abs(Sum_Area-Area_Polygon)<=Tol_9) then
    Yes_Point_in = .True.
    return
endif

return 
end SUBROUTINE Tool_Yes_Point_in_3D_Plane_Polygon                     
