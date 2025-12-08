 
subroutine Tool_Area_3D_Plane_Polygon(num_Point,Polygon,Area)

use Global_Float_Type 
implicit none
integer,intent(in) :: num_Point
real(kind=FT),intent(in) :: Polygon(num_Point,3)
real(kind=FT),intent(out) :: Area
real(kind=FT) total(3)
real(kind=FT) vi1(3),vi2(3),prod(3)
real(kind=FT) unit_normal_vector(3)
integer i
real(kind=FT) c_Point(3),c_Point_Next(3),c_Point_Next_Next(3)
logical c_Yes
integer Deleted_Points_Flag(num_Point)
integer Num_Deleted_Points,num_New_Points
real(kind=FT),ALLOCATABLE:: Polygon_New(:,:)
integer j


if (num_Point < 3) then
    Area = ZR
    return 
endif

Deleted_Points_Flag(1:num_Point) = 0
do i=1, num_Point
    c_Point = Polygon(i,1:3)
    if (i==num_Point) then
        c_Point_Next = Polygon(1,1:3)
        c_Point_Next_Next = Polygon(2,1:3)
    elseif (i==(num_Point-1)) then
        c_Point_Next = Polygon(num_Point,1:3)
        c_Point_Next_Next = Polygon(1,1:3)
    else
        c_Point_Next = Polygon(i+1,1:3)
        c_Point_Next_Next  = Polygon(i+2,1:3)
    endif
    call Tool_Yes_Point_on_Line_Segment_3D(c_Point,c_Point_Next_Next,c_Point_Next,c_Yes)
    if(c_Yes)then
        if (i==num_Point) then
            Deleted_Points_Flag(1) = 1
        elseif (i==(num_Point-1)) then
            Deleted_Points_Flag(num_Point) = 1
        else
            Deleted_Points_Flag(i+1) = 1
        endif
    endif
enddo
Num_Deleted_Points = sum(Deleted_Points_Flag)
num_New_Points     = num_Point - Num_Deleted_Points

allocate(Polygon_New(num_New_Points,3))
j = 0 
do i=1, num_Point
    if (Deleted_Points_Flag(i)==0) then
        j = j+ 1
        Polygon_New(j,1:3) = Polygon(i,1:3)
    endif
enddo


if(Num_Deleted_Points==0)  then
    total = [ZR, ZR, ZR]
    do i=1, num_Point
        vi1 = Polygon(i,1:3)
        if (i==num_Point) then
            vi2 = Polygon(1,1:3)
        else
            vi2 = Polygon(i+1,1:3)
        endif
        call Vector_Cross_Product_3(vi1, vi2,prod)  
        total(1) = total(1) + prod(1)
        total(2) = total(2) + prod(2)
        total(3) = total(3) + prod(3)
    enddo
    call unit_normal(Polygon(1,1:3),Polygon(2,1:3),Polygon(3,1:3),unit_normal_vector)
    Area = dot_product(total,unit_normal_vector)
    Area = abs(Area/TWO)
else
    total = [ZR, ZR, ZR]
    do i=1, num_New_Points
        vi1 = Polygon_New(i,1:3)
        if (i==num_New_Points) then
            vi2 = Polygon_New(1,1:3)
        else
            vi2 = Polygon_New(i+1,1:3)
        endif
        call Vector_Cross_Product_3(vi1, vi2,prod)  
        total(1) = total(1) + prod(1)
        total(2) = total(2) + prod(2)
        total(3) = total(3) + prod(3)
    enddo
    call unit_normal(Polygon_New(1,1:3),Polygon_New(2,1:3),Polygon_New(3,1:3),unit_normal_vector)
    Area = dot_product(total,unit_normal_vector)
    Area = abs(Area/TWO)
endif

return 
end SUBROUTINE Tool_Area_3D_Plane_Polygon      



subroutine unit_normal(A,B,C,unit_normal_vector)
use Global_Float_Type 
implicit none
real(kind=FT),intent(in) :: A(3), B(3), C(3)
real(kind=FT),intent(out) :: unit_normal_vector(3)
real(kind=FT) x,y,z,Matrix(3,3),magnitude

Matrix(1,1:3) = [ONE,A(2),A(3)] 
Matrix(2,1:3) = [ONE,B(2),B(3)] 
Matrix(3,1:3) = [ONE,C(2),C(3)] 
call Matrix_Det_3x3(Matrix,x) 

Matrix(1,1:3) = [A(1),ONE,A(3)] 
Matrix(2,1:3) = [B(1),ONE,B(3)] 
Matrix(3,1:3) = [C(1),ONE,C(3)] 
call Matrix_Det_3x3(Matrix,y)

Matrix(1,1:3) = [A(1),A(2),ONE] 
Matrix(2,1:3) = [B(1),B(2),ONE] 
Matrix(3,1:3) = [C(1),C(2),ONE] 
call Matrix_Det_3x3(Matrix,z)

magnitude = sqrt((x**2 + y**2 + z**2))

unit_normal_vector= [x/magnitude, y/magnitude, z/magnitude]

return 
end SUBROUTINE unit_normal
              
