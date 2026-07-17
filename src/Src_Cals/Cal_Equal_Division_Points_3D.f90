!-----------------------------------------------------------
! Brief: Compute equally spaced 3D points along the line segment AB.
!
! Parameters:
!   Input:  A                   - endpoint A of the segment (x,y,z)
!           B                   - endpoint B of the segment (x,y,z)
!           Num_Div             - number of equal divisions
!           Yes_Include_Endpoint - true to include A and B in the output
!   Output: Div_Points          - array of division points
!           Num_Div_Points      - number of points actually written
!
! Notes:   Points are ordered from A to B.
!-----------------------------------------------------------

subroutine Cal_Equal_Division_Points_3D(A,B, Num_Div,Yes_Include_Endpoint, Div_Points,Num_Div_Points)
!Get the equal diversion points of line AB in three-dimension.
!Diversion point are arranged from A to B.   
use Global_Float_Type
use Global_Crack

implicit none

real(kind=FT),intent(in)::A(3),B(3)
integer,intent(in)::Num_Div
logical,intent(in)::Yes_Include_Endpoint
real(kind=FT),intent(out)::Div_Points(Max_Num_Seg_CalP,3)
integer,intent(out)::Num_Div_Points

integer :: i
real(kind=FT) a_x,a_y,a_z,b_x,b_y,b_z

a_x = A(1)
a_y = A(2)
a_z = A(3)
b_x = B(1)
b_y = B(2)
b_z = B(3)

Num_Div_Points = 0
if (Yes_Include_Endpoint.eqv..False.) then
    do i = 1,Num_Div-1
        Div_Points(i,1) = (i*b_x+(Num_Div-i)*a_x)/Num_Div
        Div_Points(i,2) = (i*b_y+(Num_Div-i)*a_y)/Num_Div
        Div_Points(i,3) = (i*b_z+(Num_Div-i)*a_z)/Num_Div
    end do
    Num_Div_Points = Num_Div-1
else
    Div_Points(1,:)= A
    Num_Div_Points = 1
    do i = 1,Num_Div-1
        Num_Div_Points  = Num_Div_Points + 1
        Div_Points(Num_Div_Points,1)= (i*b_x+(Num_Div-i)*a_x)/Num_Div
        Div_Points(Num_Div_Points,2)= (i*b_y+(Num_Div-i)*a_y)/Num_Div
        Div_Points(Num_Div_Points,3)= (i*b_z+(Num_Div-i)*a_z)/Num_Div
    end do 
    Num_Div_Points  = Num_Div_Points + 1
    Div_Points(Num_Div_Points,:) = B
end if

return 
end subroutine Cal_Equal_Division_Points_3D              
