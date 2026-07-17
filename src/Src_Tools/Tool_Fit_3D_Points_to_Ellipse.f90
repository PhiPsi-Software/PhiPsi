!-----------------------------------------------------------
! Brief: Fit 3D points to an ellipse and snap to it
!
! Parameters:
!   Input:  num_Point   - number of points to fit
!   Input:  In_Points   - input point coordinates
!   Output: Out_Points  - nearest points on the fitted ellipse
!
! Notes:   Approximate semi-axes from min/max distances and plane normal.
!-----------------------------------------------------------

subroutine Tool_Fit_3D_Points_to_Ellipse(num_Point, In_Points,Out_Points)
!     Fit 3D spatial points into an ellipse, and then output the points on the fitted ellipse.
!     Ref: Parametric Equation of a 3D Ellipse: \theory_documents\035 3D Ellipse Parametric Equation-2022-05-17.png
!     Ref: General equation of a 2D ellipse: \theory_documents\035.1 General equation of a 2D ellipse-2022-05-17.png
!     Approach: (1) Calculate the average coordinates to find the center coordinates of the ellipse.
!     (2) Calculate the average distance from the points to the center of the circle, which is the radius.
!     (3) Calculate the outward normal vector of the plane.
!     (4) Loop through each point, find the point on the circle closest to each point, and mark it as the new coordinate.
!     ---------------Reference materials for obtaining out-of-plane normal vectors-----------------
!     Theory Ref:https://www.freesion.com/article/7157945105/
!                 or
!     \theory_documents\026 Optimal Spatial Circle Fitting for 3D Discrete Points and Implementation-2021-10-30.pdf
!                 
!     Added on 2022-05-17.
!     Based on Tool_Fit_3D_Points_to_Circle.f.

!......................
! Variable Declaration
!......................
use Global_Float_Type     

implicit none
integer,intent(in)::num_Point
real(kind=FT),intent(in)::In_Points(num_Point,3)
real(kind=FT),intent(out)::Out_Points(num_Point,3)
real(kind=FT) Center(3),n_Vector(3)
real(kind=FT) :: Tool_Function_2Point_Dis_3D
real(kind=FT) :: Nearest_Point(3)
real(kind=FT) all_Dis(num_Point),vector_U(3),vector_V(3)
real(kind=FT) a,b
integer :: Farthest_Point

integer :: i_P


Center(1) = sum(In_Points(1:num_Point,1))/dble(num_Point)
Center(2) = sum(In_Points(1:num_Point,2))/dble(num_Point)
Center(3) = sum(In_Points(1:num_Point,3))/dble(num_Point)

all_Dis(1:num_Point) = ZR
do i_P = 1,num_Point
    all_Dis(i_P) = Tool_Function_2Point_Dis_3D(In_Points(i_P,1:3), Center)
enddo

a = maxval(all_Dis)
b = minval(all_Dis)

Farthest_Point = maxloc(all_Dis,1)
vector_U = In_Points(Farthest_Point,1:3)-Center(1:3)
call Vector_Normalize(3,vector_U)   

call Tool_Get_Normal_Vector_of_Points_3D( In_Points,num_Point,n_Vector)

call Vector_Normalize(3,n_Vector)   
call Vector_Cross_Product_3(n_Vector,vector_U,vector_V)   
call Vector_Normalize(3,vector_V)   


do i_P = 1,num_Point
    call Tool_Nearest_Point_from_Point_to_Ellipse_3D( In_Points(i_P,1:3),Center,a,b,vector_U,vector_V, Nearest_Point)
    Out_Points(i_P,1:3) = Nearest_Point(1:3)
enddo

end subroutine Tool_Fit_3D_Points_to_Ellipse             
