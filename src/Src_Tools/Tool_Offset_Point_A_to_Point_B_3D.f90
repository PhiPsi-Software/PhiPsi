!-----------------------------------------------------------
! Brief: Shift point A toward point B by a signed distance in 3D.
!
! Parameters:
!   Input:  A       - Starting point coordinates (3-vector)
!   Input:  B       - Target point coordinates (3-vector)
!   Input:  delta_L - Signed offset distance along AB
!   Output: new_A   - New location of point A after offset
!
! Notes:   Uses direction cosines of the AB vector to compute
!   the offset coordinates; positive delta_L moves A toward B.
!-----------------------------------------------------------

subroutine Tool_Offset_Point_A_to_Point_B_3D(A,B,delta_L,new_A)
! Offset point A towards B by delta_L.
!2022-05-05.
!
!                A                       B
!
!                o-----------------------o
!
!                delta_L > 0:
!
!                |<---delta_L--->|
!                                o-------o
!              
!              
! Reference: http://wenku.baidu.com/view/2613ac1a227916888486d7cb.html, Vectors in Three-Dimensional
! Space 2011

use Global_Float_Type      
implicit none
real(kind=FT),intent(in) ::A(3),B(3),delta_L
real(kind=FT),intent(out)::new_A(3)
real(kind=FT) Vector_AB(3),norm_Vector_AB
real(kind=FT) theta1,theta2,theta3

Vector_AB = B(1:3) - A(1:3) 

norm_Vector_AB = sqrt(Vector_AB(1)**2 + Vector_AB(2)**2 + Vector_AB(3)**2)
theta1 = acos(Vector_AB(1)/norm_Vector_AB)
theta2 = acos(Vector_AB(2)/norm_Vector_AB) 
theta3 = acos(Vector_AB(3)/norm_Vector_AB) 

new_A(1)= A(1)+delta_L*cos(theta1)
new_A(2)= A(2)+delta_L*cos(theta2)
new_A(3)= A(3)+delta_L*cos(theta3)

return 
end subroutine Tool_Offset_Point_A_to_Point_B_3D                         
