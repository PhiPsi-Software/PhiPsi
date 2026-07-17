!-----------------------------------------------------------
! Brief: Compute a 2D point C at distance L from A along AB.
!
! Parameters:
!   Input:  A     - Start point (2-vector)
!   Input:  B     - Reference point (2-vector)
!   Input:  L     - Distance from A to C
!   Output: Out_C - Resulting point C coordinates (2-vector)
!   Output: Statu - Status flag: -1 if C is between A and B,
!                   1 if beyond B, 0 if C coincides with B
!   Output: L_BC  - Resulting distance from C to B
!
! Notes:   The status flag distinguishes whether the new point
!   falls short of, exactly on, or past the second reference.
!-----------------------------------------------------------

subroutine Tool_Point_Giv_2P_and_L(A,B,L,Out_C,Statu,L_BC)
!     Calculate the coordinates of the point, where point C is at a distance L from A, along the direction of AB.

!     case1:Statu = -1
!  
!     *-------------------*----------*
!     A                   C          B

!     case2:Statu = 1
!  
!     *---------*---------*
!     A         B         C

!     case3: Status = 0, point B and point C coincide
!  
!     *-------------------*
!     A                  B(C)

use Global_Float_Type     
implicit none
real(kind=FT),intent(in)::A(2),B(2),L
real(kind=FT),intent(out)::Out_C(2),L_BC
integer,intent(out)::Statu
real(kind=FT) L_AB,theta_AB

L_AB = sqrt((A(1)-B(1))**2+(A(2)-B(2))**2)
if(L_AB > L)then
    Statu = -1
elseif(L_AB < L)then
    Statu = 1
else
    Statu = 0
endif
theta_AB = atan2(B(2)-A(2),B(1)-A(1))

Out_C(1) = A(1)+L*cos(theta_AB)
Out_C(2) = A(2)+L*sin(theta_AB)

L_BC = sqrt((Out_C(1)-B(1))**2+(Out_C(2)-B(2))**2)

return 
end subroutine Tool_Point_Giv_2P_and_L                          
