 
subroutine Tool_Get_3D_Rect_by_n_Vector_Center_LVector_LW_Ratio(n_Vector,&
           c_Crack_Center,L,W,LongSide_Vecotr,Out_Coors)

use Global_Float_Type 
implicit none
real(kind=FT),intent(in) :: n_Vector(3),c_Crack_Center(3),L,W
real(kind=FT),intent(in) :: LongSide_Vecotr(3)
real(kind=FT),intent(out) :: Out_Coors(4,3)

real(kind=FT) W_Vector(3),P_A(3),P_B(3)
real(kind=FT) P1(3),P2(3),P3(3),P4(3)
real(kind=FT) LongSide_Vecotr_Norm(3)

LongSide_Vecotr_Norm = LongSide_Vecotr
call Vector_Normalize(3,LongSide_Vecotr_Norm)     

call Vector_Cross_Product_3(n_Vector,LongSide_Vecotr,W_Vector)   
call Vector_Normalize(3,W_Vector)    

P_A = c_Crack_Center + LongSide_Vecotr_Norm*L/TWO

P_B = c_Crack_Center - LongSide_Vecotr_Norm*L/TWO

P1 = P_A - W_Vector*W/TWO

P2 = P1 - LongSide_Vecotr_Norm*L

P3 = P2 + W_Vector*W

P4 = P_A + W_Vector*W/TWO

Out_Coors(1,1:3) = P1
Out_Coors(2,1:3) = P2
Out_Coors(3,1:3) = P3
Out_Coors(4,1:3) = P4

return 
end SUBROUTINE Tool_Get_3D_Rect_by_n_Vector_Center_LVector_LW_Ratio                     
