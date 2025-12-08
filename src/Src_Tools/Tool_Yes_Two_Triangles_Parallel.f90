 
subroutine Tool_Yes_Two_Triangles_Parallel(Tri_1,Tri_2,Yes_Parallel)   

use Global_Float_Type
use Global_Common
implicit none
real(kind=FT),intent(in)::Tri_1(3,3),Tri_2(3,3)
logical,intent(out)::Yes_Parallel            
real(kind=FT) Tri_1_n_Vector(3), Tri_2_n_Vector(3)
real(kind=FT) Tri_1_a_Vector(3), Tri_1_b_Vector(3)
real(kind=FT) Tri_2_a_Vector(3), Tri_2_b_Vector(3)
real(kind=FT) c_angle
Yes_Parallel  = .False.
     
Tri_1_a_Vector(1:3) = Tri_1(2,1:3) - Tri_1(1,1:3)
Tri_1_b_Vector(1:3) = Tri_1(3,1:3) - Tri_1(1,1:3)
call Vector_Cross_Product_3(Tri_1_a_Vector,Tri_1_b_Vector,Tri_1_n_Vector)   

Tri_2_a_Vector(1:3) = Tri_2(2,1:3) - Tri_2(1,1:3)
Tri_2_b_Vector(1:3) = Tri_2(3,1:3) - Tri_2(1,1:3)
call Vector_Cross_Product_3(Tri_2_a_Vector,Tri_2_b_Vector,Tri_2_n_Vector)   

call Tool_Angle_of_Vectors_a_and_b_3D(Tri_1_n_Vector,Tri_2_n_Vector,c_angle,2)

if(abs(c_angle)<=Tol_8) then
    Yes_Parallel  = .True.
    return
endif

if(abs(c_angle-pi)<=Tol_8) then
    Yes_Parallel  = .True.
    return
endif

return 
end SUBROUTINE Tool_Yes_Two_Triangles_Parallel       
