 
subroutine Tool_Dis_Point_to_3D_Tri_only_Dis(Point,Tri_P1,Tri_P2,Tri_P3, Distance)

use Global_Float_Type     
use Global_Common
implicit none
real(kind=FT),intent(in)::Point(3),Tri_P1(3),Tri_P2(3),Tri_P3(3)
real(kind=FT),intent(out):: Distance
real(kind=FT) P1P2(3),P1P3(3),P1P0(3),Np(3),cosa
real(kind=FT) abs_P1P0,abs_Np
real(kind=FT) B_i(3),C_i(3),X_i(3)
real(kind=FT) c_DET
integer c_Orient

Distance   =  ZR

P1P2 = Tri_P2 - Tri_P1
P1P3 = Tri_P3 - Tri_P1
Np(1) = P1P2(2) * P1P3(3) - P1P2(3) * P1P3(2)
Np(2) = P1P2(3) * P1P3(1) - P1P2(1) * P1P3(3)
Np(3) = P1P2(1) * P1P3(2) - P1P2(2) * P1P3(1)

P1P0     = Point - Tri_P1
abs_P1P0 = sqrt(P1P0(1)**2 + P1P0(2)**2 + P1P0(3)**2)
if(abs_P1P0<=Tol_10) abs_P1P0=Tol_10
abs_Np   = sqrt(Np(1)**2   + Np(2)**2   + Np(3)**2)
cosa     = dot_product(P1P0,Np)/abs_P1P0/abs_Np 


B_i = Tri_P2-Tri_P1
C_i = Tri_P3-Tri_P1
X_i = Point-Tri_P1

c_DET = B_i(1)*C_i(2)*X_i(3) +B_i(2)*C_i(3)*X_i(1) +  &
          C_i(1)*X_i(2)*B_i(3) -B_i(1)*C_i(3)*X_i(2) -  &
          B_i(2)*C_i(1)*X_i(3) -B_i(3)*C_i(2)*X_i(1)

if(c_DET>ZR)then
      c_Orient =  1
elseif(c_DET<ZR)then
      c_Orient = -1
else
      c_Orient =  0
endif


Distance = dble(c_Orient)*abs(abs_P1P0*cosa)
#ifndef Silverfrost
if (isnan(Distance)) then
    print *, '    Error-2023061501 :: Distance is NAN in Tool_Dis_Point_to_3D_Tri_only_Dis.f90!'
    call Warning_Message('S',Keywords_Blank)
endif
#endif  
return 
end SUBROUTINE Tool_Dis_Point_to_3D_Tri_only_Dis                  
