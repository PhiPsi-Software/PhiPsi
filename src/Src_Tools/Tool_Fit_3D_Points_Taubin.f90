 
subroutine Tool_Fit_3D_Points_Taubin(i_C,num_Point,In_Points,Out_Points)

use Global_Float_Type     
use Global_Model

implicit none
integer,intent(in)::i_C,num_Point
real(kind=FT),intent(in)::In_Points(num_Point,3)
real(kind=FT),intent(out)::Out_Points(num_Point,3)
real(kind=FT) point(num_Point+1,3)
integer n_Loops,i_Loop,i_point
real(kind=FT) lambda,mu,c_point(3),pre_point(3),nex_point(3)
real(kind=FT) pre_L,nex_L,weight_pre,weight_nex
real(kind=FT) L1_vector(num_Point,3),L2_vector(num_Point,3)


n_Loops =    20


lambda  =  0.330D0
mu      = -0.331D0

point(1:num_Point,1:3) = In_Points(1:num_Point,1:3)
point(num_point+1,1:3) = In_Points(1,1:3)
L1_vector(1:num_Point,1:3) = ZR
L2_vector(1:num_Point,1:3) = ZR

      
do i_Loop=1,n_Loops
    do i_point=1,num_point
        c_point = point(i_point,1:3)
        nex_point = point(i_point+1,1:3)
        if (i_point>=2) then
            pre_point = point(i_point-1,1:3)
        else
            pre_point = point(num_point,1:3)
        endif
		

        pre_L = sqrt((pre_point(1)-c_point(1))**2+(pre_point(2)-c_point(2))**2+(pre_point(3)-c_point(3))**2)
        nex_L = sqrt((nex_point(1)-c_point(1))**2+(nex_point(2)-c_point(2))**2+(nex_point(3)-c_point(3))**2)
        weight_pre = ONE/pre_L
        weight_nex = ONE/nex_L
        L1_vector(i_point,1:3) = (weight_pre*pre_point+weight_nex*nex_point)/(weight_pre+weight_nex)- c_point
    enddo
    
    do i_point=1,num_point
        point(i_point,1:3) = point(i_point,1:3) + lambda*L1_vector(i_point,1:3)
    enddo
    point(num_point+1,1:3) = point(1,1:3)
    
    do i_point=1,num_point
        c_point = point(i_point,1:3)
        nex_point = point(i_point+1,1:3)
        if (i_point>=2) then
            pre_point = point(i_point-1,1:3)
        else
            pre_point = point(num_point,1:3)
        endif
		

        pre_L = sqrt((pre_point(1)-c_point(1))**2+(pre_point(2)-c_point(2))**2+(pre_point(3)-c_point(3))**2)
        nex_L = sqrt((nex_point(1)-c_point(1))**2+(nex_point(2)-c_point(2))**2+(nex_point(3)-c_point(3))**2)
        weight_pre = ONE/pre_L
        weight_nex = ONE/nex_L
        L2_vector(i_point,1:3) = (weight_pre*pre_point+weight_nex*nex_point)/(weight_pre+weight_nex)- c_point
		
        
    enddo
    
    do i_point=1,num_point
        point(i_point,1:3) = point(i_point,1:3)  + mu*L2_vector(i_point,1:3)
    enddo
    point(num_point+1,1:3) = point(1,1:3)
    
enddo

Out_Points(1:num_Point,1:3)  = point(1:num_Point,1:3) 


return 
end SUBROUTINE Tool_Fit_3D_Points_Taubin             
