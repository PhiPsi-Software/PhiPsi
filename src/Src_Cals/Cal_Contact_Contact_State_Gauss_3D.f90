 
SUBROUTINE Cal_Contact_Contact_State_Gauss_3D(iter,ifra,Counter_Iter,i_NR_P, & 
                c_DISP,Yes_Contact,c_Elem_Conta_Sta,CT_State_Gauss) 

use Global_Float_Type      
use Global_Crack_Common
use Global_Crack_3D
use Global_Model
use Global_Common
use Global_Contact
use Global_HF
use Global_Elem_Area_Vol
use Global_Crack
use omp_lib

implicit none
integer,intent(in)::iter,ifra,Counter_Iter,i_NR_P
real(kind=FT),intent(in)::c_DISP(Total_FD)
logical,intent(out)::Yes_Contact
integer,intent(out)::CT_State_Gauss(num_Crack,Max_Max_N_FluEl_3D)
integer,intent(out)::c_Elem_Conta_Sta(Num_Elem,num_Crack)
integer i_C
integer i_CT_Elem
real(kind=FT) ori_n(3),ori_CT_elem
real(kind=FT) c_Aperture
real(kind=FT) c_Center(3),T(3,3)
integer Solid_El
real(kind=FT) Relative_Disp(3)




CT_State_Gauss(1:num_Crack,1:Max_Max_N_FluEl_3D)  = 0
c_Elem_Conta_Sta(1:Num_Elem,1:num_Crack)   = 0
Yes_Contact = .False.
  
  
  
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i_C,i_CT_Elem,c_Center,  & 
!$OMP            ori_n,T,Solid_El,c_Aperture,Relative_Disp)   & 
!$OMP            SCHEDULE(STATIC) 
do i_C=1,num_Crack
      do i_CT_Elem = 1,Cracks_FluidEle_num_3D(i_C)
          c_Center = Cracks_FluidEle_Centroid_3D(i_C)%row(i_CT_Elem,1:3)
          ori_n    = Cracks_FluidEle_Vector_3D(i_C)%row(i_CT_Elem,1:3)
          Solid_El = Cracks_FluidEle_EleNum_3D(i_C)%row(i_CT_Elem)  
          T = Cracks_FluidEle_LCS_T_3D(i_C)%row(i_CT_Elem,1:3,1:3)
          call Cal_Point_Aperture_3D(i_C,c_Center,c_DISP,ori_n,T,c_Aperture,10)
          if(Key_Crack_Aperture_Method==1) then
              call Cal_Crack_Point_Aperture_3D(c_DISP,i_C,c_Center,Relative_Disp,i_CT_Elem,0,0,0)
              c_Aperture   = (Relative_Disp(1)*ori_n(1)+ Relative_Disp(2)*ori_n(2)+ Relative_Disp(3)*ori_n(3))/1.0D0
          elseif(Key_Crack_Aperture_Method==2) then
              call Cal_Point_Aperture_3D(i_C,c_Center,c_DISP,ori_n,T,c_Aperture,10)
          endif
          if (c_Aperture < Contact_Aper_Tol) then 
              CT_State_Gauss(i_C,i_CT_Elem) = 1
              !$omp critical
              c_Elem_Conta_Sta(Solid_El,i_C) = 1
              !$omp end critical
              Yes_Contact = .True.
          end if    
      enddo
enddo
!$omp end parallel do  

RETURN
END SUBROUTINE Cal_Contact_Contact_State_Gauss_3D
