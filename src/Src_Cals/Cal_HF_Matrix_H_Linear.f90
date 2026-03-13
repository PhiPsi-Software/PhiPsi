!     ================================================= !
!             ____  _       _   ____  _____   _         !
!            |  _ \| |     |_| |  _ \|  ___| |_|        !
!            | |_) | |___   _  | |_) | |___   _         !
!            |  _ /|  _  | | | |  _ /|___  | | |        !
!            | |   | | | | | | | |    ___| | | |        !
!            |_|   |_| |_| |_| |_|   |_____| |_|        !
!     ================================================= !
!     PhiPsi:     a general-purpose computational       !
!                 mechanics program written in Fortran. !
!     Website:    http://phipsi.top                     !
!     Author:     Shi Fang, Huaiyin Institute of        !
!                 Technology, Huaian, JiangSu, China    !
!     Email:      shifang@hyit.edu.cn                   !
!     ------------------------------------------------- !
!     Please cite the following papers:                 !
!     (1)Shi F., Lin C. Modeling fluid-driven           !
!        propagation of 3D complex crossing fractures   !
!        with the extended finite element method.       !
!        Computers and Geotechnics, 2024, 172, 106482.  !
!     (2)Shi F., Wang D., Li H. An XFEM-based approach  !
!        for 3D hydraulic fracturing simulation         !
!        considering crack front segmentation. Journal  !
!        of Petroleum Science and Engineering, 2022,    !
!        214, 110518.                                   !
!     (3)Shi F., Wang D., Yang Q. An XFEM-based         !
!        numerical strategy to model three-dimensional  !
!        fracture propagation regarding crack front     !
!        segmentation. Theoretical and Applied Fracture !
!        Mechanics, 2022, 118, 103250.                  !
!     (4)Shi F., Liu J. A fully coupled hydromechanical !
!        XFEM model for the simulation of 3D non-planar !
!        fluid-driven fracture propagation. Computers   !
!        and Geotechnics, 2021, 132: 103971.            !
!     (5)Shi F., Wang X.L., Liu C., Liu H., Wu H.A. An  !
!        XFEM-based method with reduction technique     !
!        for modeling hydraulic fracture propagation    !
!        in formations containing frictional natural    !
!        fractures. Engineering Fracture Mechanics,     !
!        2017, 173: 64-90.                              !
!     ------------------------------------------------- !
 
!subroutine Cal_HF_Matrix_H_Linear(ifra,Counter,Matrix_H,Last_Cr_CalP_Aper,total_time)

subroutine Cal_HF_Matrix_H_Linear(ifra,Counter,Matrix_H,Last_Cr_CalP_Aper,total_time)
! Matrix H used in the Newton-Raphson iteration algorithm for hydraulic fracturing.
! For linear HF elements
      
!----------------------------
!Read public variable module
!----------------------------
use Global_Float_Type
use Global_Crack
use Global_Crack_Common
use Global_HF

implicit none

integer,intent(in)::ifra,Counter
real(kind=FT),intent(in)::total_time
real(kind=FT),intent(in)::Last_Cr_CalP_Aper(Max_Num_Cr,Max_Num_Cr_CalP)
real(kind=FT),intent(out)::Matrix_H(num_Tol_CalP_Water,num_Tol_CalP_Water)

      
integer i_C,i_HF_Elem,i_HF_GAUSS    
integer Num_Div_Points 

real(kind=FT) kesi_HF(4),Weight_HF(4)
real(kind=FT) N_HF(2)
real(kind=FT) P_N_HF(2)
real(kind=FT) HF_Length
integer Local_P(2)                          
real(kind=FT) detJ,w1,w2,w,k
real(kind=FT) tem_H(2,2),tem(2,2)
integer num_HF_Gauss,E_count,HF_node_1,HF_node_2
real(kind=FT) Old_c1,Old_c2,Old_c,Old_Vis
integer Num_Cr_has_Fluid
real(kind=FT) Length_Crack,Length_Out,c_Time,c_CalP_x,c_CalP_y
integer i
real(kind=FT), allocatable :: Last_Cr_CalP_Conc(:,:)
num_HF_Gauss = 2
kesi_HF(1)  = -0.577350269189626D0
kesi_HF(2)  =  0.577350269189626D0
Weight_HF(1)=  ONE
Weight_HF(2)=  ONE

allocate(Last_Cr_CalP_Conc(Max_Num_Cr,Max_Num_Cr_CalP))
Last_Cr_CalP_Conc(1:Max_Num_Cr,1:Max_Num_Cr_CalP) =ZR
if(ifra==1)then
elseif(ifra > 1)then
     if (Key_Propp_Trans==1) then
        Last_Cr_CalP_Conc = Map_L_Cracks_CalP_Conc
     endif
endif
  
Matrix_H(1:num_Tol_CalP_Water,1:num_Tol_CalP_Water) = ZR
 E_count = 0
 Num_Cr_has_Fluid = 0
 do i_C=1,num_Crack
     if (Cracks_HF_State(i_C) == 1) then   
         Num_Cr_has_Fluid = Num_Cr_has_Fluid+1
         Num_Div_Points = Cracks_CalP_Num(i_C)
         do i_HF_Elem = 1,Num_Div_Points - 1
             
             E_count = E_count + 1
             tem_H(1:2,1:2) = ZR   
             HF_node_1 = (Num_Cr_has_Fluid-1) + E_count            
             HF_node_2 = (Num_Cr_has_Fluid-1) + E_count + 1     
             
             Local_P(1) =  HF_node_1
             Local_P(2) =  HF_node_2
             HF_Length = Cracks_HF_Ele_L(i_C,i_HF_Elem)
             detJ   = HF_Length/TWO
             w1 = Last_Cr_CalP_Aper(i_C,i_HF_Elem)
             w2 = Last_Cr_CalP_Aper(i_C,i_HF_Elem+1)


             Old_Vis = Viscosity
             
             if(Key_Visco_Type==1)then
                 Old_Vis = Viscosity
             elseif(Key_Visco_Type==2)then
                 Old_c1 = Last_Cr_CalP_Conc(i_C,i_HF_Elem)
                 Old_c2 = Last_Cr_CalP_Conc(i_C,i_HF_Elem+1)
                 Old_c =  HLF*(Old_c1 + Old_c2)
                 Old_Vis = Viscosity*(ONE-Old_c/Max_c)**(-Viscosity_Par_m)
             endif       
             
             if(Old_Vis > Visco_Zoom_Factor*Viscosity)then
                 Old_Vis=Visco_Zoom_Factor*Viscosity
             endif
             
             if(Key_Visco_Type==3)then
                 call Tool_Length_of_Crack(i_C,Length_Crack)
                 c_CalP_x = Cracks_CalP_Coors(i_C,i_HF_Elem,1) 
                 c_CalP_y = Cracks_CalP_Coors(i_C,i_HF_Elem,2) 
                 call Tool_Length_CalP_to_InjP(i_C,Length_Out,c_CalP_x,c_CalP_y,i_HF_Elem)
                 c_Time=(total_time/Length_Crack**2)*Length_Out**2
                 
                 if(c_Time<=1.0D-5)then
                     c_Time = 1.0D-5
                 endif
                 Old_Vis =Viscosity_td_m*c_Time**Viscosity_td_n
             end if 
             
             do i_HF_GAUSS  = 1,num_HF_Gauss
                 N_HF(1) = (ONE-kesi_HF(i_HF_GAUSS))/TWO      
                 N_HF(2) = (ONE+kesi_HF(i_HF_GAUSS))/TWO   
                 
                 w = N_HF(1)*w1 + N_HF(2)*w2
                 k = w**3/12.0D0/Old_Vis
                 
                 P_N_HF(1) = -ONE/HF_Length  
                 P_N_HF(2) =  ONE/HF_Length  
                 
                 tem(1,1) = P_N_HF(1)*P_N_HF(1)
                 tem(1,2) = P_N_HF(1)*P_N_HF(2)
                 tem(2,1) = P_N_HF(2)*P_N_HF(1)
                 tem(2,2) = P_N_HF(2)*P_N_HF(2)
                 
                 tem_H = tem_H + tem*Weight_HF(i_HF_GAUSS)*detJ*k  

             end do
      
             Matrix_H(Local_P(1:2),Local_P(1:2)) =  Matrix_H(Local_P(1:2),Local_P(1:2)) +tem_H
         end do  
     end if
 end do  
deallocate(Last_Cr_CalP_Conc)      
return 
end SUBROUTINE Cal_HF_Matrix_H_Linear        
