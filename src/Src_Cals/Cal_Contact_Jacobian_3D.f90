 
subroutine Cal_Contact_Jacobian_3D(iter,ifra,Counter_Iter,i_NR_P, &
     c_Total_Freedom,c_num_freeDOF,ori_globalK, &
     freeDOF,Kn,Kn_Gauss,Kt1_Gauss,Kt2_Gauss, &
     CT_State_Gauss,CT_Jacobian)


use Global_Float_Type
use Global_Common
use Global_Model
use Global_Crack_Common
use Global_Crack_3D
use Global_HF
use Global_Material

implicit none
integer, intent(in)::iter,ifra,Counter_Iter,i_NR_P
integer, intent(in)::c_Total_Freedom,c_num_freeDOF
real(kind=FT),intent(in)::ori_globalK(c_Total_Freedom,c_Total_Freedom),Kn
real(kind=FT),intent(inout)::  &
     Kt1_Gauss(num_Crack,Max_Max_N_FluEl_3D),     &
     Kt2_Gauss(num_Crack,Max_Max_N_FluEl_3D),    &
     Kn_Gauss(num_Crack,Max_Max_N_FluEl_3D)
integer, intent(in)::freeDOF(c_Total_Freedom)
integer,intent(in)::CT_State_Gauss(num_Crack,Max_Max_N_FluEl_3D)
real(kind=FT),intent(out)::CT_Jacobian(c_Total_Freedom,c_Total_Freedom)
integer i_N,i_C
integer c_NN(8)
real(kind=FT) c_kesi,c_yita
real(kind=FT) detJ
integer CT_GAUSS_Elem
integer num_CT_Gauss,i_CT_Elem
integer Local_W(60),num_Local_W,Num_Div_Points  
real(kind=FT) x1,y1,x2,y2,ori_CT_elem
real(kind=FT) CT_Length,tem_n(8),N_W(3,60),N(3,24)
real(kind=FT) D_ep(3,3),D_ep_o(3,3),T(3,3)
real(kind=FT) Coor_Tip(2),r,Coor_AB(2,2),Orient_Tip,ele_Sign,c_Sign
real(kind=FT) tem_val1,tem_val2,tem_val3,tem_val4,c_kt,c_kn
integer j_E,c_Adj_Ele,c_Junc_Ele
real(kind=FT) c_Center(3),ori_n(3),Flu_Ele_Area,c_kt1,c_kt2
integer Solid_El
real(kind=FT) Kesi,Yita,Zeta,N_HF(3),Weight_HF
real(kind=FT)BaseLine_A(3),BaseLine_B(3)
real(kind=FT)BaseLine_x_Vec(3),BaseLine_y_Vec(3),BaseLine_z_Vec(3)
real(kind=FT) c_ThetaX,c_ThetaY,c_ThetaZ,c_T_Matrx(3,3)
real(kind=FT) r_Node,Foot_Point_Gauss(3),Foot_Point_Node(3)
real(kind=FT) Node_Point(3),theta_Node
real(kind=FT) BaseLine_Mid(3),Gauss_Coor_Local(3),Node_Coor_Local(3)    
integer ref_elem
real(kind=FT) Tip_T_Matrx(3,3),r_Gauss,theta_Gauss
integer c_Cr_Location
CT_Jacobian(1:c_Total_Freedom,1:c_Total_Freedom) = ZR
CT_Jacobian = ori_globalK


do i_C=1,num_Crack
  do i_CT_Elem = 1,Cracks_FluidEle_num_3D(i_C)
      c_Center = Cracks_FluidEle_Centroid_3D(i_C)%row(i_CT_Elem,1:3)
      ori_n    = Cracks_FluidEle_Vector_3D(i_C)%row(i_CT_Elem,1:3)
      T=Cracks_FluidEle_LCS_T_3D(i_C)%row(i_CT_Elem,1:3,1:3)
      Solid_El = Cracks_FluidEle_EleNum_3D(i_C)%row(i_CT_Elem)     
      Flu_Ele_Area = Cracks_FluidEle_Area_3D(i_C)%row(i_CT_Elem)  

      if(CT_State_Gauss(i_C,i_CT_Elem) /= 0)then
          c_kt1 = ZR
          c_kt2 = ZR
          c_kn = Kn_Gauss(i_C,i_CT_Elem)
          
          
          D_ep_o(1,1:3) = [c_kt1,    ZR,       ZR]
          D_ep_o(2,1:3) = [ZR,    c_kt2,       ZR]
          D_ep_o(3,1:3) = [ZR,       ZR,     c_kn]
          
          call Cal_KesiYita_by_Coor_3D(c_Center,Solid_El,Kesi,Yita,Zeta)
          c_NN(1:8)    = G_NN(1:8,Solid_El)
          call Cal_N_3D(Kesi,Yita,Zeta,N)     
          tem_n(1) = N(1,1);    tem_n(2) = N(1,4)
          tem_n(3) = N(1,7);    tem_n(4) = N(1,10)      
          tem_n(5) = N(1,13);   tem_n(6) = N(1,16)
          tem_n(7) = N(1,19);   tem_n(8) = N(1,22)    
          N_HF(1:3) = ONE/THR      
          D_ep = MATMUL(MATMUL(transpose(T),D_ep_o),T)
          num_Local_W = 0
          N_W(1:3,1:60) =ZR
          do i_N = 1,8
              if(Enriched_Node_Type_3D(c_NN(i_N),i_C) ==2)then
                  num_Local_W = num_Local_W+1
                  Local_W(3*num_Local_W-2)=3*c_POS_3D(c_NN(i_N),i_C)-2
                  Local_W(3*num_Local_W-1)=3*c_POS_3D(c_NN(i_N),i_C)-1     
                  Local_W(3*num_Local_W)  =3*c_POS_3D(c_NN(i_N),i_C)
                  N_W(1,3*num_Local_W-2)   = TWO*tem_n(i_N)
                  N_W(2,3*num_Local_W-1)   = TWO*tem_n(i_N)
                  N_W(3,3*num_Local_W)     = TWO*tem_n(i_N)     
              elseif(Enriched_Node_Type_3D(c_NN(i_N),i_C) ==3)then

              elseif(Enriched_Node_Type_3D(c_NN(i_N),i_C)==1) then
                  
                  if ((Elem_Type_3D(Solid_El,i_C).eq.1 )) then
                      ref_elem=Solid_El
                  else
                      ref_elem=Ele_Num_Tip_Enriched_Node_3D(c_NN(i_N))%row(i_C)
                  endif                
   
                  call Vector_Location_Int_v2(Solid_El_Max_num_Crs,&
                         Solid_El_Crs(ref_elem,1:Solid_El_Max_num_Crs),i_C,c_Cr_Location)

                  BaseLine_A = Solid_El_Tip_BaseLine(ref_elem)%row(c_Cr_Location,1,1:3)
                  BaseLine_B = Solid_El_Tip_BaseLine(ref_elem)%row(c_Cr_Location,2,1:3)
                  BaseLine_Mid = (BaseLine_A+BaseLine_B)/TWO
                  BaseLine_x_Vec=Solid_El_Tip_BaseLine_x_Vec(ref_elem)%row(c_Cr_Location,1:3)
                  BaseLine_y_Vec=Solid_El_Tip_BaseLine_y_Vec(ref_elem)%row(c_Cr_Location,1:3)
                  BaseLine_z_Vec=Solid_El_Tip_BaseLine_z_Vec(ref_elem)%row(c_Cr_Location,1:3)
                  Tip_T_Matrx(1:3,1:3)=Solid_El_Tip_BaseLine_T_Matrix(ref_elem)%row(c_Cr_Location,1:3,1:3)

                  Gauss_Coor_Local= MATMUL(Tip_T_Matrx,c_Center-BaseLine_Mid)
                  r_Gauss = sqrt(Gauss_Coor_Local(1)**2 + Gauss_Coor_Local(2)**2)
                  theta_Gauss = atan2(Gauss_Coor_Local(2),Gauss_Coor_Local(1))    

                  num_Local_W = num_Local_W+1
                  Local_W(3*num_Local_W-2)=3*c_POS_3D(c_NN(i_N),i_C)-2
                  Local_W(3*num_Local_W-1)=3*c_POS_3D(c_NN(i_N),i_C)-1     
                  Local_W(3*num_Local_W)  =3*c_POS_3D(c_NN(i_N),i_C) 

                  N_W(1,3*num_Local_W-2)= TWO*sqrt(r_Gauss)*tem_n(i_N)
                  N_W(2,3*num_Local_W-1)= TWO*sqrt(r_Gauss)*tem_n(i_N)
                  N_W(3,3*num_Local_W)  = TWO*sqrt(r_Gauss)*tem_n(i_N)     
              endif
          enddo
          
          CT_Jacobian(Local_W(1:3*num_Local_W),Local_W(1:3*num_Local_W))                &
             = CT_Jacobian(Local_W(1:3*num_Local_W),Local_W(1:3*num_Local_W)) +         &
          MATMUL(MATMUL(transpose(N_W(1:3,1:3*num_Local_W)),D_ep),N_W(1:3,1:3*num_Local_W))*Flu_Ele_Area
      endif
  enddo
enddo


return 
end SUBROUTINE Cal_Contact_Jacobian_3D         
