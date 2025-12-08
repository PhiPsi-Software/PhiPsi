 
SUBROUTINE Force_Vector_3D(c_Total_FD,isub,Lambda,globalF)


use Global_Float_Type
use Global_Common
use Global_Model
use Global_Elem_Area_Vol
use Global_Material
use Global_Stress
use Global_Crack_Common
use Global_Crack_3D
use Global_Filename
use Global_Surface_Load
use omp_lib
   
implicit none
integer,intent(in)::isub,c_Total_FD
real(kind=FT),intent(in)::Lambda
real(kind=FT),intent(out)::globalF(c_Total_FD)
real(kind=FT) h,c_density
integer i,i_E,i_N,cur_Node,mat_num,c_mat_type,c_node
real(kind=FT) c_thick,c_D(6,6),c_B(6,24)
real(kind=FT) c_X_NODES(8),c_Y_NODES(8),c_Z_NODES(8)
integer c_NN(8) 
integer c_Num_Gauss_Point,i_G
real(kind=FT) c_kesi,c_yita,c_zeta
integer G_Counter
real(kind=FT) c_Stress(6),c_temp(24),detJ
real(kind=FT) kesi(Num_Gauss_P_FEM_3D),yita(Num_Gauss_P_FEM_3D),zeta(Num_Gauss_P_FEM_3D),weight(Num_Gauss_P_FEM_3D)
real(kind=FT) globalF_T(c_Total_FD)
integer local(24)
real(kind=FT) J(3,3),dNdkesi(8,3),N(3,24),Inverse_J(3,3),dNdx(8,3)    
integer i_C,i_FluidEle,num_Flu_Nodes
real(kind=FT) Flu_Ele_Area,Flu_Ele_Vector(3),Flu_Ele_Centroid(3)
real(kind=FT) Weight_HF
integer Local_W(60),num_Local_W                           
real(kind=FT) N_W(3,60),ori_n(3)
real(kind=FT) HF_GAUSS_x,HF_GAUSS_y,HF_GAUSS_z
integer HF_GAUSS_Elem
real(kind=FT) c_S_xx,c_S_yy,c_S_zz,c_S_xy,c_S_yz,c_S_xz
real(kind=FT) Contact_Force_x,Contact_Force_y,Contact_Force_z
integer Num_Div_Points 
real(kind=FT) tem_n(8)
real(kind=FT) N_HF(100)  
real(kind=FT) c_F(3)
real(kind=FT) c_T_Alpha,c_TStress(6)
integer ref_elem
real(kind=FT) Tip_T_Matrx(3,3),r_Gauss,theta_Gauss,T(3,3)   
integer Solid_El
real(kind=FT)BaseLine_A(3),BaseLine_B(3)
real(kind=FT)BaseLine_x_Vec(3),BaseLine_y_Vec(3),BaseLine_z_Vec(3)      
real(kind=FT) BaseLine_Mid(3),Gauss_Coor_Local(3),Node_Coor_Local(3)    
real(kind=FT) c_Center(3)
integer Yes_Node_Dealed(num_Node)
integer num_t,num_h,num_j,cnt,tem,i_f
integer c_LocXFEM(3)
real(kind=FT) c_Beta(3),c_Alpha 
integer i_Node
real(kind=FT) c_Penalty_CS
integer c_Cr_Location
real(kind=FT) c_q_vector(3),c_Force_Vector(24) 
real(kind=FT),POINTER::globalF_thread(:,:)
integer i_Thread,c_thread,max_threads
integer i_SL,i_SL_El
character*200 temp_name
character*200 temp_name_low_case,temp_File_Surface_Load,low_temp_File_Surface_Load
LOGICAL alive
LOGICAL alive_low
integer c_Num_Sur_Elem
integer Tool_Count_Lines
logical Flag_Blank
real(kind=FT),ALLOCATABLE::Temp_DATA(:,:)
real(kind=FT) tem_result(24)
real(kind=FT) Tri_P1(3),Tri_P2(3),Tri_P3(3),Tri_Center(3),c_Np(3)
integer c_Solid_Element
real(kind=FT) c_Solid_Center(3),Vector_Ele_Center_Tri_Center(3)
real(kind=FT) c_angle,c_SL_El_Area
integer i_Node_4, Local_Node_Number(4)
integer c_Node_4
real(kind=FT) kesi_4(4),yita_4(4),zeta_4(4),weight_4(4)
real(kind=FT) L_vector(3)
real(kind=FT) c_Force_24x8(24,8),N8(8),c_Force_24(24)
real(kind=FT) p_vector(4)
real(kind=FT) c_Force_Vector_Total(12)
real(kind=FT) c_weight
real(kind=FT) N_3x12(3,12)
real(kind=FT) N1,N2,N3,N4
real(kind=FT) N_4(4),c_Force_12(12) 
real(kind=FT) c_Force_12x4(12,4)
real(kind=FT) c_Force_Vector_12(12)
real(kind=FT) SURFACE_X_NODES(4),SURFACE_Y_NODES(4),SURFACE_Z_NODES(4)
real(kind=FT) Mapped_SURFACE_X_NODES(4),Mapped_SURFACE_Y_NODES(4)
integer local12(12)
integer Position4(4)
integer Face_ID
real(kind=FT) c_SL_globalF(c_Total_FD)

real(kind=FT) c_Stress_6(6)

print *,'    Constructing global force vector...'        

globalF(1:c_Total_FD) = ZR
do i = 1,Num_Foc_x
  cur_Node = int(Foc_x(i,1))
  globalF(3*cur_Node-2) = Lambda*Foc_x(i,2)                   
end do

do i = 1,Num_Foc_y
  cur_Node = int(Foc_y(i,1))
  globalF(3*cur_Node-1) =   Lambda*Foc_y(i,2)                  
end do

do i = 1,Num_Foc_z
  cur_Node = int(Foc_z(i,1))
  globalF(3*cur_Node) =   Lambda*Foc_z(i,2)                  
end do     


if(Key_Gravity==1) then
  call Cal_Gauss_Points_3D_8nodes(Num_Gauss_P_FEM_3D,kesi,yita,zeta,weight)  
  do i_E = 1,Num_Elem
      mat_num    = Elem_Mat(i_E)
      c_mat_type = Material_Type(mat_num)
      if (c_mat_type ==1) then
          c_density = Material_Para(mat_num,3)
      else
          c_density = Material_Para(mat_num,3); 
      end if
      
      c_q_vector = -c_density*g_X_Y_Z(1:3)
      
      c_NN    = G_NN(1:8,i_E)
      c_X_NODES = G_X_NODES(1:8,i_E)
      c_Y_NODES = G_Y_NODES(1:8,i_E)    
      c_Z_NODES = G_Z_NODES(1:8,i_E)  
      local=[c_NN(1)*3-2,c_NN(1)*3-1,c_NN(1)*3,c_NN(2)*3-2,c_NN(2)*3-1,c_NN(2)*3,   &
             c_NN(3)*3-2,c_NN(3)*3-1,c_NN(3)*3,c_NN(4)*3-2,c_NN(4)*3-1,c_NN(4)*3,   &
             c_NN(5)*3-2,c_NN(5)*3-1,c_NN(5)*3,c_NN(6)*3-2,c_NN(6)*3-1,c_NN(6)*3,   &
             c_NN(7)*3-2,c_NN(7)*3-1,c_NN(7)*3,c_NN(8)*3-2,c_NN(8)*3-1,c_NN(8)*3]
             
      do i_G = 1,Num_Gauss_P_FEM_3D 
          c_kesi = kesi(i_G)
          c_yita = yita(i_G)
          c_zeta = zeta(i_G)
          call Cal_N_detJ_3D(c_kesi,c_yita,c_zeta,c_X_NODES,c_Y_NODES,c_Z_NODES,detJ,N)
          c_Force_Vector(1:24) = MATMUL(TRANSPOSE(N),c_q_vector) 
          globalF(local) = globalF(local) + c_Force_Vector(1:24)*weight(i_G)*detJ     
      end do
  end do
end if



if(Key_Thermal_Stress==1)then
    call Cal_Gauss_Points_3D_8nodes(Num_Gauss_P_FEM_3D,kesi,yita,zeta,weight)          
    do i_E = 1,Num_Elem
      if(Key_XA/=2) then
          c_D     = D(Elem_Mat(i_E),1:6,1:6)  
          c_T_Alpha = T_Alpha(Elem_Mat(i_E))
      else
          c_D       = Elem_D_XA(i_E,1:6,1:6)
          c_T_Alpha = Elem_TEC_XA(i_E)
      endif

      c_NN    = G_NN(1:8,i_E)
      c_X_NODES = G_X_NODES(1:8,i_E)
      c_Y_NODES = G_Y_NODES(1:8,i_E)    
      c_Z_NODES = G_Z_NODES(1:8,i_E)  
      local=[c_NN(1)*3-2,c_NN(1)*3-1,c_NN(1)*3,c_NN(2)*3-2,c_NN(2)*3-1,c_NN(2)*3,   &
             c_NN(3)*3-2,c_NN(3)*3-1,c_NN(3)*3,c_NN(4)*3-2,c_NN(4)*3-1,c_NN(4)*3,   &
             c_NN(5)*3-2,c_NN(5)*3-1,c_NN(5)*3,c_NN(6)*3-2,c_NN(6)*3-1,c_NN(6)*3,   &
             c_NN(7)*3-2,c_NN(7)*3-1,c_NN(7)*3,c_NN(8)*3-2,c_NN(8)*3-1,c_NN(8)*3]
      do i_G = 1,Num_Gauss_P_FEM_3D 
          c_kesi = kesi(i_G)
          c_yita = yita(i_G)
          c_zeta = zeta(i_G)
          call Cal_N_dNdkesi_J_detJ_3D(c_kesi,c_yita,c_zeta,c_X_NODES,c_Y_NODES,c_Z_NODES,detJ,J,N,dNdkesi)    
          call Matrix_Inverse_3x3(J,Inverse_J)
          dNdx = MATMUL(dNdkesi,Inverse_J)
          c_B(1:6,1:24) = ZR
          c_B(1,1:24:3)   =  dNdx(1:8,1)
          c_B(2,2:24:3)   =  dNdx(1:8,2)
          c_B(3,3:24:3)   =  dNdx(1:8,3)
          c_B(4,1:24:3)   =  dNdx(1:8,2)
          c_B(4,2:24:3)   =  dNdx(1:8,1)
          c_B(5,2:24:3)   =  dNdx(1:8,3)
          c_B(5,3:24:3)   =  dNdx(1:8,2)
          c_B(6,1:24:3)   =  dNdx(1:8,3)
          c_B(6,3:24:3)   =  dNdx(1:8,1)    
          
          
          c_TStress = c_T_Alpha*Elem_T_for_Stress(i_E)*MATMUL(c_D,[ONE,ONE,ONE,ZR,ZR,ZR])
          
          c_temp(1:24) =  MATMUL(transpose(c_B(1:6,1:24)),c_TStress)
          globalF(local) = globalF(local) +c_temp(1:24)*weight(i_G)*detJ     
      end do
    enddo
endif
      
      
if(Key_InSitu_Strategy==2)then
    if(abs(InSitu_x) < 0.0001D6 .and. abs(InSitu_y) < 0.0001D6 .and. abs(InSitu_z) < 0.0001D6  )then
       print *,'    #--#--#--#--#--#--#--#--#--#--#--#--#--#--##--#--#--#--#--#--#--#--#--#--#--#--#--#--#'
       print *,'    Warning :: cannot find InSitu_x, InSitu_y, InSitu_z (warning in Force_Vector_3D.f)!'
       print *,'    #--#--#--#--#--#--#--#--#--#--#--#--#--#--##--#--#--#--#--#--#--#--#--#--#--#--#--#--#'
       goto 199
    endif
      
    call Cal_Gauss_Points_3D_8nodes(Num_Gauss_P_FEM_3D,kesi,yita,zeta,weight)          
    G_Counter = 0
     
    do i_E = 1,Num_Elem
      c_D     = D(Elem_Mat(i_E),1:6,1:6)     
      c_NN    = G_NN(1:8,i_E)
      c_X_NODES = G_X_NODES(1:8,i_E)
      c_Y_NODES = G_Y_NODES(1:8,i_E)    
      c_Z_NODES = G_Z_NODES(1:8,i_E)        
      local=[c_NN(1)*3-2,c_NN(1)*3-1,c_NN(1)*3,c_NN(2)*3-2,c_NN(2)*3-1,c_NN(2)*3, &
             c_NN(3)*3-2,c_NN(3)*3-1,c_NN(3)*3,c_NN(4)*3-2,c_NN(4)*3-1,c_NN(4)*3, &
             c_NN(5)*3-2,c_NN(5)*3-1,c_NN(5)*3,c_NN(6)*3-2,c_NN(6)*3-1,c_NN(6)*3, &
             c_NN(7)*3-2,c_NN(7)*3-1,c_NN(7)*3,c_NN(8)*3-2,c_NN(8)*3-1,c_NN(8)*3]
      do i_G = 1,Num_Gauss_P_FEM_3D 
          G_Counter = G_Counter +1
          c_kesi = kesi(i_G)
          c_yita = yita(i_G)
          c_zeta = zeta(i_G)
          call Cal_N_dNdkesi_J_detJ_3D(c_kesi,c_yita,c_zeta,c_X_NODES,c_Y_NODES,c_Z_NODES,detJ,J,N,dNdkesi)    
          call Matrix_Inverse_3x3(J,Inverse_J)
          dNdx = MATMUL(dNdkesi,Inverse_J)
          c_B(1:6,1:24)   =  ZR
          c_B(1,1:24:3)   =  dNdx(1:8,1)
          c_B(2,2:24:3)   =  dNdx(1:8,2)
          c_B(3,3:24:3)   =  dNdx(1:8,3)
          c_B(4,1:24:3)   =  dNdx(1:8,2)
          c_B(4,2:24:3)   =  dNdx(1:8,1)
          c_B(5,2:24:3)   =  dNdx(1:8,3)
          c_B(5,3:24:3)   =  dNdx(1:8,2)
          c_B(6,1:24:3)   =  dNdx(1:8,3)
          c_B(6,3:24:3)   =  dNdx(1:8,1)    
          c_Stress = [InSitu_Strs_Gaus_xx(i_E,i_G),InSitu_Strs_Gaus_yy(i_E,i_G),InSitu_Strs_Gaus_zz(i_E,i_G), &
                      InSitu_Strs_Gaus_xy(i_E,i_G),InSitu_Strs_Gaus_yz(i_E,i_G),InSitu_Strs_Gaus_xz(i_E,i_G)]
          c_temp(1:24) =  MATMUL(transpose(c_B(1:6,1:24)),c_Stress)
          globalF(local) = globalF(local) - c_temp(1:24)*weight(i_G)*detJ     
      end do
    enddo
     
    if(Key_Analysis_Type==1 .and. num_Crack/=0)then    
      do i_C=1,num_Crack
        if(Key_Crack_Inner_Pressure==1  .and. (Crack_Pressure(i_C)>Tol_10)) then
          CONTINUE
        else
          cycle
        endif

        do i_FluidEle = 1,Cracks_FluidEle_num_3D(i_C)
          num_Flu_Nodes=Cracks_FluidEle_num_CalP_3D(i_C)%row(i_FluidEle)
          Flu_Ele_Area = Cracks_FluidEle_Area_3D(i_C)%row(i_FluidEle)  
          detJ  = Flu_Ele_Area 
          Flu_Ele_Vector(1:3) = Cracks_FluidEle_Vector_3D(i_C)%row(i_FluidEle,1:3) 
          ori_n = Flu_Ele_Vector(1:3)
          Flu_Ele_Centroid(1:3) = Cracks_FluidEle_Centroid_3D(i_C)%row(i_FluidEle,1:3)   
          c_Center = Cracks_FluidEle_Centroid_3D(i_C)%row(i_FluidEle,1:3)     
          HF_GAUSS_Elem = Cracks_FluidEle_EleNum_3D(i_C)%row(i_FluidEle)
          call Cal_KesiYita_by_Coor_3D(Flu_Ele_Centroid,HF_GAUSS_Elem,c_kesi,c_yita,c_zeta) 
 

          T=Cracks_FluidEle_LCS_T_3D(i_C)%row(i_FluidEle,1:3,1:3)
          Solid_El = Cracks_FluidEle_EleNum_3D(i_C)%row(i_FluidEle)        
        
          c_S_xx=sum(InSitu_Strs_Gaus_xx(HF_GAUSS_Elem,1:Num_Gauss_P_FEM_3D))/Num_Gauss_P_FEM_3D
          c_S_yy=sum(InSitu_Strs_Gaus_yy(HF_GAUSS_Elem,1:Num_Gauss_P_FEM_3D))/Num_Gauss_P_FEM_3D
          c_S_zz=sum(InSitu_Strs_Gaus_zz(HF_GAUSS_Elem,1:Num_Gauss_P_FEM_3D))/Num_Gauss_P_FEM_3D
          c_S_xy=sum(InSitu_Strs_Gaus_xy(HF_GAUSS_Elem,1:Num_Gauss_P_FEM_3D))/Num_Gauss_P_FEM_3D         
          c_S_yz=sum(InSitu_Strs_Gaus_yz(HF_GAUSS_Elem,1:Num_Gauss_P_FEM_3D))/Num_Gauss_P_FEM_3D
          c_S_xz=sum(InSitu_Strs_Gaus_xz(HF_GAUSS_Elem,1:Num_Gauss_P_FEM_3D))/Num_Gauss_P_FEM_3D   
          Contact_Force_x = (c_S_xx*ori_n(1)+c_S_xy*ori_n(2)+c_S_xz*ori_n(3))
          Contact_Force_y = (c_S_xy*ori_n(1)+c_S_yy*ori_n(2)+c_S_yz*ori_n(3))
          Contact_Force_z = (c_S_xz*ori_n(1)+c_S_yz*ori_n(2)+c_S_zz*ori_n(3))         
          c_NN(1:8)    = G_NN(1:8,HF_GAUSS_Elem)
          call Cal_N_3D(c_kesi,c_yita,c_zeta,N)     
          tem_n(1) = N(1,1);    tem_n(2) = N(1,4)
          tem_n(3) = N(1,7);    tem_n(4) = N(1,10)      
          tem_n(5) = N(1,13);   tem_n(6) = N(1,16)
          tem_n(7) = N(1,19);   tem_n(8) = N(1,22)    
          N_HF(1:num_Flu_Nodes) = ONE/num_Flu_Nodes   
          num_Local_W = 0
          N_W(1:3,1:60) =ZR
          do i_N = 1,8
            if (Enriched_Node_Type_3D(c_NN(i_N),i_C) ==2) then
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

              call Vector_Location_Int_v2(Solid_El_Max_num_Crs,Solid_El_Crs(ref_elem,1:Solid_El_Max_num_Crs),i_C,c_Cr_Location)
              
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

              N_W(1,3*num_Local_W-2)=TWO*sqrt(r_Gauss)*tem_n(i_N)
              N_W(2,3*num_Local_W-1)=TWO*sqrt(r_Gauss)*tem_n(i_N)
              N_W(3,3*num_Local_W)  =TWO*sqrt(r_Gauss)*tem_n(i_N)   
            endif
          end do
          c_F =[Contact_Force_x,Contact_Force_y,Contact_Force_z]
          
          globalF(Local_W(1:3*num_Local_W)) = globalF(Local_W(1:3*num_Local_W))+&
                   MATMUL(transpose(N_W(1:3,1:3*num_Local_W)),c_F)*Weight_HF*detJ      
          
        enddo
      enddo
    endif 
endif
      
      
if(Key_InSitu_Strategy==4)then
    max_threads = omp_get_max_threads()
    ALLOCATE(globalF_thread(1:c_Total_FD,1:max_threads))
    call Cal_Gauss_Points_3D_8nodes(Num_Gauss_P_FEM_3D,kesi,yita,zeta,weight)          
    
    
    globalF_thread = ZR   
    !$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(c_thread,i_E,i_G,c_D,           &
    !$OMP          c_NN,c_X_NODES,c_Y_NODES,c_Z_NODES,local,                &
    !$OMP          c_kesi,c_yita,c_zeta,detJ,J,N,dNdkesi,Inverse_J,         &
    !$OMP          dNdx,c_B,c_Stress,c_temp) 
    c_thread = omp_get_thread_num()+1
    !$OMP DO SCHEDULE(static)  
    do i_E = 1,Num_Elem
      if(Key_XA/=2) then
          c_D     = D(Elem_Mat(i_E),1:6,1:6)     
      else
          c_D     = Elem_D_XA(i_E,1:6,1:6)
      endif
      c_NN    = G_NN(1:8,i_E)
      c_X_NODES = G_X_NODES(1:8,i_E)
      c_Y_NODES = G_Y_NODES(1:8,i_E)    
      c_Z_NODES = G_Z_NODES(1:8,i_E)        
      local=[c_NN(1)*3-2,c_NN(1)*3-1,c_NN(1)*3,c_NN(2)*3-2,c_NN(2)*3-1,c_NN(2)*3,   &
             c_NN(3)*3-2,c_NN(3)*3-1,c_NN(3)*3,c_NN(4)*3-2,c_NN(4)*3-1,c_NN(4)*3,   &
             c_NN(5)*3-2,c_NN(5)*3-1,c_NN(5)*3,c_NN(6)*3-2,c_NN(6)*3-1,c_NN(6)*3,   &
             c_NN(7)*3-2,c_NN(7)*3-1,c_NN(7)*3,c_NN(8)*3-2,c_NN(8)*3-1,c_NN(8)*3]
      do i_G = 1,Num_Gauss_P_FEM_3D 
          c_kesi = kesi(i_G)
          c_yita = yita(i_G)
          c_zeta = zeta(i_G)
          call Cal_N_dNdkesi_J_detJ_3D(c_kesi,c_yita,c_zeta,c_X_NODES,c_Y_NODES,c_Z_NODES,detJ,J,N,dNdkesi)    
          call Matrix_Inverse_3x3(J,Inverse_J)
          dNdx = MATMUL(dNdkesi,Inverse_J)
          c_B(1:6,1:24)   =  ZR
          c_B(1,1:24:3)   =  dNdx(1:8,1)
          c_B(2,2:24:3)   =  dNdx(1:8,2)
          c_B(3,3:24:3)   =  dNdx(1:8,3)
          c_B(4,1:24:3)   =  dNdx(1:8,2)
          c_B(4,2:24:3)   =  dNdx(1:8,1)
          c_B(5,2:24:3)   =  dNdx(1:8,3)
          c_B(5,3:24:3)   =  dNdx(1:8,2)
          c_B(6,1:24:3)   =  dNdx(1:8,3)
          c_B(6,3:24:3)   =  dNdx(1:8,1)   
          
          c_Stress = MATMUL(c_D,[ InSitu_Strain_Gaus_xx(i_E,i_G),InSitu_Strain_Gaus_yy(i_E,i_G),InSitu_Strain_Gaus_zz(i_E,i_G),    &
                                  InSitu_Strain_Gaus_xy(i_E,i_G),InSitu_Strain_Gaus_yz(i_E,i_G),InSitu_Strain_Gaus_xz(i_E,i_G)]) 
                                  
          c_temp(1:24) =  MATMUL(transpose(c_B(1:6,1:24)),c_Stress)
          globalF_thread(local,c_thread) = globalF_thread(local,c_thread) + c_temp(1:24)*weight(i_G)*detJ         
      end do
      
    enddo
    !$omp end do 
    !$omp end parallel       
    DO i_Thread = 1,omp_get_max_threads()
        globalF  =  globalF  + globalF_thread(:,i_Thread)
    ENDDO     
    if(ASSOCIATED(globalF_thread)) DEALLOCATE(globalF_thread)
endif
 
if(Key_InStress_for_Mat==1)then
  call Cal_Gauss_Points_3D_8nodes(Num_Gauss_P_FEM_3D,kesi,yita,zeta,weight)   
  do i_E = 1,Num_Elem
      c_D     = D(Elem_Mat(i_E),1:6,1:6)  
      if (any(Mat_Number_of_InStress == Elem_Mat(i_E)))then   
          c_NN    = G_NN(1:8,i_E)
          c_X_NODES = G_X_NODES(1:8,i_E)
          c_Y_NODES = G_Y_NODES(1:8,i_E)    
          c_Z_NODES = G_Z_NODES(1:8,i_E)    
          
          local=[c_NN(1)*3-2,c_NN(1)*3-1,c_NN(1)*3,c_NN(2)*3-2,c_NN(2)*3-1,c_NN(2)*3,   &
                 c_NN(3)*3-2,c_NN(3)*3-1,c_NN(3)*3,c_NN(4)*3-2,c_NN(4)*3-1,c_NN(4)*3,   &
                 c_NN(5)*3-2,c_NN(5)*3-1,c_NN(5)*3,c_NN(6)*3-2,c_NN(6)*3-1,c_NN(6)*3,   &
                 c_NN(7)*3-2,c_NN(7)*3-1,c_NN(7)*3,c_NN(8)*3-2,c_NN(8)*3-1,c_NN(8)*3]
          do i_G = 1,Num_Gauss_P_FEM_3D  
              c_kesi = kesi(i_G)
              c_yita = yita(i_G)
              c_zeta = zeta(i_G)
              call Cal_N_dNdkesi_J_detJ_3D(c_kesi,c_yita,c_zeta,c_X_NODES,c_Y_NODES,c_Z_NODES,detJ,J,N,dNdkesi)    
              call Matrix_Inverse_3x3(J,Inverse_J)
              dNdx = MATMUL(dNdkesi,Inverse_J)
              c_B(1:6,1:24)   =  ZR
              c_B(1,1:24:3)   =  dNdx(1:8,1)
              c_B(2,2:24:3)   =  dNdx(1:8,2)
              c_B(3,3:24:3)   =  dNdx(1:8,3)
              c_B(4,1:24:3)   =  dNdx(1:8,2)
              c_B(4,2:24:3)   =  dNdx(1:8,1)
              c_B(5,2:24:3)   =  dNdx(1:8,3)
              c_B(5,3:24:3)   =  dNdx(1:8,2)
              c_B(6,1:24:3)   =  dNdx(1:8,3)
              c_B(6,3:24:3)   =  dNdx(1:8,1)   
              
              c_Stress_6 = [Mat_InStress_x,Mat_InStress_y,Mat_InStress_z,ZR,ZR,ZR]
                                      
              c_temp(1:24) =  MATMUL(transpose(c_B(1:6,1:24)),c_Stress_6)
              
              globalF(local) = globalF(local) + c_temp(1:24)*weight(i_G)*detJ 
          end do   
      end if
  enddo 
endif   
   
199 continue 
  
if (Num_Surface_Loads>=1) then
#ifndef Silverfrost
    do i_SL = 1,Num_Surface_Loads
        temp_name = trim(trim(Full_Pathname))//'.'//trim(trim(File_Surface_Load(i_SL)))
        inquire(file=temp_name, exist=alive)
        
        temp_File_Surface_Load = trim(trim(File_Surface_Load(i_SL)))
        call Tool_chrpak_s_low(temp_File_Surface_Load)
        temp_name_low_case = trim(trim(Full_Pathname))//'.'//trim(trim(temp_File_Surface_Load))
        inquire(file=temp_name_low_case, exist=alive_low)
        
        if((alive.EQV..FALSE.) .and. (alive_low.EQV..FALSE.))then
            print *, "    Error-2023012101 :: Can not find File_Surface_Load files!" 
            print *, "                        In Force_Vector_3D.f90"      
            print *, "    Missing file:",temp_name
            call Warning_Message('S',Keywords_Blank) 
        endif
    enddo
    
    if(.not. allocated(Surface_Load_Elements_Nodes)) then
        allocate(Surface_Load_Elements_Nodes(Num_Surface_Loads))
        allocate(Surface_Load_Elements_Normal(Num_Surface_Loads))
        allocate(Surface_Load_Elements_Area(Num_Surface_Loads))
    endif
    
    do i_SL = 1,Num_Surface_Loads
        temp_name = trim(trim(Full_Pathname))//'.'//trim(trim(File_Surface_Load(i_SL)))
        inquire(file=temp_name, exist=alive)
        
        temp_File_Surface_Load = trim(trim(File_Surface_Load(i_SL)))
        call Tool_chrpak_s_low(temp_File_Surface_Load)
        temp_name_low_case = trim(trim(Full_Pathname))//'.'//trim(trim(temp_File_Surface_Load))
        inquire(file=temp_name_low_case, exist=alive_low)
        
        c_SL_globalF(1:c_Total_FD) = ZR
        
        if(alive)then
            c_Num_Sur_Elem = Tool_Count_Lines(temp_name)
        elseif(alive_low)then
            c_Num_Sur_Elem = Tool_Count_Lines(temp_name_low_case)
        endif
        
        if(.not. allocated(Surface_Load_Elements_Nodes(i_SL)%row)) then
            allocate(Surface_Load_Elements_Nodes(i_SL)%row(c_Num_Sur_Elem,5))
            allocate(Surface_Load_Elements_Normal(i_SL)%row(c_Num_Sur_Elem,3))
            allocate(Surface_Load_Elements_Area(i_SL)%row(c_Num_Sur_Elem))
            
            ALLOCATE(Temp_DATA(c_Num_Sur_Elem,6))
            
            if(alive)then
                Call Tool_Read_File(temp_name,'blnk',c_Num_Sur_Elem,6, &
                                Temp_DATA(1:c_Num_Sur_Elem,1:6),Flag_Blank)
            elseif(alive_low)then
                Call Tool_Read_File(temp_name_low_case,'blnk',c_Num_Sur_Elem,6, &
                                Temp_DATA(1:c_Num_Sur_Elem,1:6),Flag_Blank)
            endif
        
            Surface_Load_Elements_Nodes(i_SL)%row(1:c_Num_Sur_Elem,1:5) = int(Temp_DATA(1:c_Num_Sur_Elem,1:5)) 
            Surface_Load_Elements_Area(i_SL)%row(1:c_Num_Sur_Elem) = Temp_DATA(1:c_Num_Sur_Elem,6)
            deallocate(Temp_DATA)
            
            do i_SL_El=1,c_Num_Sur_Elem
                c_Solid_Element = Surface_Load_Elements_Nodes(i_SL)%row(i_SL_El,1)
                c_Solid_Center = Elem_Centroid(c_Solid_Element,1:3)
                Tri_P1 = Coor(Surface_Load_Elements_Nodes(i_SL)%row(i_SL_El,2),1:3)
                Tri_P2 = Coor(Surface_Load_Elements_Nodes(i_SL)%row(i_SL_El,3),1:3)
                Tri_P3 = Coor(Surface_Load_Elements_Nodes(i_SL)%row(i_SL_El,4),1:3)
                call Tool_Normal_vector_of_3D_Tri(Tri_P1,Tri_P2,Tri_P3,c_Np)
                call Vector_Normalize(3,c_Np) 
                Tri_Center = (Tri_P1+Tri_P1+Tri_P1)/THR
                Vector_Ele_Center_Tri_Center = Tri_Center - c_Solid_Center
                call Vector_Normalize(3,Vector_Ele_Center_Tri_Center) 
                call Tool_Angle_of_Vectors_a_and_b_3D(c_Np,Vector_Ele_Center_Tri_Center,c_angle,3)
                if(c_angle< pi/TWO) then
                    Surface_Load_Elements_Normal(i_SL)%row(i_SL_El,1:3) = c_Np
                else
                    Surface_Load_Elements_Normal(i_SL)%row(i_SL_El,1:3) = -c_Np
                endif
            enddo
        endif
        

        do i_SL_El=1,c_Num_Sur_Elem
        
            c_Force_Vector_Total(1:12) = ZR
            
            c_Solid_Element = Surface_Load_Elements_Nodes(i_SL)%row(i_SL_El,1)
            
            
            c_NN    = G_NN(1:8,c_Solid_Element)
            
            c_X_NODES = G_X_NODES(1:8,c_Solid_Element)
            c_Y_NODES = G_Y_NODES(1:8,c_Solid_Element)    
            c_Z_NODES = G_Z_NODES(1:8,c_Solid_Element) 
            
            local=[c_NN(1)*3-2,c_NN(1)*3-1,c_NN(1)*3,c_NN(2)*3-2,c_NN(2)*3-1,c_NN(2)*3,   &
                   c_NN(3)*3-2,c_NN(3)*3-1,c_NN(3)*3,c_NN(4)*3-2,c_NN(4)*3-1,c_NN(4)*3,   &
                   c_NN(5)*3-2,c_NN(5)*3-1,c_NN(5)*3,c_NN(6)*3-2,c_NN(6)*3-1,c_NN(6)*3,   &
                   c_NN(7)*3-2,c_NN(7)*3-1,c_NN(7)*3,c_NN(8)*3-2,c_NN(8)*3-1,c_NN(8)*3]
                   
            c_SL_El_Area = Surface_Load_Elements_Area(i_SL)%row(i_SL_El)
            
            do i_Node_4=1,4
                c_Node_4 = Surface_Load_Elements_Nodes(i_SL)%row(i_SL_El,i_Node_4+1)
                call Vector_Location_Int_v2(8,c_NN,c_Node_4,Local_Node_Number(i_Node_4))   
            enddo
            
            if(sum(Local_Node_Number(1:4))==10) then
                Face_ID = 1
                Position4(1:4) = [4,3,2,1]
            elseif(sum(Local_Node_Number(1:4))==26) then
                Face_ID = 2
                Position4(1:4) = [5,6,7,8]
            elseif(sum(Local_Node_Number(1:4))==14) then
                Face_ID = 3
                Position4(1:4) = [1,2,6,5]
            elseif(sum(Local_Node_Number(1:4))==22) then
                Face_ID = 4
                Position4(1:4) = [3,4,8,7]   
            elseif(sum(Local_Node_Number(1:4))==18) then
                if(any(Local_Node_Number == 2)) then
                    Face_ID = 5
                    Position4(1:4) = [2,3,7,6]
                else
                    Face_ID = 6
                    Position4(1:4) = [4,1,5,8]
                endif
            endif

            call Cal_Gauss_Points_QUAD(4,kesi_4,yita_4,weight_4)
            
            L_vector = Surface_Load_Elements_Normal(i_SL)%row(i_SL_El,1:3)
            
            
            SURFACE_X_NODES(1) = c_X_NODES(Position4(1))
            SURFACE_X_NODES(2) = c_X_NODES(Position4(2))
            SURFACE_X_NODES(3) = c_X_NODES(Position4(3))
            SURFACE_X_NODES(4) = c_X_NODES(Position4(4))
            SURFACE_Y_NODES(1) = c_Y_NODES(Position4(1))
            SURFACE_Y_NODES(2) = c_Y_NODES(Position4(2))
            SURFACE_Y_NODES(3) = c_Y_NODES(Position4(3))
            SURFACE_Y_NODES(4) = c_Y_NODES(Position4(4))            
            SURFACE_Z_NODES(1) = c_Z_NODES(Position4(1))
            SURFACE_Z_NODES(2) = c_Z_NODES(Position4(2))
            SURFACE_Z_NODES(3) = c_Z_NODES(Position4(3))
            SURFACE_Z_NODES(4) = c_Z_NODES(Position4(4))  
            
            call Tool_Map_3D_Points_to_2D(4,SURFACE_X_NODES(1:4),SURFACE_Y_NODES(1:4),SURFACE_Z_NODES(1:4),&
                                          Mapped_SURFACE_X_NODES(1:4),Mapped_SURFACE_Y_NODES(1:4))

            
            do i_G = 1,4
                c_kesi   = kesi_4(i_G)
                c_yita   = yita_4(i_G)
                c_weight = weight_4(i_G)
                
                N1 = (ONE-c_kesi)*(ONE-c_yita)/FOU
                N2 = (ONE+c_kesi)*(ONE-c_yita)/FOU
                N3 = (ONE+c_kesi)*(ONE+c_yita)/FOU
                N4 = (ONE-c_kesi)*(ONE+c_yita)/FOU
                N_3x12(1,1:12) = [N1,ZR,ZR,N2,ZR,ZR,N3,ZR,ZR,N4,ZR,ZR]
                N_3x12(2,1:12) = [ZR,N1,ZR,ZR,N2,ZR,ZR,N3,ZR,ZR,N4,ZR]
                N_3x12(3,1:12) = [ZR,ZR,N1,ZR,ZR,N2,ZR,ZR,N3,ZR,ZR,N4]
                
                call Cal_detJ(c_kesi,c_yita,Mapped_SURFACE_X_NODES,Mapped_SURFACE_Y_NODES,detJ)
                
                N_4(1:4)= [N1,N2,N3,N4]
                
                c_Force_12(1:12) = MATMUL(TRANSPOSE(N_3x12),L_vector)
                
                call Vectors_Multi(c_Force_12(1:12),12,N_4(1:4),4,c_Force_12x4(1:12,1:4))   
                
                c_Force_12x4(1:12,1:4)= c_Force_12x4(1:12,1:4)*c_weight*detJ
                
                p_vector(1:4) =  -Surface_Pressure(i_SL)
                
                c_Force_Vector_12(1:12)  = MATMUL(c_Force_12x4(1:12,1:4),p_vector(1:4))
                
                c_Force_Vector_Total(1:12)  = c_Force_Vector_Total(1:12)  + c_Force_Vector_12(1:12)
                
                local12(1:12) = [c_NN(Position4(1))*3-2,c_NN(Position4(1))*3-1,c_NN(Position4(1))*3,&
                                 c_NN(Position4(2))*3-2,c_NN(Position4(2))*3-1,c_NN(Position4(2))*3,&
                                 c_NN(Position4(3))*3-2,c_NN(Position4(3))*3-1,c_NN(Position4(3))*3,&
                                 c_NN(Position4(4))*3-2,c_NN(Position4(4))*3-1,c_NN(Position4(4))*3]
                globalF(local12) = globalF(local12) + c_Force_Vector_12(1:12)
                
                c_SL_globalF(local12) = c_SL_globalF(local12) + c_Force_Vector_12(1:12)
            enddo
            
            
        enddo
        
        call Save_Surface_Load_3D(isub,i_SL,c_SL_globalF,c_Total_FD)
    
    enddo
#endif  
#ifdef Silverfrost
      print *,'    ERROR :: Silverfrost compiler failed to compile codes if (Num_Surface_Loads>=1)!'
      print *,'             In Force_Vector_3D.F90.'
      call Warning_Message('S',Keywords_Blank)
#endif  
endif


RETURN
END SUBROUTINE Force_Vector_3D
