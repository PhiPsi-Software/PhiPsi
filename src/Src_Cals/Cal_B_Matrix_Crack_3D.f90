 
subroutine Cal_B_Matrix_Crack_3D(kesi,yita,zeta,i_C,i_E,i_G, &
                     c_NN,c_X_NODES,c_Y_NODES,c_Z_NODES,tem_B,num_tem_B)  

use Global_Float_Type
use Global_Crack_3D
use Global_Model
use Global_Filename
use Global_Common
use Global_Material
use Global_Elem_Area_Vol
use Global_INTERFACE_D3_Get_Signed_Dis_to_Crack_Mesh
use Global_INTERFACE_D3_Get_Signed_Dis_to_Crack_Mesh_for_InPlane
use Global_Inter_Cal_N_dNdkesi_J_detJ_3D
use Global_Inter_Tool_Cal_Dis_Point_to_3D_Quad

implicit none

integer,intent(in)::i_C,i_E,i_G
real(kind=FT),intent(in)::c_X_NODES(8),c_Y_NODES(8),c_Z_NODES(8)
integer,intent(in)::c_NN(8)
real(kind=FT),intent(in)::kesi,yita,zeta
real(kind=FT),intent(out)::tem_B(6,MDOF_3D)
integer,intent(out)::num_tem_B        
real(kind=FT) detJ, J(3,3), Inverse_J(3,3)
real(kind=FT) N(3,24),dNdkesi(8,3),dNdx(8,3)
integer mat_num,c_mat_type
real(kind=FT) B_FEM(6,24),B_XFEM(6,MDOF_3D),BI_enr(6,3)
integer num_B_XFEM
integer i_N,ref_elem
real(kind=FT) Global_coor_Gauss(3)
real(kind=FT) H_i,H_gp
real(kind=FT) c_Min_Signed_Dis    
logical c_Yes_Node_PER_in_FS
real(kind=FT) c_PER_Node_to_FS(3) 
logical Yes_Found_Min_Signed_Dis
real(kind=FT) c_Distance_Node_to_FS
real(kind=FT)c_Distance
real(kind=FT),ALLOCATABLE::BI_enr_Tip(:,:)
real(kind=FT) aa,bb,cc,B_enr(6,3)
integer i_F,num_BI_enr_Tip
real(kind=FT) c_N(8),r_Gauss,theta_Gauss,omega
real(kind=FT) dFdx(4),dFdy(4),dFdz(4),F_Gauss(4),F_Node(4)
real(kind=FT)BaseLine_A(3),BaseLine_B(3)
real(kind=FT)BaseLine_x_Vec(3),BaseLine_y_Vec(3),BaseLine_z_Vec(3)
real(kind=FT) c_ThetaX,c_ThetaY,c_ThetaZ,c_T_Matrx(3,3)
real(kind=FT) r_Node
real(kind=FT) Node_Point(3),theta_Node
real(kind=FT) BaseLine_Mid(3),Gauss_Coor_Local(3),Node_Coor_Local(3),x0(3)
integer neg_C,Jun_elem
real(kind=FT) distance_Node_M,distance_Gauss_M,distance_Node_sm,distance_Gauss_sm
real(kind=FT) distance_x0,H_gp_M,H_gp_sm,H_i_M,H_i_sm,H_x0, J_gp, J_i  
real(kind=FT) Check_Ball_R
logical Yes_Found_Min_Signed_Dis_GM,Yes_Found_Min_Signed_Dis_Gsm
logical Yes_Found_Min_Signed_Dis_x0
real(kind=FT) c_n_Vector(3)
integer c_maxval_Enriched_Node_Type,c_minval_Enriched_Node_Type
real(kind=FT) c_Signed_Dis_v2
integer c_Cr_Location
integer :: temp_vector(Solid_El_Max_num_Crs)

tem_B(1:6,1:MDOF_3D) = ZR

c_maxval_Enriched_Node_Type =maxval(Enriched_Node_Type_3D(c_NN,i_C))
c_minval_Enriched_Node_Type =minval(Enriched_Node_Type_3D(c_NN,i_C))

if(c_maxval_Enriched_Node_Type .eq. 0 .and. c_minval_Enriched_Node_Type .eq. 0) then
  if(i_C >=2) then
      num_tem_B = 0
      return
  end if
endif

call Cal_N_dNdkesi_J_detJ_3D(kesi,yita,zeta,c_X_NODES,c_Y_NODES,c_Z_NODES,detJ,J,N,dNdkesi)    
c_N(1:8) = N(1,1:24:3)

call Matrix_Inverse_3x3(J,Inverse_J)

dNdx = MATMUL(dNdkesi,Inverse_J)


B_FEM(1:6,1:24) = ZR
if (EleGaus_yes_FEM_asemd(i_E,i_G).eqv. .False.) then
  B_FEM(1,1:24:3)   =  dNdx(:,1)
  B_FEM(2,2:24:3)   =  dNdx(:,2)
  B_FEM(3,3:24:3)   =  dNdx(:,3)
  B_FEM(4,1:24:3)   =  dNdx(:,2)
  B_FEM(4,2:24:3)   =  dNdx(:,1)
  B_FEM(5,2:24:3)   =  dNdx(:,3)
  B_FEM(5,3:24:3)   =  dNdx(:,2)
  B_FEM(6,1:24:3)   =  dNdx(:,3)
  B_FEM(6,3:24:3)   =  dNdx(:,1)
  EleGaus_yes_FEM_asemd(i_E,i_G)=.True.
end if

if(c_maxval_Enriched_Node_Type .eq. 0 .and. c_minval_Enriched_Node_Type .eq. 0) then
  if (i_C.eq.1) then
      tem_B(1:6,1:24) = B_FEM
      num_tem_B = 24
      return
  end if
  return
endif


allocate(BI_enr_Tip(6,Num_F_Functions*3))

if(c_maxval_Enriched_Node_Type .gt. 0 .or. c_minval_Enriched_Node_Type .gt. 0)then 
  Global_coor_Gauss(1) = DOT_PRODUCT(N(1,1:24:3),c_X_NODES(1:8))
  Global_coor_Gauss(2) = DOT_PRODUCT(N(1,1:24:3),c_Y_NODES(1:8))      
  Global_coor_Gauss(3) = DOT_PRODUCT(N(1,1:24:3),c_Z_NODES(1:8)) 

  
  if(Key_XA/=2) then
      mat_num    = Elem_Mat(i_E)
      c_mat_type = Material_Type(mat_num)  
  else
      c_mat_type = 1
  endif
      
  B_XFEM(1:6,1:MDOF_3D) = ZR
  num_B_XFEM = 0
  
  Yes_Found_Min_Signed_Dis = .False.
  c_Min_Signed_Dis  = Con_Big_20
  
  Check_Ball_R  = 3.0D0*Ave_Elem_L
  if(Key_InPlane_Growth== 0 ) then
      call D3_Get_Signed_Dis_to_Crack_Mesh(Global_coor_Gauss,i_C,&
                   Check_Ball_R,c_Distance,c_Signed_Dis_v2,      &
                    c_Yes_Node_PER_in_FS,c_PER_Node_to_FS(1:3),  &
                    Yes_Found_Min_Signed_Dis,c_n_Vector)
  elseif(Key_InPlane_Growth== 1) then
      call D3_Get_Signed_Dis_to_Crack_Mesh_for_InPlane_Growth(Global_coor_Gauss,i_C,&
                               Check_Ball_R,c_Distance,                      &
                               c_Yes_Node_PER_in_FS,c_PER_Node_to_FS(1:3),   &
                               Yes_Found_Min_Signed_Dis,c_n_Vector)
  endif
  
  do i_N = 1,8
      if (Enriched_Node_Type_3D(c_NN(i_N),i_C).eq.2)then
          if ((Elem_Type_3D(i_E,i_C).eq.2 ).or. (Elem_Type_3D(i_E,i_C).eq.3)) then
              if  (Elem_Type_3D(i_E,i_C) .eq. 3) then
       
              else
                 c_Distance_Node_to_FS = Dis_Node_to_FS(c_NN(i_N))%row(i_C)
                 
                 call Cal_Sign(c_Distance_Node_to_FS,H_i)
                

                 if(Yes_Found_Min_Signed_Dis .eqv. .True.) then
                     call Cal_Sign(c_Distance,H_gp)
                 else
                     H_gp = H_i
                 endif
              end if
                          
              
              BI_enr(1,1:3) = [dNdx(i_N,1)*(H_gp-H_i), ZR, ZR]
              BI_enr(2,1:3) = [ZR, dNdx(i_N,2)*(H_gp-H_i), ZR]
              BI_enr(3,1:3) = [ZR, ZR, dNdx(i_N,3)*(H_gp-H_i)]
              BI_enr(4,1:3) = [dNdx(i_N,2)*(H_gp-H_i), dNdx(i_N,1)*(H_gp-H_i),ZR]
              BI_enr(5,1:3) = [ZR,dNdx(i_N,3)*(H_gp-H_i), dNdx(i_N,2)*(H_gp-H_i)]
              BI_enr(6,1:3) = [dNdx(i_N,3)*(H_gp-H_i),ZR,dNdx(i_N,1)*(H_gp-H_i)]      
          else
              BI_enr(1:6,1:3)=ZR
          end if 
          
          B_XFEM(1:6,num_B_XFEM+1:num_B_XFEM+3) = BI_enr
          num_B_XFEM = num_B_XFEM + 3      
     
      elseif (Enriched_Node_Type_3D(c_NN(i_N),i_C).eq.3)then
        Jun_elem = 0 
        if(allocated(Node_Jun_elem_3D(c_NN(i_N))%row)) then
            Jun_elem = Node_Jun_elem_3D(c_NN(i_N))%row(i_C) 
        endif
        
        neg_C = 0
        if(allocated(Jun_Ele_Negative_Cr_Num_3D(Jun_elem)%row)) then
            neg_C = Jun_Ele_Negative_Cr_Num_3D(Jun_elem)%row(i_C) 
        endif
        
        x0 = ZR
        if(allocated(Coors_Junction_3D(Jun_elem)%row)) then
            x0 = Coors_Junction_3D(Jun_elem)%row(i_C,1:3) 
        endif        

        distance_Node_M = Dis_Node_to_FS(c_NN(i_N))%row(neg_C)
        
        Check_Ball_R  = 3.0D0*Ave_Elem_L  
        if(Key_InPlane_Growth== 0 ) then
            call D3_Get_Signed_Dis_to_Crack_Mesh(Global_coor_Gauss,neg_C,Check_Ball_R, &
                       c_Distance,c_Signed_Dis_v2,c_Yes_Node_PER_in_FS,            &
                       c_PER_Node_to_FS(1:3),Yes_Found_Min_Signed_Dis_GM,c_n_Vector)
        elseif(Key_InPlane_Growth== 1) then
            call D3_Get_Signed_Dis_to_Crack_Mesh_for_InPlane_Growth(Global_coor_Gauss,neg_C,&
                                       Check_Ball_R,c_Distance,                      &
                                       c_Yes_Node_PER_in_FS,c_PER_Node_to_FS(1:3),   &
                                       Yes_Found_Min_Signed_Dis,c_n_Vector)
        endif                       
        
        distance_Gauss_M = c_Distance
        distance_Node_sm = Dis_Node_to_FS(c_NN(i_N))%row(i_C)
        
        if(Key_InPlane_Growth== 0 ) then
            call D3_Get_Signed_Dis_to_Crack_Mesh(Global_coor_Gauss,i_C,Check_Ball_R,              &
                       c_Distance,c_Signed_Dis_v2,c_Yes_Node_PER_in_FS,c_PER_Node_to_FS(1:3), &
                       Yes_Found_Min_Signed_Dis_Gsm,c_n_Vector)
        elseif(Key_InPlane_Growth== 1) then
            call D3_Get_Signed_Dis_to_Crack_Mesh_for_InPlane_Growth(Global_coor_Gauss,i_C,&
                                       Check_Ball_R,c_Distance,                      &
                                       c_Yes_Node_PER_in_FS,c_PER_Node_to_FS(1:3),   &
                                       Yes_Found_Min_Signed_Dis,c_n_Vector)
        endif                         
        
        distance_Gauss_sm = c_Distance
        if(Key_InPlane_Growth== 0 ) then
            call D3_Get_Signed_Dis_to_Crack_Mesh(x0,neg_C,Check_Ball_R,c_Distance,c_Signed_Dis_v2, &
                       c_Yes_Node_PER_in_FS,c_PER_Node_to_FS(1:3),                             &
                       Yes_Found_Min_Signed_Dis_x0,c_n_Vector)    
        elseif(Key_InPlane_Growth== 1) then
            call D3_Get_Signed_Dis_to_Crack_Mesh_for_InPlane_Growth(x0,neg_C,        &
                                       Check_Ball_R,c_Distance,                      &
                                       c_Yes_Node_PER_in_FS,c_PER_Node_to_FS(1:3),   &
                                       Yes_Found_Min_Signed_Dis,c_n_Vector)
        endif                         
        
        distance_x0 = c_Distance         
        
        call Cal_Sign(distance_Gauss_M,H_gp_M)
        call Cal_Sign(distance_Node_M,H_i_M)
        call Cal_Sign(distance_Gauss_sm,H_gp_sm)
        call Cal_Sign(distance_Node_sm,H_i_sm)
        call Cal_Sign(distance_x0,H_x0)
        
        if (H_gp_M*H_x0 > ZR) then
            if (H_gp_sm >= ZR) then
                J_gp=  ONE
            else
                J_gp= -ONE
            end if
        else
            J_gp=  ZR
        end if
        if (H_i_M*H_x0 > ZR) then
            if (H_i_sm >= ZR)then
                J_i=  ONE
            else
                J_i= -ONE
            end if
        else
            J_i=  ZR
        end if  
          
        BI_enr(1,1:3) = [dNdx(i_N,1)*(J_gp-J_i), ZR, ZR]
        BI_enr(2,1:3) = [ZR, dNdx(i_N,2)*(J_gp-J_i), ZR]
        BI_enr(3,1:3) = [ZR, ZR, dNdx(i_N,3)*(J_gp-J_i)]
        BI_enr(4,1:3) = [dNdx(i_N,2)*(J_gp-J_i), dNdx(i_N,1)*(J_gp-J_i),ZR]
        BI_enr(5,1:3) = [ZR,dNdx(i_N,3)*(J_gp-J_i), dNdx(i_N,2)*(J_gp-J_i)]
        BI_enr(6,1:3) = [dNdx(i_N,3)*(J_gp-J_i),ZR,dNdx(i_N,1)*(J_gp-J_i)]
        B_XFEM(1:6,num_B_XFEM+1:num_B_XFEM+3) = BI_enr
        num_B_XFEM = num_B_XFEM + 3  
      elseif (Enriched_Node_Type_3D(c_NN(i_N),i_C).eq.1)then
          if ((Elem_Type_3D(i_E,i_C).eq.1 )) then
              ref_elem=i_E
          else
              ref_elem=Ele_Num_Tip_Enriched_Node_3D(c_NN(i_N))%row(i_C)
          endif
          
          temp_vector = Solid_El_Crs(ref_elem,1:Solid_El_Max_num_Crs)
          call Vector_Location_Int_v2(Solid_El_Max_num_Crs, temp_vector, i_C, c_Cr_Location) 
                     
          BaseLine_A = Solid_El_Tip_BaseLine(ref_elem)%row(c_Cr_Location,1,1:3)
          BaseLine_B = Solid_El_Tip_BaseLine(ref_elem)%row(c_Cr_Location,2,1:3)
          BaseLine_Mid = (BaseLine_A+BaseLine_B)/TWO
          BaseLine_x_Vec=Solid_El_Tip_BaseLine_x_Vec(ref_elem)%row(c_Cr_Location,1:3)
          BaseLine_y_Vec=Solid_El_Tip_BaseLine_y_Vec(ref_elem)%row(c_Cr_Location,1:3)
          BaseLine_z_Vec=Solid_El_Tip_BaseLine_z_Vec(ref_elem)%row(c_Cr_Location,1:3)
          c_T_Matrx(1:3,1:3)=Solid_El_Tip_BaseLine_T_Matrix(ref_elem)%row(c_Cr_Location,1:3,1:3)
          c_ThetaX=Solid_El_Tip_BaseLine_T_theta(ref_elem)%row(c_Cr_Location,1)
          c_ThetaY=Solid_El_Tip_BaseLine_T_theta(ref_elem)%row(c_Cr_Location,2)
          c_ThetaZ=Solid_El_Tip_BaseLine_T_theta(ref_elem)%row(c_Cr_Location,3)

          Gauss_Coor_Local= MATMUL(c_T_Matrx,Global_coor_Gauss-BaseLine_Mid)
          r_Gauss = sqrt(Gauss_Coor_Local(1)**2 + Gauss_Coor_Local(2)**2)
          theta_Gauss = atan2(Gauss_Coor_Local(2),Gauss_Coor_Local(1))    
          Node_Point = [c_X_NODES(i_N),c_Y_NODES(i_N),c_Z_NODES(i_N)]    
          Node_Coor_Local= MATMUL(c_T_Matrx,Node_Point-BaseLine_Mid)
          r_Node = sqrt(Node_Coor_Local(1)**2 + Node_Coor_Local(2)**2)
          theta_Node = atan2(Node_Coor_Local(2),Node_Coor_Local(1))    

          
          call Cal_F_dFdx_dFdy_dFdz_3D(r_Gauss,theta_Gauss,c_T_Matrx,c_mat_type,F_Gauss,dFdx,dFdy,dFdz)
                          
          call Cal_F(r_Node,theta_Node,omega,c_mat_type,F_Node)
     
          num_BI_enr_Tip = 0
          
          
          BI_enr_Tip= ZR
          
          do i_F =1,Num_F_Functions                       
              aa = dNdx(i_N,1)*(F_Gauss(i_F)-F_Node(i_F)) + c_N(i_N)*dFdx(i_F)
              bb = dNdx(i_N,2)*(F_Gauss(i_F)-F_Node(i_F)) + c_N(i_N)*dFdy(i_F)
              cc = dNdx(i_N,3)*(F_Gauss(i_F)-F_Node(i_F)) + c_N(i_N)*dFdz(i_F)    

              B_enr(1,1:3) = [aa,  ZR,  ZR]
              B_enr(2,1:3) = [ZR,  bb,  ZR]
              B_enr(3,1:3) = [ZR,  ZR,  cc]
              B_enr(4,1:3) = [bb,  aa,  ZR]
              B_enr(5,1:3) = [ZR,  cc,  bb]
              B_enr(6,1:3) = [cc,  ZR,  aa]
              BI_enr_Tip(1:6,num_BI_enr_Tip+1:num_BI_enr_Tip+3)  = B_enr  
              num_BI_enr_Tip = num_BI_enr_Tip +3
          end do
          
          B_XFEM(1:6,num_B_XFEM+1:num_B_XFEM+Num_F_Functions*3) = BI_enr_Tip
          num_B_XFEM = num_B_XFEM + Num_F_Functions*3
      end if
  end do
  
  if (i_C.eq.1) then
      tem_B(1:6,1:24)             = B_FEM  
      tem_B(1:6,25:24+num_B_XFEM) = B_XFEM(1:6,1:num_B_XFEM)
      num_tem_B = 24 + num_B_XFEM   
  else
      tem_B(1:6,1:num_B_XFEM) = B_XFEM(1:6,1:num_B_XFEM)
      num_tem_B = num_B_XFEM
  end if 
  
  
  
end if   

deallocate(BI_enr_Tip)

RETURN
END SUBROUTINE Cal_B_Matrix_Crack_3D