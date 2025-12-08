 
subroutine Cal_Any_Point_Disp_KesiYita_3D(in_Elem,c_Kesi,c_Yita,c_Zeta,c_DISP,P_Disp)

use Global_Float_Type
use Global_Common
use Global_Crack_Common
use Global_Crack_3D
use Global_Model
use Global_Elem_Area_Vol
use Global_DISP
use Global_Inter_Tool_Cal_Dis_Point_to_3D_Quad
use Global_INTERFACE_D3_Get_Signed_Dis_to_Crack_Mesh

implicit none

integer,intent(in)::in_Elem
real(kind=FT),intent(in)::c_Kesi,c_Yita,c_Zeta,c_DISP(Total_FD)
real(kind=FT),intent(out)::P_Disp(3)

real(kind=FT) c_X_NODES(8),c_Y_NODES(8),c_Z_NODES(8),  &
                detJ,dNdkesi(8,3),J(3,3),N(3,24),c_N(8),Gauss_P(3)
integer c_Elem,i_C,i_N,NODES_iE(8),Num_Enri_Node
real(kind=FT) Disp_x,Disp_y,Disp_z
real(kind=FT) H_i,H_gp
real(kind=FT) c_Distance_Node,c_Distance_Gauss

real(kind=FT) c_Min_Signed_Dis
integer i_Crack_Ele
integer Crack_Node1,Crack_Node2,Crack_Node3
real(kind=FT) Point1(3),Point2(3),Point3(3)
real(kind=FT) Tri123_Distance,Tri123_PER(3)
logical Tri123_Yes_PER_in,Tri123_Yes_PER_on  
logical Yes_Found_Min_Signed_Dis
real(kind=FT) c_Distance
integer i_F
real(kind=FT) r_Gauss,theta_Gauss,omega
real(kind=FT) F_Gauss(4),F_Node(4)
real(kind=FT)BaseLine_A(3),BaseLine_B(3)
real(kind=FT)BaseLine_x_Vec(3),BaseLine_y_Vec(3),BaseLine_z_Vec(3)
real(kind=FT) c_ThetaX,c_ThetaY,c_ThetaZ,c_T_Matrx(3,3)
real(kind=FT) r_Node
real(kind=FT) Node_Point(3),theta_Node,F
integer ref_elem,mat_num,c_mat_type
real(kind=FT) BaseLine_Mid(3),Gauss_Coor_Local(3),Node_Coor_Local(3)
integer neg_C,Jun_elem
real(kind=FT)    distance_Node_M,distance_Gauss_M,   &
                   distance_Node_sm,distance_Gauss_sm,   &
                   distance_x0,H_gp_M,H_gp_sm,H_i_M,H_i_sm,  &
                   H_x0, J_gp, J_i  
real(kind=FT) Check_Ball_R,x0(3)     
real(kind=FT) c_PER_Node_to_FS(3),c_n_Vector(3)
logical c_Yes_Node_PER_in_FS     
real(kind=FT) c_Signed_Dis_v2
integer c_Cr_Location
integer c_POS_3D_c_Ele(8)
integer::temp_Solid_El_Crs(Solid_El_Max_num_Crs)
  

c_Elem = in_Elem

if(Key_XA/=2) then
      mat_num    = Elem_Mat(c_Elem)
      c_mat_type = Material_Type(mat_num) 
else
      c_mat_type = 1
endif

NODES_iE = G_NN(1:8,c_Elem)

c_X_NODES = G_X_NODES(1:8,c_Elem)
c_Y_NODES = G_Y_NODES(1:8,c_Elem)
c_Z_NODES = G_Z_NODES(1:8,c_Elem)

call Cal_Coor_by_KesiYita_3D(c_Kesi,c_Yita,c_Zeta,c_X_NODES,c_Y_NODES,c_Z_NODES,Gauss_P(1),Gauss_P(2),Gauss_P(3))

call Cal_N_dNdkesi_J_detJ_3D(c_Kesi,c_Yita,c_Zeta,c_X_NODES,c_Y_NODES,c_Z_NODES,detJ,J,N,dNdkesi)
c_N(1:8) = N(1,1:24:3)

Disp_x = c_DISP(NODES_iE(1)*3-2)*c_N(1) +  c_DISP(NODES_iE(2)*3-2)*c_N(2) +  &
           c_DISP(NODES_iE(3)*3-2)*c_N(3) +  c_DISP(NODES_iE(4)*3-2)*c_N(4) + &
           c_DISP(NODES_iE(5)*3-2)*c_N(5) +  c_DISP(NODES_iE(6)*3-2)*c_N(6) +  &
           c_DISP(NODES_iE(7)*3-2)*c_N(7) +  c_DISP(NODES_iE(8)*3-2)*c_N(8)

Disp_y = c_DISP(NODES_iE(1)*3-1)*c_N(1) +  c_DISP(NODES_iE(2)*3-1)*c_N(2) +  &
           c_DISP(NODES_iE(3)*3-1)*c_N(3) +  c_DISP(NODES_iE(4)*3-1)*c_N(4) + &
           c_DISP(NODES_iE(5)*3-1)*c_N(5) +  c_DISP(NODES_iE(6)*3-1)*c_N(6) +  &
           c_DISP(NODES_iE(7)*3-1)*c_N(7) +  c_DISP(NODES_iE(8)*3-1)*c_N(8)

Disp_z = c_DISP(NODES_iE(1)*3)*c_N(1) + c_DISP(NODES_iE(2)*3)*c_N(2) +  &
           c_DISP(NODES_iE(3)*3)*c_N(3) + c_DISP(NODES_iE(4)*3)*c_N(4) + &
           c_DISP(NODES_iE(5)*3)*c_N(5) + c_DISP(NODES_iE(6)*3)*c_N(6) +  &
           c_DISP(NODES_iE(7)*3)*c_N(7) + c_DISP(NODES_iE(8)*3)*c_N(8)
  
if(Num_Crack /=0)then
    do i_C=1,Num_Crack
      
      c_POS_3D_c_Ele(1:8) = c_POS_3D(NODES_iE,i_C)
      if(sum(c_POS_3D_c_Ele(1:8))==0) cycle

      Yes_Found_Min_Signed_Dis = .False.
      c_Min_Signed_Dis = 1.0D20
      do i_Crack_Ele=1,Crack3D_Meshed_Ele_num(i_C)
          Crack_Node1 = Crack3D_Meshed_Ele(i_C)%row(i_Crack_Ele,1)
          Crack_Node2 = Crack3D_Meshed_Ele(i_C)%row(i_Crack_Ele,2)
          Crack_Node3 = Crack3D_Meshed_Ele(i_C)%row(i_Crack_Ele,3)
          Point1(1:3) =Crack3D_Meshed_Node(i_C)%row(Crack_Node1,1:3)
          Point2(1:3) =Crack3D_Meshed_Node(i_C)%row(Crack_Node2,1:3)
          Point3(1:3) =Crack3D_Meshed_Node(i_C)%row(Crack_Node3,1:3)
          call Tool_Dis_Point_to_3D_Tri       &
            (Gauss_P,Point1,Point2,Point3, &
              Tri123_Distance,Tri123_PER,Tri123_Yes_PER_in,Tri123_Yes_PER_on)
 
         if(Tri123_Yes_PER_in .or. Tri123_Yes_PER_on) then
           if(abs(Tri123_Distance)<c_Min_Signed_Dis) then
              c_Distance_Gauss = Tri123_Distance
              c_Min_Signed_Dis = abs(Tri123_Distance) 
              Yes_Found_Min_Signed_Dis = .True.
           endif
         endif
      enddo

           
      
      do i_N = 1,8
          if (Enriched_Node_Type_3D(NODES_iE(i_N),i_C) ==2) then
              Num_Enri_Node = c_POS_3D(NODES_iE(i_N),i_C)
              if ((Elem_Type_3D(c_Elem,i_C) == 2)) then
                  c_Distance_Node =Dis_Node_to_FS(NODES_iE(i_N))%row(i_C)
                  
                  call Cal_Sign(c_Distance_Node,H_i)
                  if(Yes_Found_Min_Signed_Dis .eqv. .True.) then
                     call Cal_Sign(c_Distance_Gauss,H_gp)
                  else
                     H_gp = H_i
                  endif   
                  Disp_x = Disp_x+(H_gp-H_i)*c_N(i_N)*c_DISP(Num_Enri_Node*3-2)
                  Disp_y = Disp_y+(H_gp-H_i)*c_N(i_N)*c_DISP(Num_Enri_Node*3-1)
                  Disp_z = Disp_z+(H_gp-H_i)*c_N(i_N)*c_DISP(Num_Enri_Node*3)
              end if
          elseif (Enriched_Node_Type_3D(NODES_iE(i_N),i_C).eq.3)then
            Num_Enri_Node = c_POS_3D(NODES_iE(i_N),i_C)
            
            Jun_elem = 0 
            if(allocated(Node_Jun_elem_3D(NODES_iE(i_N))%row)) then
                Jun_elem = Node_Jun_elem_3D(NODES_iE(i_N))%row(i_C) 
            endif
            
            neg_C = 0
            if(allocated(Jun_Ele_Negative_Cr_Num_3D(Jun_elem)%row)) then
                neg_C = Jun_Ele_Negative_Cr_Num_3D(Jun_elem)%row(i_C) 
            endif
            
            x0 = ZR
            if(allocated(Coors_Junction_3D(Jun_elem)%row)) then
                x0 = Coors_Junction_3D(Jun_elem)%row(i_C,1:3) 
            endif
            

            distance_Node_M = Dis_Node_to_FS(NODES_iE(i_N))%row(neg_C)
            
            Check_Ball_R  = 3.0D0*Node_Max_L(NODES_iE(i_N))
            call D3_Get_Signed_Dis_to_Crack_Mesh(Gauss_P,neg_C,Check_Ball_R,  &
                        c_Distance,c_Signed_Dis_v2,c_Yes_Node_PER_in_FS,  &
                        c_PER_Node_to_FS(1:3),Yes_Found_Min_Signed_Dis,c_n_Vector)
            distance_Gauss_M = c_Distance
            
            distance_Node_sm = Dis_Node_to_FS(NODES_iE(i_N))%row(i_C)
            
            call D3_Get_Signed_Dis_to_Crack_Mesh(Gauss_P,i_C,Check_Ball_R,  &
                        c_Distance,c_Signed_Dis_v2,c_Yes_Node_PER_in_FS,     &
                        c_PER_Node_to_FS(1:3),Yes_Found_Min_Signed_Dis,c_n_Vector)
            distance_Gauss_sm = c_Distance
            call D3_Get_Signed_Dis_to_Crack_Mesh(x0,neg_C,Check_Ball_R,  &
                        c_Distance,c_Signed_Dis_v2,c_Yes_Node_PER_in_FS,  &
                        c_PER_Node_to_FS(1:3),Yes_Found_Min_Signed_Dis,c_n_Vector)
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
            
            
            Disp_x = Disp_x+(J_gp-J_i)*c_N(i_N)*c_DISP(Num_Enri_Node*3-2)
            Disp_y = Disp_y+(J_gp-J_i)*c_N(i_N)*c_DISP(Num_Enri_Node*3-1)
            Disp_z = Disp_z+(J_gp-J_i)*c_N(i_N)*c_DISP(Num_Enri_Node*3)                

          elseif (Enriched_Node_Type_3D(NODES_iE(i_N),i_C).eq.1)then
              if ((Elem_Type_3D(in_Elem,i_C).eq.1 )) then
                  ref_elem=in_Elem
              else
                  ref_elem=Ele_Num_Tip_Enriched_Node_3D(NODES_iE(i_N))%row(i_C)
              endif
              
              
              temp_Solid_El_Crs = Solid_El_Crs(ref_elem, 1:Solid_El_Max_num_Crs)
              call Vector_Location_Int_v2(Solid_El_Max_num_Crs, temp_Solid_El_Crs, i_C, c_Cr_Location)  
              
              BaseLine_A = Solid_El_Tip_BaseLine(ref_elem)%row(c_Cr_Location,1,1:3)
              BaseLine_B = Solid_El_Tip_BaseLine(ref_elem)%row(c_Cr_Location,2,1:3)
              BaseLine_Mid = (BaseLine_A+BaseLine_B)/TWO
              c_T_Matrx(1:3,1:3)=Solid_El_Tip_BaseLine_T_Matrix(ref_elem)%row(c_Cr_Location,1:3,1:3)
         
              Gauss_Coor_Local= MATMUL(c_T_Matrx,Gauss_P-BaseLine_Mid)
              r_Gauss = sqrt(Gauss_Coor_Local(1)**2 + Gauss_Coor_Local(2)**2)
              theta_Gauss = atan2(Gauss_Coor_Local(2),Gauss_Coor_Local(1))    
 
              Node_Point = [c_X_NODES(i_N),c_Y_NODES(i_N),c_Z_NODES(i_N)]    
              Node_Coor_Local= MATMUL(c_T_Matrx,Node_Point-BaseLine_Mid)
              r_Node = sqrt(Node_Coor_Local(1)**2 + Node_Coor_Local(2)**2)
              theta_Node = atan2(Node_Coor_Local(2),Node_Coor_Local(1))    

              call Cal_F(r_Gauss,theta_Gauss,omega,c_mat_type,F_Gauss)
              
              call Cal_F(r_Node,theta_Node,omega,c_mat_type,F_Node)
 
              
              do i_F =1,Num_F_Functions
                  Num_Enri_Node=c_POS_3D(NODES_iE(i_N),i_C)+i_F-1 
                  F = F_Gauss(i_F)-F_Node(i_F)
                  Disp_x=Disp_x+F*c_N(i_N)*c_DISP(Num_Enri_Node*3-2)
                  Disp_y=Disp_y+F*c_N(i_N)*c_DISP(Num_Enri_Node*3-1)
                  Disp_z=Disp_z+F*c_N(i_N)*c_DISP(Num_Enri_Node*3)     
              end do                  
          end if
      end do
    end do
endif
    
P_Disp(1) = Disp_x
P_Disp(2) = Disp_y
P_Disp(3) = Disp_z
  
return 
end SUBROUTINE Cal_Any_Point_Disp_KesiYita_3D            
