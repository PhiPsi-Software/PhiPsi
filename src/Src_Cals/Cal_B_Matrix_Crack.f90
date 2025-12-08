 
subroutine Cal_B_Matrix_Crack(kesi,yita,i_C,i_E,i_G,c_NN,c_X_NODES,c_Y_NODES,tem_B,num_tem_B)  

use Global_Float_Type
use Global_Crack
use Global_Crack_Common
use Global_Model
use Global_Filename
use Global_Common
use Global_Material

implicit none

integer,intent(in)::i_C,i_E,i_G
real(kind=FT),intent(in)::c_X_NODES(4),c_Y_NODES(4)
integer,intent(in)::c_NN(4)
real(kind=FT),intent(in)::kesi,yita
real(kind=FT),intent(out)::tem_B(3,80)
integer,intent(out)::num_tem_B        
real(kind=FT) detJ, J(2,2), Inverse_J(2,2)
real(kind=FT) N(2,8),dNdkesi(4,2),dNdx(4,2),Coor_AB(2,2)
integer mat_num,c_mat_type
real(kind=FT) B_FEM(3,8),B_XFEM(3,60),BI_enr(3,2)
integer num_B_XFEM
integer i_N,ref_elem,i_F
logical Flag_Check_Kink_Tip,Yes_Gauss_in_BCD,Yes_in_KAB
real(kind=FT) distance_Vertex,Global_coor_Gauss(2)
real(kind=FT) Tri_KAB(3,2),Closed_Tri_KAB(4,2)
real(kind=FT) H_i,omega, Coor_Tip(2),Trans_Matrix(2,2)
real(kind=FT) H_Vertex,distance_Gauss,H_gp,distance_Node
integer Jun_elem,i_tem,num_points_crack
real(kind=FT) Coor_CD(2,2),x0(2),distance_Node_M,distance_Gauss_M,&
   distance_Node_sm,distance_Gauss_sm,distance_x0,&
    H_gp_M,H_gp_sm,H_i_M,H_i_sm,H_x0, J_gp, J_i,c_Segment(2)
integer find_element(8)
real(kind=FT) Coor_A(2),Coor_B(2),Coor_C(2),coor_D(2),len_AB,Line_AB(2,2),&
    new_Line_AB(2,2),tip_x,tip_y,Tri_BCD(3,2),Closed_Tri_BCD(4,2),&
    angle_AB_BC,Line_BC(2,2),Line_BGauss(2,2),angle_AB_BGauss, &
    len_BGauss,thete_Doc
real(kind=FT) Global_coor_Gauss_Doc(2),xp_Gauss(2),r_Gauss,&
            theta_Gauss,xp_Node(2),r_Node,theta_Node
real(kind=FT) dFdx(4),dFdy(4),F_Gauss(4),F_Node(4)
real(kind=FT) c_N(4),aa,bb,B_enr(3,2),BI_enr_Tip(3,8)
integer num_BI_enr_Tip
integer c_Hole
real(kind=FT) x0_Hole,y0_Hole,R0_Hole
real(kind=FT) a_Hole,b_Hole,theta_Hole
tem_B(1:3,1:80) = ZR

call Cal_N_dNdkesi_J_detJ(kesi,yita,c_X_NODES,c_Y_NODES,detJ,J,N,dNdkesi)    

call Matrix_Inverse_2x2(J,Inverse_J)

dNdx = MATMUL(dNdkesi,Inverse_J)

Global_coor_Gauss(1) = DOT_PRODUCT(N(1,1:7:2),c_X_NODES(1:4))
Global_coor_Gauss(2) = DOT_PRODUCT(N(1,1:7:2),c_Y_NODES(1:4))      


mat_num    = Elem_Mat(i_E)

c_mat_type = Material_Type(mat_num)     

B_FEM(1:3,1:8) = ZR

if (EleGaus_yes_FEM_asemd(i_E,i_G) .eqv. .False.) then
  B_FEM(1,1:8:2)   =  dNdx(:,1)
  B_FEM(2,2:8:2)   =  dNdx(:,2)
  B_FEM(3,1:8:2)   =  dNdx(:,2)
  B_FEM(3,2:8:2)   =  dNdx(:,1)
  EleGaus_yes_FEM_asemd(i_E,i_G)=.True.
end if

if(maxval(Enriched_Node_Type(c_NN,i_C)).eq.0 .and.minval(Enriched_Node_Type(c_NN,i_C)).eq.0) then
  if (i_C.eq.1) then

      tem_B(1:3,1:8) = B_FEM
      num_tem_B = 8
  else
      num_tem_B = 0
  end if
elseif(maxval(Enriched_Node_Type(c_NN,i_C)).gt.0 .or. &
     minval(Enriched_Node_Type(c_NN,i_C)).gt.0) then 
  B_XFEM(1:3,1:60) = ZR
  num_B_XFEM = 0
  do i_N = 1,4
      Flag_Check_Kink_Tip =  .False.
      Yes_Gauss_in_BCD    =  .False.
      if (Enriched_Node_Type(c_NN(i_N),i_C).eq.2)then
          ref_elem = 0
          if ((Elem_Type(i_E,i_C).eq.2 ) .or.  (Elem_Type(i_E,i_C).eq.3)) then
              ref_elem = i_E
              Coor_AB(1,1:2)=[Coors_Element_Crack(ref_elem,i_C,1),Coors_Element_Crack(ref_elem,i_C,2)]
              Coor_AB(2,1:2)=[Coors_Element_Crack(ref_elem,i_C,3),Coors_Element_Crack(ref_elem,i_C,4)]  
              if  (Elem_Type(i_E,i_C) .eq. 3) then
                  
                  
                  call Cal_Signed_Distance(Coor_AB,Coors_Vertex(ref_elem,:),distance_Vertex)
                  call Cal_Sign(distance_Vertex,H_Vertex)
                  call Cal_Signed_Distance(Coor_AB,Global_coor_Gauss,distance_Gauss)
                  call Cal_Sign(distance_Gauss,H_gp)
                  if(Key_Heaviside_Value==0)then
                      if (H_gp<=ZR) H_gp = ZR
                  endif
                  if (H_Vertex*H_gp <= 0) then
                      H_gp = H_gp
                  else
                      Tri_KAB(1,:) =Coor_AB(1,:)
                      Tri_KAB(2,:) =Coor_AB(2,:)
                      Tri_KAB(3,:) =Coors_Vertex(ref_elem,:)
                      Closed_Tri_KAB(1:3,:) =  Tri_KAB
                      Closed_Tri_KAB(4,:) =  Tri_KAB(1,:)
                      call Tool_Yes_In_Poly(Global_coor_Gauss(1),Global_coor_Gauss(2),&
                           Closed_Tri_KAB(:,1), Closed_Tri_KAB(:,2),4,Yes_in_KAB)

                      if (Yes_in_KAB .eqv. .True.) then
                          H_gp = -H_gp
                      end if
                  end if       
                  call Cal_Signed_Distance(Coor_AB,[c_X_NODES(i_N),c_Y_NODES(i_N)],distance_Node)
                  call Cal_Sign(distance_Node,H_i)  
                  if(Key_Heaviside_Value==0)then
                      if (H_i<=ZR) H_i = ZR
                  endif
              else
                  call Cal_Signed_Distance(Coor_AB,Global_coor_Gauss,distance_Gauss)
                  call Cal_Sign(distance_Gauss,H_gp)
                  if(Key_Heaviside_Value==0)then
                      if (H_gp<=ZR) H_gp = ZR
                  endif
                  call Cal_Signed_Distance(Coor_AB,[c_X_NODES(i_N),c_Y_NODES(i_N)],distance_Node)
                  call Cal_Sign(distance_Node,H_i)
                  if(Key_Heaviside_Value==0)then
                      if (H_i<=ZR) H_i = ZR
                  endif
              end if
              BI_enr(1,:) = [dNdx(i_N,1)*(H_gp-H_i), ZR]
              BI_enr(2,:) = [ZR,dNdx(i_N,2)*(H_gp-H_i)]
              BI_enr(3,:) = [dNdx(i_N,2)*(H_gp-H_i),dNdx(i_N,1)*(H_gp-H_i)]

          else
              BI_enr(1:3,1:2)=ZR
          end if
          B_XFEM(1:3,num_B_XFEM+1:num_B_XFEM+2) = BI_enr
          num_B_XFEM = num_B_XFEM + 2                  
      elseif (Enriched_Node_Type(c_NN(i_N),i_C).eq.3)then
          Jun_elem = Node_Jun_elem(c_NN(i_N),i_C)
          Coor_AB(1,:) = [Coors_Element_Crack(Jun_elem,i_C,1),Coors_Element_Crack(Jun_elem,i_C,2)]
          Coor_AB(2,:) = [Coors_Element_Crack(Jun_elem,i_C,3),Coors_Element_Crack(Jun_elem,i_C,4)]  
          Coor_CD(1,:) = [Coors_Junction(Jun_elem,i_C,1),Coors_Junction(Jun_elem,i_C,2)]
          Coor_CD(2,:) = [Coors_Junction(Jun_elem,i_C,3), Coors_Junction(Jun_elem,i_C,4)]
          x0 = [Coors_Junction(Jun_elem,i_C,1),Coors_Junction(Jun_elem,i_C,2)]
          
          call Cal_Signed_Distance(Coor_AB,[c_X_NODES(i_N),c_Y_NODES(i_N)],distance_Node_M)    
          call Cal_Signed_Distance(Coor_AB,Global_coor_Gauss,distance_Gauss_M)                      
          call Cal_Signed_Distance(Coor_CD,[c_X_NODES(i_N),c_Y_NODES(i_N)],distance_Node_sm)    
          call Cal_Signed_Distance(Coor_CD,Global_coor_Gauss,distance_Gauss_sm) 
          call Cal_Signed_Distance(Coor_AB,x0,distance_x0) 
          
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
          
          BI_enr(1,:)  = [dNdx(i_N,1)*(J_gp-J_i),     ZR]
          BI_enr(2,:)  = [ZR,      dNdx(i_N,2)*(J_gp-J_i)]
          BI_enr(3,:)  = [dNdx(i_N,2)*(J_gp-J_i),dNdx(i_N,1)*(J_gp-J_i)]
          B_XFEM(1:3,num_B_XFEM+1:num_B_XFEM+2) = BI_enr
          num_B_XFEM = num_B_XFEM + 2
      elseif (Enriched_Node_Type(c_NN(i_N),i_C).eq.6 .and. num_Circ_Hole>=1)then 
          Jun_elem = Node_Jun_elem(c_NN(i_N),i_C)
          c_Hole   = Node_Jun_Hole(c_NN(i_N),i_C) 
          x0_Hole =  Hole_Coor(c_Hole,1)
          y0_Hole =  Hole_Coor(c_Hole,2)
          R0_Hole =  Hole_Coor(c_Hole,3) 
          
          Coor_CD(1,:) = [Coors_Junction(Jun_elem,i_C,1),Coors_Junction(Jun_elem,i_C,2)]
          Coor_CD(2,:) = [Coors_Junction(Jun_elem,i_C,3), Coors_Junction(Jun_elem,i_C,4)]
          x0 = [Coors_Junction(Jun_elem,i_C,1),Coors_Junction(Jun_elem,i_C,2)]
          
          call Tool_Signed_Distance_Point_to_Circle(x0_Hole,y0_Hole,R0_Hole,[c_X_NODES(i_N),c_Y_NODES(i_N)],distance_Node_M)
          call Tool_Signed_Distance_Point_to_Circle(x0_Hole,y0_Hole,R0_Hole,Global_coor_Gauss,distance_Gauss_M)

          call Cal_Signed_Distance(Coor_CD,[c_X_NODES(i_N),c_Y_NODES(i_N)],distance_Node_sm)    
          call Cal_Signed_Distance(Coor_CD,Global_coor_Gauss,distance_Gauss_sm) 
        
          call Tool_Signed_Distance_Point_to_Circle(x0_Hole,y0_Hole,R0_Hole,x0,distance_x0)
          
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
          
          BI_enr(1,:)  = [dNdx(i_N,1)*(J_gp-J_i),     ZR]
          BI_enr(2,:)  = [ZR,      dNdx(i_N,2)*(J_gp-J_i)]
          BI_enr(3,:)  = [dNdx(i_N,2)*(J_gp-J_i),dNdx(i_N,1)*(J_gp-J_i)]
          B_XFEM(1:3,num_B_XFEM+1:num_B_XFEM+2) = BI_enr
          num_B_XFEM = num_B_XFEM + 2
      elseif (Enriched_Node_Type(c_NN(i_N),i_C).eq.6 .and. num_Ellip_Hole>=1)then 
          Jun_elem = Node_Jun_elem(c_NN(i_N),i_C)
          c_Hole   = Node_Jun_Hole(c_NN(i_N),i_C) 
          x0_Hole =  Ellip_Hole_Coor(c_Hole,1)
          y0_Hole =  Ellip_Hole_Coor(c_Hole,2)
          a_Hole  =  Ellip_Hole_Coor(c_Hole,3) 
          b_Hole  =  Ellip_Hole_Coor(c_Hole,4) 
          theta_Hole =  Ellip_Hole_Coor(c_Hole,5)
                  
          Coor_CD(1,:) = [Coors_Junction(Jun_elem,i_C,1),Coors_Junction(Jun_elem,i_C,2)]
          Coor_CD(2,:) = [Coors_Junction(Jun_elem,i_C,3), Coors_Junction(Jun_elem,i_C,4)]
          x0 = [Coors_Junction(Jun_elem,i_C,1),Coors_Junction(Jun_elem,i_C,2)]     
          
          call Tool_Signed_Distance_Point_to_Ellipse(x0_Hole,y0_Hole,a_Hole,b_Hole,theta_Hole,&
            [c_X_NODES(i_N),c_Y_NODES(i_N)],distance_Node_M)          
          call Tool_Signed_Distance_Point_to_Ellipse(x0_Hole,y0_Hole,a_Hole,b_Hole,theta_Hole,&
            Global_coor_Gauss,distance_Gauss_M)       

          call Cal_Signed_Distance(Coor_CD,[c_X_NODES(i_N),c_Y_NODES(i_N)],distance_Node_sm)    
          call Cal_Signed_Distance(Coor_CD,Global_coor_Gauss,distance_Gauss_sm) 
        
          call Tool_Signed_Distance_Point_to_Ellipse(x0_Hole,y0_Hole,a_Hole,b_Hole,theta_Hole,x0,distance_x0) 

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
          
          BI_enr(1,:)  = [dNdx(i_N,1)*(J_gp-J_i),     ZR]
          BI_enr(2,:)  = [ZR,      dNdx(i_N,2)*(J_gp-J_i)]
          BI_enr(3,:)  = [dNdx(i_N,2)*(J_gp-J_i),dNdx(i_N,1)*(J_gp-J_i)]
          B_XFEM(1:3,num_B_XFEM+1:num_B_XFEM+2) = BI_enr
          num_B_XFEM = num_B_XFEM + 2                  
      elseif (Enriched_Node_Type(c_NN(i_N),i_C).eq.4)then

      elseif (Enriched_Node_Type(c_NN(i_N),i_C).eq.1)then
          ref_elem = 0
          if (Elem_Type(i_E,i_C) .eq. 1) then
              ref_elem = i_E
          else
              find_element = Node_Elements(c_NN(i_N),:)
              
              do i_tem=1,num_Node_Elements(c_NN(i_N))
                  if (Elem_Type(find_element(i_tem),i_C) .eq.1) then
                      ref_elem = find_element(i_tem)
                      Flag_Check_Kink_Tip = .True.
                  end if
              end do
              
              if (ref_elem ==0) then  
                  ref_elem=Ele_Num_Tip_Enriched_Node(i_C,c_NN(i_N))
                  Flag_Check_Kink_Tip = .True.
              end if  
              if (ref_elem ==0) then
                  print *,'    WARNING::in Cal_B_Matrix_Crack!'   
              end if                       
          end if
          Coor_AB(1,:) = [Coors_Element_Crack(ref_elem,i_C,1),Coors_Element_Crack(ref_elem,i_C,2)]
          Coor_AB(2,:) = [Coors_Element_Crack(ref_elem,i_C,3),Coors_Element_Crack(ref_elem,i_C,4)]  
          c_Segment = Coor_AB(2,:) - Coor_AB(1,:)     
          omega   = atan2(c_Segment(2),c_Segment(1))            
                      
          Coor_Tip = [Coor_AB(2,1),Coor_AB(2,2)]
          Trans_Matrix(1,:) = [  cos(omega),sin(omega)]
          Trans_Matrix(2,:) = [ -sin(omega),cos(omega)]
          
          num_points_crack = Each_Cr_Poi_Num(i_C)
          if ((Flag_Check_Kink_Tip .eqv. .True.) .and. (num_points_crack > 2)) then
              if ((Elem_Type(i_E,i_C) .eq. 2) .or.(Elem_Type(i_E,i_C) .eq. 3))then
                   coor_A = Coor_Tip
                   tip_x  = x_cr_tip_nodes(i_C,c_NN(i_N))
                   tip_y  = y_cr_tip_nodes(i_C,c_NN(i_N))
                   if((Crack_Coor(i_C,1,1) .eq. tip_x) .and.  (Crack_Coor(i_C,1,2) .eq. tip_y)) then
                       coor_B =Crack_Coor(i_C,2,:)
                   elseif((Crack_Coor(i_C,num_points_crack,1) .eq. tip_x) .and.&
                          (Crack_Coor(i_C,num_points_crack,2) .eq. tip_y)) then    
                       coor_B = Crack_Coor(i_C,num_points_crack-1,:)
                   end if
                   if (coor_A(1).eq. Crack_Coor(i_C,1,1).and. coor_A(2) .eq.Crack_Coor(i_C,1,2)) then
                       coor_C = Crack_Coor(i_C,3,:)
                   elseif(coor_A(1).eq. Crack_Coor(i_C,num_points_crack,1) .and.&
                          coor_A(2).eq. Crack_Coor(i_C,num_points_crack,2)) then
                       coor_C =  Crack_Coor(i_C,num_points_crack-2,:)
                   end if 
                   if((Crack_Coor(i_C,1,1) .eq. tip_x) .and.(Crack_Coor(i_C,1,2) .eq. tip_y)) then
                       if(abs(Crack_Arc_Tip_A_B_C_x(i_C,1,3)) >=Tol_12)then
                        coor_A(1)=Crack_Arc_Tip_A_B_C_x(i_C,1,1)
                        coor_A(2)=Crack_Arc_Tip_A_B_C_y(i_C,1,1)
                        coor_B(1)=Crack_Arc_Tip_A_B_C_x(i_C,1,2)
                        coor_B(2)=Crack_Arc_Tip_A_B_C_y(i_C,1,2)     
                        coor_C(1)=Crack_Arc_Tip_A_B_C_x(i_C,1,3)
                        coor_C(2)=Crack_Arc_Tip_A_B_C_y(i_C,1,3)  
                       endif
                   elseif((Crack_Coor(i_C,num_points_crack,1).eq. tip_x) .and.&
                          (Crack_Coor(i_C,num_points_crack,2).eq. tip_y)) then    
                       if(abs(Crack_Arc_Tip_A_B_C_x(i_C,2,3))>=Tol_12)then
                        coor_A(1)=Crack_Arc_Tip_A_B_C_x(i_C,2,1)
                        coor_A(2)=Crack_Arc_Tip_A_B_C_y(i_C,2,1)
                        coor_B(1)=Crack_Arc_Tip_A_B_C_x(i_C,2,2)
                        coor_B(2)=Crack_Arc_Tip_A_B_C_y(i_C,2,2)     
                        coor_C(1)=Crack_Arc_Tip_A_B_C_x(i_C,2,3)
                        coor_C(2)=Crack_Arc_Tip_A_B_C_y(i_C,2,3)  
                       endif
                   end if
                   len_AB = sqrt((coor_A(1)-coor_B(1))**2+(coor_A(2)-coor_B(2))**2)
                   Line_AB(1,:) = coor_A
                   Line_AB(2,:) = coor_B
                   call Tool_Shorten_or_Extend_Line(Line_AB,TWO*len_AB,'B',new_Line_AB,coor_D)    
                   Tri_BCD(1,:) = coor_B
                   Tri_BCD(2,:) = coor_C
                   Tri_BCD(3,:) = coor_D
                   Closed_Tri_BCD(1:3,:) =  Tri_BCD
                   Closed_Tri_BCD(4,:)   =  Tri_BCD(1,:)
                   call Tool_Yes_In_Poly(Global_coor_Gauss(1),Global_coor_Gauss(2),&
                                         Closed_Tri_BCD(:,1), Closed_Tri_BCD(:,2),4,Yes_Gauss_in_BCD)        
                   if (Yes_Gauss_in_BCD) then
                       Line_BC(1,:) = coor_B
                       Line_BC(2,:) = coor_C                               
                       call Cal_Angle_of_AB_and_BC(Line_AB,Line_BC,angle_AB_BC)
                       Line_BGauss(1,:) = coor_B
                       Line_BGauss(2,:) = Global_coor_Gauss
                       call Cal_Angle_of_AB_and_BC(Line_AB,Line_BGauss,angle_AB_BGauss)
                       len_BGauss  = sqrt((coor_B(1)-Global_coor_Gauss(1))**2.0+(coor_B(2)-Global_coor_Gauss(2))**2.0)   
                       if (angle_AB_BGauss >= angle_AB_BC) then
                        thete_Doc= pi/TWO*(angle_AB_BGauss-angle_AB_BC)/(1.5D0*pi-angle_AB_BC)
                       else
                        thete_Doc= pi/TWO*(angle_AB_BGauss-angle_AB_BC)/(angle_AB_BC-HLF*pi)
                       end if
                       Global_coor_Gauss_Doc =[-len_AB-len_BGauss*cos(thete_Doc), -len_BGauss*sin(thete_Doc)]
                  end if
              end if                   
          end if
          if (Yes_Gauss_in_BCD.eqv..False.) then
              xp_Gauss= MATMUL(Trans_Matrix,Global_coor_Gauss-Coor_Tip)
              r_Gauss = sqrt(xp_Gauss(1)**2+xp_Gauss(2)**2)
              theta_Gauss= atan2(xp_Gauss(2),xp_Gauss(1))
          elseif(Yes_Gauss_in_BCD.eqv..True.) then
              r_Gauss     = sqrt(Global_coor_Gauss_Doc(1)**2 + Global_coor_Gauss_Doc(2)**2)
              theta_Gauss = atan2(Global_coor_Gauss_Doc(2),Global_coor_Gauss_Doc(1))
          end if 

          if ((theta_Gauss >pi  .or. theta_Gauss < -pi)) then
              print *,'    **************************************************************' 
              print *, '    Error :: Angle is wrong when calculates r and theta!'
              print *,'    **************************************************************'     
          end if       
          call Cal_F_dFdx_dFdy(r_Gauss,theta_Gauss,omega,c_mat_type,F_Gauss,dFdx,dFdy)
          xp_Node= MATMUL(Trans_Matrix,[c_X_NODES(i_N),c_Y_NODES(i_N)]-Coor_Tip)
          r_Node = sqrt(xp_Node(1)**2+xp_Node(2)**2)
          theta_Node= atan2(xp_Node(2),xp_Node(1))            
          call Cal_F(r_Node,theta_Node,omega,c_mat_type,F_Node)     
          

          num_BI_enr_Tip = 0
          if (Key_TipEnrich==1) then   
            do i_F =1,4
              c_N(1) = N(1,1)
              c_N(2) = N(1,3)
              c_N(3) = N(1,5)
              c_N(4) = N(1,7)     
              aa = dNdx(i_N,1)*(F_Gauss(i_F)-F_Node(i_F)) + c_N(i_N)*dFdx(i_F)
              bb = dNdx(i_N,2)*(F_Gauss(i_F)-F_Node(i_F)) + c_N(i_N)*dFdy(i_F)

              B_enr(1,:) = [aa,  ZR]
              B_enr(2,:) = [ZR,  bb]
              B_enr(3,:) = [bb,     aa]
              BI_enr_Tip(1:3,num_BI_enr_Tip+1:num_BI_enr_Tip+2) = B_enr  
              num_BI_enr_Tip = num_BI_enr_Tip +2
            end do
            B_XFEM(1:3,num_B_XFEM+1:num_B_XFEM+8) = BI_enr_Tip
            num_B_XFEM = num_B_XFEM + 8  
          elseif (Key_TipEnrich==2) then   
          
              c_N(1) = N(1,1)
              c_N(2) = N(1,3)
              c_N(3) = N(1,5)
              c_N(4) = N(1,7)     
              aa = dNdx(i_N,1)*(F_Gauss(1)-F_Node(1)) + c_N(i_N)*dFdx(1)
              bb = dNdx(i_N,2)*(F_Gauss(1)-F_Node(1)) + c_N(i_N)*dFdy(1)
              B_enr(1,:) = [aa,  ZR]
              B_enr(2,:) = [ZR,  bb]
              B_enr(3,:) = [bb,     aa] 
              num_BI_enr_Tip = num_BI_enr_Tip +2
              B_XFEM(1:3,num_B_XFEM+1:num_B_XFEM+2) = B_enr
              num_B_XFEM = num_B_XFEM + 2  
          elseif (Key_TipEnrich==3) then   
          elseif (Key_TipEnrich==4) then   
              c_N(1) = N(1,1)
              c_N(2) = N(1,3)
              c_N(3) = N(1,5)
              c_N(4) = N(1,7)     
              aa = dNdx(i_N,1)*(F_Gauss(1)-F_Node(1)) + c_N(i_N)*dFdx(1)
              bb = dNdx(i_N,2)*(F_Gauss(1)-F_Node(1)) + c_N(i_N)*dFdy(1)
              B_enr(1,:) = [aa,  ZR]
              B_enr(2,:) = [ZR,  bb]
              B_enr(3,:) = [bb,     aa] 
              num_BI_enr_Tip = num_BI_enr_Tip +2
              B_XFEM(1:3,num_B_XFEM+1:num_B_XFEM+2) = B_enr
              num_B_XFEM = num_B_XFEM + 2  
          endif 
      end if
  end do
  if (i_C.eq.1) then
      tem_B(1:3,1:8)            = B_FEM
      tem_B(1:3,9:8+num_B_XFEM) = B_XFEM(1:3,1:num_B_XFEM)
      num_tem_B = 8 + num_B_XFEM
  else
      tem_B(1:3,1:num_B_XFEM) = B_XFEM(1:3,1:num_B_XFEM)
      num_tem_B = num_B_XFEM
  end if
end if   

RETURN
END SUBROUTINE Cal_B_Matrix_Crack