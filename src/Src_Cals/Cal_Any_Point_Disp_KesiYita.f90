 
subroutine Cal_Any_Point_Disp_KesiYita(elem,Kesi,Yita,i_G,c_DISP,P_Disp)

use Global_Float_Type
use Global_Common
use Global_Crack
use Global_Crack_Common
use Global_Model
use Global_Elem_Area_Vol
use Global_DISP
use Global_Inclusion
use Global_Cross

implicit none
integer,intent(in)::elem,i_G
real(kind=FT),intent(in)::Kesi,Yita,c_DISP(Total_FD)
real(kind=FT),intent(out)::P_Disp(2)
real(kind=FT) X_NODES(4),Y_NODES(4),&
              detJ,dNdkesi(4,2),J(2,2),N(2,8),&
              Coor_AB(2,2),Point(2)
integer c_Elem,i_C,i_N,NODES_iE(4),Num_Enri_Node,ref_elem
real(kind=FT) c_Segment(2),&
              omega, Coor_Tip(2),Trans_Matrix(2,2),  &
              tem_N(4),coor_A(2),coor_B(2),&
              coor_C(2),len_AB,line_AB(2,2),&
              Tri_BCD(3,2),Closed_Tri_BCD(4,2),&
              tip_x,tip_y,xp_Gauss(2),r_Gauss,&
              Line_BC(2,2),thete_Doc,Global_coor_Gauss_Doc(2),&
              len_BGauss,Line_BGauss(2,2),angle_AB_BGauss,&
              angle_AB_BC,new_Line_AB(2,2),coor_D(2)
integer find_element(8)
integer i_tem,num_points_crack,mat_num,c_mat_type
logical Flag_Check_Kink_Tip,Yes_Gauss_in_BCD
real(kind=FT) theta_Gauss,Disp_x,Disp_y
logical Yes_in_KAB
real(kind=FT) Coor_CD(2,2),x0(2),&
              distance_Node_M,distance_Gauss_M,&
              distance_Node_sm,distance_Gauss_sm,&
              distance_x0,H_gp_M,H_gp_sm,H_i_M,H_i_sm,&
              H_x0, J_gp, J_i
real(kind=FT) distance_Vertex,&
        Tri_KAB(3,2),Closed_Tri_KAB(4,2),H_i, &
        H_Vertex,distance_Gauss,H_gp,distance_Node, &
        xp_Node(2),r_Node,theta_Node,F,&
        F_Gauss(4),F_Node(4)
integer i_F
integer Jun_elem,i_H
real(kind=FT) c_Hole_x,c_Hole_y,c_Hole_r
real(kind=FT) Hl_i,c_G_x,c_G_y,Tool_Function_2Point_Dis
real(kind=FT) c_Dis_G,c_Dis_N,Hl_gp
integer c_NODE,i_Incl
real(kind=FT) c_Incl_x,c_Incl_y,c_Incl_r
real(kind=FT) Zeta_Node(4),Zeta_Node_abs(4)
real(kind=FT) Zm,Za,c_N(4),Incl_gp
integer n_Vertex
integer i_Cross
integer c_Cross_elem
real(kind=FT) Cross_gp,Cross_i
integer c_Hole
real(kind=FT) x0_Hole,y0_Hole,R0_Hole
real(kind=FT) a_Hole,b_Hole,theta_Hole

c_Elem = elem

NODES_iE = G_NN(1:4,c_Elem)

X_NODES = G_X_NODES(:,c_Elem)
Y_NODES = G_Y_NODES(:,c_Elem)

call Cal_Coor_by_KesiYita(Kesi,Yita,X_NODES,Y_NODES,point(1),point(2))
c_G_x = point(1)
c_G_y = point(2)

call Cal_N_dNdkesi_J_detJ(Kesi,Yita,X_NODES,Y_NODES,detJ,J,N,dNdkesi)
tem_N(1:4) = N(1,1:7:2)

Disp_x = c_DISP(NODES_iE(1)*2-1)*N(1,1) + c_DISP(NODES_iE(2)*2-1)*N(1,3) + &
       c_DISP(NODES_iE(3)*2-1)*N(1,5) + c_DISP(NODES_iE(4)*2-1)*N(1,7)

Disp_y = c_DISP(NODES_iE(1)*2)*N(1,1) + c_DISP(NODES_iE(2)*2)*N(1,3) + &
       c_DISP(NODES_iE(3)*2)*N(1,5) + c_DISP(NODES_iE(4)*2)*N(1,7)

mat_num    = Elem_Mat(c_Elem)

c_mat_type = Material_Type(mat_num)

if(Num_Crack /=0)then
    do i_C=1,Num_Crack
      if (maxval(Enriched_Node_Type(NODES_iE,i_C)).eq. 0) then
          goto 1000
      end if

      do i_N = 1,4
          Flag_Check_Kink_Tip =  .False.
          Yes_Gauss_in_BCD    =  .False.
          if (Enriched_Node_Type(NODES_iE(i_N),i_C) ==2) then
              Num_Enri_Node = c_POS(NODES_iE(i_N),i_C)

              if ((Elem_Type(c_Elem,i_C) == 2) .or.  (Elem_Type(c_Elem,i_C) == 3)) then
                  ref_elem = c_Elem
                  Coor_AB(1,:)=[Coors_Element_Crack(ref_elem,i_C,1),Coors_Element_Crack(ref_elem,i_C,2)]
                  Coor_AB(2,:)=[Coors_Element_Crack(ref_elem,i_C,3),Coors_Element_Crack(ref_elem,i_C,4)]

                  if  (Elem_Type(c_Elem,i_C) .eq. 3) then
                      call Cal_Signed_Distance(Coor_AB,Coors_Vertex(ref_elem,:),distance_Vertex)
                      call Cal_Sign(distance_Vertex,H_Vertex)
                      call Cal_Signed_Distance(Coor_AB,Point,distance_Gauss)
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
                          Closed_Tri_KAB(4,:)   =  Tri_KAB(1,:)
                          call Tool_Yes_In_Poly(Point(1),Point(2),Closed_Tri_KAB(:,1), Closed_Tri_KAB(:,2),4,Yes_in_KAB)
                          if (Yes_in_KAB .eqv. .True.) then
                              H_gp = -H_gp
                          end if
                      end if
                      call Cal_Signed_Distance(Coor_AB,[X_NODES(i_N),Y_NODES(i_N)],distance_Node)
                      call Cal_Sign(distance_Node,H_i)
                      if(Key_Heaviside_Value==0)then
                          if (H_i<=ZR) H_i = ZR
                      endif

                  else
                      call Cal_Signed_Distance(Coor_AB,Point,distance_Gauss)
                      call Cal_Sign(distance_Gauss,H_gp)
                      if(Key_Heaviside_Value==0)then
                          if (H_gp<=ZR) H_gp = ZR
                      endif
                      call Cal_Signed_Distance(Coor_AB,[X_NODES(i_N),Y_NODES(i_N)],distance_Node)
                      call Cal_Sign(distance_Node,H_i)
                      if(Key_Heaviside_Value==0)then
                          if (H_i<=ZR) H_i = ZR
                      endif
                  end if
                  Disp_x = Disp_x+(H_gp-H_i)*tem_N(i_N)*c_DISP(Num_Enri_Node*2-1)
                  Disp_y = Disp_y+(H_gp-H_i)*tem_N(i_N)*c_DISP(Num_Enri_Node*2)
              end if
          elseif (Enriched_Node_Type(NODES_iE(i_N),i_C) ==3) then
              Jun_elem = Node_Jun_elem(NODES_iE(i_N),i_C)
              Num_Enri_Node = c_POS(NODES_iE(i_N),i_C)
              Coor_AB(1,:) = [Coors_Element_Crack(Jun_elem,i_C,1),Coors_Element_Crack(Jun_elem,i_C,2)]
              Coor_AB(2,:) = [Coors_Element_Crack(Jun_elem,i_C,3),Coors_Element_Crack(Jun_elem,i_C,4)]
              Coor_CD(1,:) = [Coors_Junction(Jun_elem,i_C,1),Coors_Junction(Jun_elem,i_C,2)]
              Coor_CD(2,:) = [Coors_Junction(Jun_elem,i_C,3), Coors_Junction(Jun_elem,i_C,4)]
              x0 = [Coors_Junction(Jun_elem,i_C,1),Coors_Junction(Jun_elem,i_C,2)]
              call Cal_Signed_Distance(Coor_AB,[X_NODES(i_N),Y_NODES(i_N)],distance_Node_M)
              call Cal_Signed_Distance(Coor_AB,Point,distance_Gauss_M)
              call Cal_Signed_Distance(Coor_CD,[X_NODES(i_N),Y_NODES(i_N)],distance_Node_sm)
              call Cal_Signed_Distance(Coor_CD,Point,distance_Gauss_sm)
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
              Disp_x = Disp_x + (J_gp-J_i)*tem_N(i_N)*c_DISP(Num_Enri_Node*2-1)
              Disp_y = Disp_y + (J_gp-J_i)*tem_N(i_N)*c_DISP(Num_Enri_Node*2)
          elseif (Enriched_Node_Type(NODES_iE(i_N),i_C) ==6 .and. num_Circ_Hole>=1)then
              Jun_elem = Node_Jun_elem(NODES_iE(i_N),i_C)
              c_Hole   = Node_Jun_Hole(NODES_iE(i_N),i_C)
              x0_Hole =  Hole_Coor(c_Hole,1)
              y0_Hole =  Hole_Coor(c_Hole,2)
              R0_Hole =  Hole_Coor(c_Hole,3)
              Num_Enri_Node = c_POS(NODES_iE(i_N),i_C)

              Coor_CD(1,:) = [Coors_Junction(Jun_elem,i_C,1),Coors_Junction(Jun_elem,i_C,2)]
              Coor_CD(2,:) = [Coors_Junction(Jun_elem,i_C,3), Coors_Junction(Jun_elem,i_C,4)]
              x0 = [Coors_Junction(Jun_elem,i_C,1), Coors_Junction(Jun_elem,i_C,2)]

              call Tool_Signed_Distance_Point_to_Circle(x0_Hole,y0_Hole,R0_Hole,[X_NODES(i_N),Y_NODES(i_N)],distance_Node_M)
              call Tool_Signed_Distance_Point_to_Circle(x0_Hole,y0_Hole,R0_Hole,Point,distance_Gauss_M)

              call Cal_Signed_Distance(Coor_CD,[X_NODES(i_N),Y_NODES(i_N)],distance_Node_sm)
              call Cal_Signed_Distance(Coor_CD,Point,distance_Gauss_sm)

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
              Disp_x = Disp_x + (J_gp-J_i)*tem_N(i_N)*c_DISP(Num_Enri_Node*2-1)
              Disp_y = Disp_y + (J_gp-J_i)*tem_N(i_N)*c_DISP(Num_Enri_Node*2)
          elseif (Enriched_Node_Type(NODES_iE(i_N),i_C) ==6 .and. num_Ellip_Hole>=1)then
              Jun_elem = Node_Jun_elem(NODES_iE(i_N),i_C)
              c_Hole   = Node_Jun_Hole(NODES_iE(i_N),i_C)
              x0_Hole =  Ellip_Hole_Coor(c_Hole,1)
              y0_Hole =  Ellip_Hole_Coor(c_Hole,2)
              a_Hole  =  Ellip_Hole_Coor(c_Hole,3)
              b_Hole  =  Ellip_Hole_Coor(c_Hole,4)
              theta_Hole =  Ellip_Hole_Coor(c_Hole,5)

              Num_Enri_Node = c_POS(NODES_iE(i_N),i_C)

              Coor_CD(1,:) = [Coors_Junction(Jun_elem,i_C,1),Coors_Junction(Jun_elem,i_C,2)]
              Coor_CD(2,:) = [Coors_Junction(Jun_elem,i_C,3),Coors_Junction(Jun_elem,i_C,4)]
              x0 = [Coors_Junction(Jun_elem,i_C,1),Coors_Junction(Jun_elem,i_C,2)]

              call Tool_Signed_Distance_Point_to_Ellipse(x0_Hole,y0_Hole,a_Hole,b_Hole,&
                                                         theta_Hole,[X_NODES(i_N),Y_NODES(i_N)],distance_Node_M)
              call Tool_Signed_Distance_Point_to_Ellipse(x0_Hole,y0_Hole,a_Hole,b_Hole,theta_Hole,Point,distance_Gauss_M)

              call Cal_Signed_Distance(Coor_CD,[X_NODES(i_N),Y_NODES(i_N)],distance_Node_sm)
              call Cal_Signed_Distance(Coor_CD,Point,distance_Gauss_sm)

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
              Disp_x = Disp_x + (J_gp-J_i)*tem_N(i_N)*c_DISP(Num_Enri_Node*2-1)
              Disp_y = Disp_y + (J_gp-J_i)*tem_N(i_N)*c_DISP(Num_Enri_Node*2)
          elseif (Enriched_Node_Type(NODES_iE(i_N),i_C) ==1) then
              ref_elem = 0
              if (Elem_Type(c_Elem,i_C) .eq. 1) then
                  ref_elem = c_Elem
              else
                  find_element = Node_Elements(NODES_iE(i_N),:)

                  do i_tem=1,num_Node_Elements(NODES_iE(i_N))
                      if (Elem_Type(find_element(i_tem),i_C) .eq.1) then
                          ref_elem = find_element(i_tem)
                          Flag_Check_Kink_Tip = .True.
                      end if
                  end do
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
                    if ((Elem_Type(c_Elem,i_C) .eq. 2 ) .or.(Elem_Type(c_Elem,i_C) .eq. 3))then
                        coor_A = Coor_Tip
                        tip_x = x_cr_tip_nodes(i_C,NODES_iE(i_N))
                        tip_y = y_cr_tip_nodes(i_C,NODES_iE(i_N))
                        if((Crack_Coor(i_C,1,1).eq. tip_x).and. (Crack_Coor(i_C,1,2) .eq. tip_y)) then
                             coor_B =Crack_Coor(i_C,2,:)
                        elseif((Crack_Coor(i_C,num_points_crack,1) .eq. tip_x) .and. (Crack_Coor(i_C,num_points_crack,2) .eq. &
                            tip_y)) then
                           coor_B =Crack_Coor(i_C,num_points_crack-1,:)
                        end if
                        len_AB = sqrt((coor_A(1)-coor_B(1))**2 +(coor_A(2)-coor_B(2))**2)
                        if (coor_A(1).eq. Crack_Coor(i_C,1,1).and.coor_A(2) .eq.Crack_Coor(i_C,1,2)) then
                            coor_C = Crack_Coor(i_C,3,:)
                        elseif(coor_A(1).eq.Crack_Coor(i_C,num_points_crack,1) .and. &
                               coor_A(2).eq. Crack_Coor(i_C,num_points_crack, 2)) then
                            coor_C=Crack_Coor(i_C,num_points_crack-2,:)
                        end if
                       Line_AB(1,:) = coor_A
                       Line_AB(2,:) = coor_B
                       call Tool_Shorten_or_Extend_Line(Line_AB,2*len_AB,'B',new_Line_AB,coor_D)
                       Tri_BCD(1,:) = coor_B
                       Tri_BCD(2,:) = coor_C
                       Tri_BCD(3,:) = coor_D
                       Closed_Tri_BCD(1:3,:) =  Tri_BCD
                       Closed_Tri_BCD(4,:)   =  Tri_BCD(1,:)
                       call Tool_Yes_In_Poly(Point(1),Point(2),Closed_Tri_BCD(:,1), Closed_Tri_BCD(:,2),4,Yes_Gauss_in_BCD)
                       if (Yes_Gauss_in_BCD) then
                           Line_BC(1,:) = coor_B
                           Line_BC(2,:) = coor_C
                           call Cal_Angle_of_AB_and_BC(Line_AB,Line_BC,angle_AB_BC)
                           Line_BGauss(1,:) = coor_B
                           Line_BGauss(2,:) = Point
                           call Cal_Angle_of_AB_and_BC(Line_AB,Line_BGauss,angle_AB_BGauss)
                           len_BGauss  = sqrt((coor_B(1)-Point(1))**2.0+(coor_B(2)-Point(2))**2.0)
                           if (angle_AB_BGauss >= angle_AB_BC) then
                               thete_Doc= pi/TWO* (angle_AB_BGauss-angle_AB_BC)/(1.5D0*pi-angle_AB_BC)
                           else
                               thete_Doc= pi/TWO*(angle_AB_BGauss-angle_AB_BC)/(angle_AB_BC-HLF*pi)
                           end if
                           Global_coor_Gauss_Doc = [-len_AB-len_BGauss*cos(thete_Doc),-len_BGauss*sin(thete_Doc)]
                        end if
                    end if
              end if
        if (Yes_Gauss_in_BCD.eqv..False.) then
                  xp_Gauss= MATMUL(Trans_Matrix,Point-Coor_Tip)
                  r_Gauss     =  sqrt(xp_Gauss(1)**2+xp_Gauss(2)**2)
                  theta_Gauss = atan2(xp_Gauss(2),xp_Gauss(1))
              else
                  r_Gauss     = sqrt(Global_coor_Gauss_Doc(1)**2 + Global_coor_Gauss_Doc(2)**2)
                  theta_Gauss = atan2(Global_coor_Gauss_Doc(2),Global_coor_Gauss_Doc(1))
              end if
              call Cal_F(r_Gauss,theta_Gauss,omega,c_mat_type,F_Gauss)
              xp_Node= MATMUL(Trans_Matrix,[X_NODES(i_N),Y_NODES(i_N)]-Coor_Tip)
              r_Node = sqrt(xp_Node(1)**2+xp_Node(2)**2)
              theta_Node= atan2(xp_Node(2),xp_Node(1))
              call Cal_F(r_Node,theta_Node,omega,c_mat_type,F_Node)

              if (Key_TipEnrich==1) then
                  do i_F =1,4
                      Num_Enri_Node=c_POS(NODES_iE(i_N),i_C)+i_F-1
                      F = F_Gauss(i_F)-F_Node(i_F)
                      Disp_x = Disp_x + F*tem_N(i_N)*c_DISP(Num_Enri_Node*2-1)
                      Disp_y = Disp_y + F*tem_N(i_N)*c_DISP(Num_Enri_Node*2)
                  end do

              elseif (Key_TipEnrich==2) then
                  Num_Enri_Node=c_POS(NODES_iE(i_N),i_C)
                  F = F_Gauss(1)-F_Node(1)
                  Disp_x = Disp_x + F*tem_N(i_N)*c_DISP(Num_Enri_Node*2-1)
                  Disp_y = Disp_y + F*tem_N(i_N)*c_DISP(Num_Enri_Node*2)
              elseif (Key_TipEnrich==3) then
              elseif (Key_TipEnrich==4) then
                  Num_Enri_Node=c_POS(NODES_iE(i_N),i_C)
                  F = F_Gauss(1)-F_Node(1)
                  Disp_x = Disp_x + F*tem_N(i_N)* c_DISP(Num_Enri_Node*2-1)
                  Disp_y = Disp_y + F*tem_N(i_N)* c_DISP(Num_Enri_Node*2)
              endif
          end if
      end do
    1000   continue
    end do
endif

if(Num_Hole /=0)then
    do i_H =1,Num_Hole
      if (maxval(Enriched_Node_Type_Hl(NODES_iE,i_H)).eq. 0) then
          goto 2000
      end if
      do i_N = 1,4
          if (Enriched_Node_Type_Hl(NODES_iE(i_N),i_H).eq.1)then
              Num_Enri_Node = c_POS_Hl(NODES_iE(i_N),i_H)
            if (Elem_Type_Hl(c_Elem,i_H).eq.1) then
                  c_NODE  = Elem_Node(c_Elem,i_N)
                  c_Hole_x = Hole_Coor(i_H,1)
                  c_Hole_y = Hole_Coor(i_H,2)
                  c_Hole_r = Hole_Coor(i_H,3)
                  c_Dis_G=Tool_Function_2Point_Dis([c_G_x,c_G_y],[c_Hole_x,c_Hole_y])
                  c_Dis_N=Tool_Function_2Point_Dis(Coor(c_NODE,1:2),[c_Hole_x,c_Hole_y])
                  if (c_Dis_G<=c_Hole_r)then
                      if(Key_Hole_Value==0)then
                          Hl_gp = ZR
                      elseif(Key_Hole_Value==-1)then
                          Hl_gp = -ONE
                      endif
                  else
                      Hl_gp = ONE
                  endif
                  if (c_Dis_N<=c_Hole_r)then
                      if(Key_Hole_Value==0)then
                          Hl_i = ZR
                      elseif(Key_Hole_Value==-1)then
                          Hl_i = -ONE
                      endif
                  else
                      Hl_i = ONE
                  endif
                  Disp_x = Disp_x+(Hl_gp-Hl_i)*tem_N(i_N)*c_DISP(Num_Enri_Node*2-1)
                  Disp_y = Disp_y+(Hl_gp-Hl_i)*tem_N(i_N)*c_DISP(Num_Enri_Node*2)
              endif
          endif
      end do
    2000    continue
    enddo
endif

if(Num_Cross /=0)then
    do i_Cross =1,Num_Cross
      if(maxval(Enriched_Node_Type_Cross(NODES_iE,i_Cross)).eq.0)then
          goto 2100
      end if
      do i_N = 1,4
          if (Enriched_Node_Type_Cross(NODES_iE(i_N),i_Cross).eq.1)then
              Num_Enri_Node = c_POS_Cross(NODES_iE(i_N),i_Cross)
              c_Cross_elem = Node_Cross_elem(NODES_iE(i_N),i_Cross)
            if (Elem_Type_Cross(c_Cross_elem,i_Cross).eq.1) then
                  Coor_AB(1,1:2)=Cross_Point_RABCD(i_Cross,2,1:2)
                  Coor_AB(2,1:2)=Cross_Point_RABCD(i_Cross,3,1:2)
                  Coor_CD(1,1:2)=Cross_Point_RABCD(i_Cross,4,1:2)
                  Coor_CD(2,1:2)=Cross_Point_RABCD(i_Cross,5,1:2)
                  call Cal_Signed_Distance(Coor_AB,[X_NODES(i_N),Y_NODES(i_N)],distance_Node_M)
                  call Cal_Signed_Distance(Coor_AB,[c_G_x,c_G_y],distance_Gauss_M)
                  call Cal_Signed_Distance(Coor_CD,[X_NODES(i_N),Y_NODES(i_N)],distance_Node_sm)
                  call Cal_Signed_Distance(Coor_CD,[c_G_x,c_G_y],distance_Gauss_sm)
                  if(Key_Heaviside_Value==-1)then
                      call Cal_Sign(distance_Gauss_M*distance_Gauss_sm,Cross_gp)
                      call Cal_Sign(distance_Node_M*distance_Node_sm,Cross_i)
                  elseif(Key_Heaviside_Value==0)then
                      call Cal_Sign_1_and_0(distance_Gauss_M*distance_Gauss_sm,Cross_gp)
                      call Cal_Sign_1_and_0(distance_Node_M*distance_Node_sm,Cross_i)
                  endif

                  Disp_x = Disp_x+(Cross_gp-Cross_i)*tem_N(i_N)*c_DISP(Num_Enri_Node*2-1)
                  Disp_y = Disp_y+(Cross_gp-Cross_i)*tem_N(i_N)*c_DISP(Num_Enri_Node*2)
              endif
          endif
      end do
    2100     continue
    enddo
endif

if(num_Circ_Incl /=0)then
    do i_Incl =1,num_Circ_Incl
      if(maxval(Enriched_Node_Type_Incl(NODES_iE,i_Incl)).eq.0)then
          goto 3000
      end if
      Zeta_Node(1:4) = ZR
      do i_N = 1,4
          c_NODE  = Elem_Node(c_Elem,i_N)
          c_Incl_x = Circ_Inclu_Coor(i_Incl,1)
          c_Incl_y = Circ_Inclu_Coor(i_Incl,2)
          c_Incl_r = Circ_Inclu_Coor(i_Incl,3)
          c_Dis_N=Tool_Function_2Point_Dis(Coor(c_NODE,1:2), [c_Incl_x,c_Incl_y])
          Zeta_Node(i_N) = c_Dis_N-c_Incl_r
          if(abs(Zeta_Node(i_N))<1.0D-8)then
              Zeta_Node(i_N) = ZR
          endif
      enddo
      Zeta_Node_abs = abs(Zeta_Node)
      c_N(1) = N(1,1)
      c_N(2) = N(1,3)
      c_N(3) = N(1,5)
      c_N(4) = N(1,7)
      Zm = Zeta_Node(1)*c_N(1) + Zeta_Node(2)*c_N(2)+Zeta_Node(3)*c_N(3) + Zeta_Node(4)*c_N(4)
      Za = Zeta_Node_abs(1)*c_N(1)+Zeta_Node_abs(2)*c_N(2)+Zeta_Node_abs(3)*c_N(3)+Zeta_Node_abs(4)*c_N(4)

      do i_N = 1,4
          if(Enriched_Node_Type_Incl(NODES_iE(i_N),i_Incl).eq.1)then
              Num_Enri_Node = c_POS_Incl(NODES_iE(i_N),i_Incl)
            if (Elem_Type_Incl(c_Elem,i_Incl).eq.1) then
                  Incl_gp = Za-abs(Zm)
                  Disp_x = Disp_x+(Incl_gp)*tem_N(i_N)*c_DISP(Num_Enri_Node*2-1)
                  Disp_y = Disp_y+(Incl_gp)*tem_N(i_N)*c_DISP(Num_Enri_Node*2)
              endif
          endif
      end do
    3000    continue
    enddo
endif

if(num_Poly_Incl /=0)then
    do i_Incl =1,num_Poly_Incl
      n_Vertex = Poly_Inclu_Edges_Num(i_Incl)
      if (maxval(Enriched_Node_Type_Incl(NODES_iE,i_Incl)).eq.0)then
          goto 4000
      end if
      Zeta_Node(1:4) = ZR
      do i_N = 1,4
          c_NODE  = Elem_Node(c_Elem,i_N)
          call Tool_Dis_Point_to_Poly(Coor(c_NODE,1),Coor(c_NODE,2),&
                     Poly_Incl_Coor_x_Cl(i_Incl,1:n_Vertex+1),&
                     Poly_Incl_Coor_y_Cl(i_Incl,1:n_Vertex+1),n_Vertex+1,Zeta_Node(i_N))
          if(abs(Zeta_Node(i_N))<1.0D-8)then
              Zeta_Node(i_N) = ZR
          endif
      enddo
      Zeta_Node_abs = abs(Zeta_Node)
      c_N(1) = N(1,1)
      c_N(2) = N(1,3)
      c_N(3) = N(1,5)
      c_N(4) = N(1,7)
      Zm = Zeta_Node(1)*c_N(1) + Zeta_Node(2)*c_N(2)+Zeta_Node(3)*c_N(3) + Zeta_Node(4)*c_N(4)
      Za = Zeta_Node_abs(1)*c_N(1)+Zeta_Node_abs(2)*c_N(2)+Zeta_Node_abs(3)*c_N(3)+Zeta_Node_abs(4)*c_N(4)

      do i_N = 1,4
          if(Enriched_Node_Type_Incl(NODES_iE(i_N),i_Incl).eq.1)then
              Num_Enri_Node = c_POS_Incl(NODES_iE(i_N),i_Incl)
            if (Elem_Type_Incl(c_Elem,i_Incl).eq.1) then
                  Incl_gp = Za-abs(Zm)
                  Disp_x = Disp_x+(Incl_gp)*tem_N(i_N)*c_DISP(Num_Enri_Node*2-1)
                  Disp_y = Disp_y+(Incl_gp)*tem_N(i_N)*c_DISP(Num_Enri_Node*2)
              endif
          endif
      end do
    4000    continue
    enddo
endif

P_Disp(1) = Disp_x
P_Disp(2) = Disp_y

return
end SUBROUTINE Cal_Any_Point_Disp_KesiYita
