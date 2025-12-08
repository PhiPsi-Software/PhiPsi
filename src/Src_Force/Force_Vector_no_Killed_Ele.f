 
      SUBROUTINE Force_Vector_no_Killed_Ele(c_Total_FD,
     &                             isub,Yes_Biot,Lambda,globalF)

      use Global_Float_Type
      use Global_Common
      use Global_Model
      use Global_Elem_Area_Vol
      use Global_Crack
      use Global_Crack_Common
      use Global_HF
      use Global_Material
      use Global_Stress
      use Global_Field_Problem
      use Global_Inclusion
            
      implicit none

      integer,intent(in)::isub,c_Total_FD
      real(kind=FT),intent(in)::Lambda
      logical,intent(in)::Yes_Biot
      real(kind=FT),intent(out)::globalF(c_Total_FD)
      real(kind=FT) th,c_density
      integer i,i_E,i_N,cur_Node,mat_num,c_mat_type,c_node
      
      
      real(kind=FT) c_thick,c_D(3,3),c_B(3,8)
      real(kind=FT) c_X_NODES(4),c_Y_NODES(4)
      integer c_NN(4) 
      integer c_Num_Gauss_Point,i_G
      
      real(kind=FT) c_kesi,c_yita
      integer local(8),i_row,i_col,nIndex,G_Counter
      real(kind=FT) c_Stress(3),c_temp(8),detJ
      
      real(kind=FT) globalF_T(c_Total_FD)
      integer i_C,i_Incl,i_H
      real(kind=FT) B(3,80),tem_B(3,80)
      integer num_B,num_tem_B
      integer:: Location_ESM(MDOF_2D)
      integer Location_ESM_C_Cr_NoFEM(60)
      integer num_Loc_ESM_C_Cr_NoFEM
      integer num_Loc_ESM
      integer::Location_ESM_C_Crack(80)
      integer num_Loc_ESM_C_Crack
      real(kind=FT) kesi_Enr(Num_Gauss_Points),
     &                 yita_Enr(Num_Gauss_Points),
     &                 weight_Enr(Num_Gauss_Points)
      real(kind=FT) kesi_N_Enr(4),
     &                 yita_N_Enr(4),weight_N_Enr(4)
      real(kind=FT) c_Inter_Force(80)
      real(kind=FT) kesi(900),yita(900),weight(900)
      real(kind=FT) Orient1,Orient2,x1,y1,x2,y2,ori_n(2),ori_CT_elem
      real(kind=FT) Div_Points(Max_Num_Cr_CalP,2),Coor_Tip(2),
     &                 r,Coor_AB(2,2),Orient_Tip,ele_Sign,c_Sign
      real(kind=FT) kesi_Enr_64(64),
     &              yita_Enr_64(64),
     &              weight_Enr_64(64)
      real(kind=FT) c_N(2,8),c_G_x,c_G_y
      logical Yes_Gauss_in_Incl
      integer c_Incl_Num,c_MatNum
      real(kind=FT) c_T_Alpha,c_TStress(3)
      real(kind=FT) c_v
      integer c_node_DOF,i_Boux_nonzero,i_Bouy_nonzero
      real(kind=FT) penalty_k,c_dis
      integer num_Ele_Killed
      integer iii
      
      print *,'    Constructing global force vector...'
      

      if (Key_Type_2D ==1) then
          if (Material_Type(1) ==1) then
              th = Material_Para(1,4)
            else
                th = Material_Para(1,4)
            end if
      elseif (Key_Type_2D ==2) then
          th =ONE
      end if
      
      globalF(1:c_Total_FD) = ZR
      do i = 1,Num_Foc_x
            cur_Node = int(Foc_x(i,1))
          globalF(2*cur_Node-1) = Lambda*Foc_x(i,2)*(1.0D0)
      end do
      do i = 1,Num_Foc_y
            cur_Node = int(Foc_y(i,1))
            globalF(2*cur_Node)=Lambda*Foc_y(i,2)*(1.0D0)
      end do

      if(Key_Gravity==1) then
          do i_E = 1,Num_Elem
                mat_num    = Elem_Mat(i_E)
                c_mat_type = Material_Type(mat_num)
              if (c_mat_type ==1) then
                  c_density = Material_Para(mat_num,3)
              else
            c_density = Material_Para(mat_num,3); 
              end if

                do i_N = 1,4
                  c_node = G_NN(i_N,i_E)
                  globalF(2*c_node-1)  = globalF(2*c_node-1)-
     &                       g_X_Y_Z(1)*c_density*th*Elem_Area(i_E)/FOU
              end do
                do i_N = 1,4
                  c_node = G_NN(i_N,i_E) 
                  globalF(2*c_node)  = globalF(2*c_node)-
     &                       g_X_Y_Z(2)*c_density*th*Elem_Area(i_E)/FOU
              end do
            end do
      end if
      
      if(Key_Thermal_Stress==1)then
          if (Key_Integral_Sol  == 2)then
              call Cal_Gauss_Points_QUAD(Num_Gauss_Points,
     &                                 kesi_Enr,
     &                                 yita_Enr,
     &                                 weight_Enr)
              call Cal_Gauss_Points_QUAD(64,kesi_Enr_64,yita_Enr_64,
     &                                 weight_Enr_64)
              call Cal_Gauss_Points_QUAD(Num_Gauss_P_FEM,
     &                                 kesi_N_Enr,
     &                                 yita_N_Enr,
     &                                 weight_N_Enr)
          elseif (Key_Integral_Sol  == 3)then
              call Cal_Gauss_Points_QUAD_for_SUBQUAD(Num_Sub_Quads,
     &                               kesi_Enr,yita_Enr,
     &                               weight_Enr)
              Num_Gauss_Points = Num_Sub_Quads*4
              Num_Gauss_P_Inc  = Num_Sub_Quads*4
              call Cal_Gauss_Points_QUAD(64,kesi_Enr_64,yita_Enr_64,
     &                             weight_Enr_64)
              call Cal_Gauss_Points_QUAD(Num_Gauss_P_FEM,kesi_N_Enr,
     &                               yita_N_Enr,weight_N_Enr)
          endif
      
          globalF_T(1:c_Total_FD)  = ZR
          EleGaus_yes_FEM_asemd(1:Num_Elem,1:Num_Gauss_Points)= .false.
          do i_E = 1,Num_Elem 
              c_thick = thick(Elem_Mat(i_E))
              c_D     = D(Elem_Mat(i_E),:,:)     
              
              if(Flag_Weibull_E)then
                  if (Key_Weibull_E(Elem_Mat(i_E)) ==1)then
                      c_D = Weibull_Elements_D_Matrix(i_E,1:3,1:3)
                  endif
              endif
          
              c_v     = v(Elem_Mat(i_E),1)

              c_NN    = G_NN(:,i_E)
              c_X_NODES = G_X_NODES(:,i_E)
              c_Y_NODES = G_Y_NODES(:,i_E)                
              Location_ESM(1:MDOF_2D)   = 0
              num_Loc_ESM           = 0
              if(num_Crack/=0)then
                do i_C =1,num_Crack
                  call Location_Element_Stiff_Matrix(i_E,i_C,
     &                                        c_POS(:,i_C),
     &                                        Location_ESM_C_Crack,
     &                                        num_Loc_ESM_C_Crack,
     &                                        Location_ESM_C_Cr_NoFEM,
     &                                        num_Loc_ESM_C_Cr_NoFEM)
                  Location_ESM(num_Loc_ESM+1:
     &                      num_Loc_ESM+num_Loc_ESM_C_Crack) = 
     &                      Location_ESM_C_Crack(1:num_Loc_ESM_C_Crack)
                  num_Loc_ESM  =  num_Loc_ESM + num_Loc_ESM_C_Crack
                end do
              endif
              if(num_Hole/=0)then
                do i_H =1,num_Hole
                  call Location_Element_Stiff_Matrix_Hl(i_E,i_H,
     &                                        c_POS_Hl(:,i_H),
     &                                        Location_ESM_C_Crack,
     &                                        num_Loc_ESM_C_Crack,
     &                                        Location_ESM_C_Cr_NoFEM,
     &                                        num_Loc_ESM_C_Cr_NoFEM)
                  Location_ESM(num_Loc_ESM+1:
     &                      num_Loc_ESM+num_Loc_ESM_C_Crack) = 
     &                      Location_ESM_C_Crack(1:num_Loc_ESM_C_Crack)
                  num_Loc_ESM  =  num_Loc_ESM + num_Loc_ESM_C_Crack
                end do
              endif
              if(num_Inclusion/=0)then
                do i_Incl =1,num_Inclusion
                  call Location_Element_Stiff_Matrix_Incl(i_E,i_Incl,
     &                                        c_POS_Incl(:,i_Incl),
     &                                        Location_ESM_C_Crack,
     &                                        num_Loc_ESM_C_Crack,
     &                                        Location_ESM_C_Cr_NoFEM,
     &                                        num_Loc_ESM_C_Cr_NoFEM)
                  Location_ESM(num_Loc_ESM+1:
     &                      num_Loc_ESM+num_Loc_ESM_C_Crack) = 
     &                      Location_ESM_C_Crack(1:num_Loc_ESM_C_Crack)
                  num_Loc_ESM  =  num_Loc_ESM + num_Loc_ESM_C_Crack
                end do
              endif
              if(Key_Integral_Sol.eq.1)then
              elseif(Key_Integral_Sol.eq.2 .or.
     &               Key_Integral_Sol.eq.3)then
                  if (num_Crack/= 0 .and. 
     &                (maxval(Enriched_Node_Type(c_NN,
     &                             1:num_Crack)).ne.0))then
                      kesi(1:Num_Gauss_Points)    = kesi_Enr
                      yita(1:Num_Gauss_Points)    = yita_Enr
                      weight(1:Num_Gauss_Points)  = weight_Enr
                      c_Num_Gauss_Point = Num_Gauss_Points
                  elseif(num_Hole/= 0 .and.
     &                (maxval(Enriched_Node_Type_Hl(
     &                      c_NN,1:num_Hole)).ne.0))then
                      kesi(1:Num_Gauss_Points)    = kesi_Enr
                      yita(1:Num_Gauss_Points)    = yita_Enr
                      weight(1:Num_Gauss_Points)  = weight_Enr
                      c_Num_Gauss_Point = Num_Gauss_Points
                  elseif(num_Inclusion/= 0 .and.
     &                (maxval(Enriched_Node_Type_Incl
     &                              (c_NN,1:num_Inclusion)).ne.0))
     &                                                            then
                      kesi(1:Num_Gauss_Points)    = kesi_Enr
                      yita(1:Num_Gauss_Points)    = yita_Enr
                      weight(1:Num_Gauss_Points)  = weight_Enr
                      c_Num_Gauss_Point = Num_Gauss_Points
                  else 
                      kesi(1:Num_Gauss_P_FEM)    = kesi_N_Enr
                      yita(1:Num_Gauss_P_FEM)    = yita_N_Enr
                      weight(1:Num_Gauss_P_FEM)  = weight_N_Enr
                      c_Num_Gauss_Point = Num_Gauss_P_FEM
                  end if 
              endif
              do i_G = 1,c_Num_Gauss_Point  
                  c_D     = D(Elem_Mat(i_E),:,:) 
                  
                  if(Flag_Weibull_E)then
                      if (Key_Weibull_E(Elem_Mat(i_E)) ==1)then
                          c_D = Weibull_Elements_D_Matrix(i_E,1:3,1:3)
                      endif
                  endif
          
                  c_v     = v(Elem_Mat(i_E),1)
                  c_T_Alpha = T_Alpha(Elem_Mat(i_E))
                  B(1:3,1:80) = ZR
                  num_B = 0
                  call Cal_N(kesi(i_G),yita(i_G),c_N)
                  call Cal_detJ(kesi(i_G),yita(i_G),
     &                          c_X_NODES,c_Y_NODES,detJ)  
                  c_G_x = DOT_PRODUCT(c_N(1,1:7:2),c_X_NODES(1:4))
                  c_G_y = DOT_PRODUCT(c_N(1,1:7:2),c_Y_NODES(1:4))   

                  if(num_Crack/=0)then
                    do i_C =1,num_Crack 
                      call Cal_B_Matrix_Crack(kesi(i_G),yita(i_G),
     &                                    i_C,i_E,i_G,
     &                                    c_NN,c_X_NODES,c_Y_NODES,
     &                                    tem_B,num_tem_B)
                      B(1:3,num_B+1:num_B+num_tem_B) =  
     &                                          tem_B(1:3,1:num_tem_B)
                      num_B = num_B + num_tem_B
                    end do
                  endif
                  if(num_Hole/=0)then
                      do i_H =1,num_Hole 
                          call Cal_B_Matrix_Hl(kesi(i_G),yita(i_G),
     &                                        i_H,i_E,i_G,
     &                                        c_NN,c_X_NODES,c_Y_NODES,
     &                                        tem_B,num_tem_B)
                          B(1:3,num_B+1:num_B+num_tem_B) =  
     &                                        tem_B(1:3,1:num_tem_B)
                          num_B = num_B + num_tem_B
                      end do
                  endif
                  if(num_Circ_Incl/=0)then
                      do i_Incl =1,num_Circ_Incl 
                          call Cal_B_Matrix_Circ_Incl(kesi(i_G),
     &                                       yita(i_G),i_Incl,i_E,i_G,
     &                                        c_NN,c_X_NODES,c_Y_NODES,
     &                                        tem_B,num_tem_B)
                          B(1:3,num_B+1:num_B+num_tem_B) =  
     &                                        tem_B(1:3,1:num_tem_B)
                          num_B = num_B + num_tem_B
                      end do
                      call Tool_Yes_Point_in_Inclusions(c_G_x,c_G_y,
     &                                Yes_Gauss_in_Incl,c_Incl_Num)
                      if(Yes_Gauss_in_Incl)then
                          c_MatNum = Circ_Inclu_Mat_Num(c_Incl_Num)
                          c_D     = D(c_MatNum,:,:)   
                          c_T_Alpha = T_Alpha(c_MatNum)
                          c_v = v(c_MatNum,1)
                      endif
                  endif
                  if(num_Poly_Incl/=0)then
                      do i_Incl =1,num_Poly_Incl
                          call Cal_B_Matrix_Poly_Incl(kesi(i_G),
     &                                       yita(i_G),i_Incl,i_E,i_G,
     &                                        c_NN,c_X_NODES,c_Y_NODES,
     &                                        tem_B,num_tem_B)
                          B(1:3,num_B+1:num_B+num_tem_B) =  
     &                                        tem_B(1:3,1:num_tem_B)
                          num_B = num_B + num_tem_B
                      end do
                      call Tool_Yes_Point_in_Inclusions(c_G_x,c_G_y,
     &                                Yes_Gauss_in_Incl,c_Incl_Num)
                      if(Yes_Gauss_in_Incl)then
                          c_MatNum = Poly_Inclu_Mat_Num(c_Incl_Num)
                          c_D     = D(c_MatNum,:,:)   
                          c_T_Alpha = T_Alpha(c_MatNum)
                          c_v = v(c_MatNum,1)
                      endif
                  endif
                  c_TStress(1:3) =ZR
                  if(Key_Type_2D==1)then
                      c_TStress = c_T_Alpha*
     &                         Thermal_Str_Temper(Elem_Mat(i_E))*
     &                         MATMUL(c_D,[ONE,ONE,ZR]) 
                  elseif(Key_Type_2D==2)then
                      c_TStress = (ONE+c_v)*
     &                         c_T_Alpha*
     &                         Thermal_Str_Temper(Elem_Mat(i_E))*
     &                         MATMUL(c_D,[ONE,ONE,ZR]) 
                  endif
                  c_Inter_Force(1:num_Loc_ESM) = 
     &                   MATMUL(transpose(B(1:3,1:num_Loc_ESM)),
     &                          c_TStress)
                  globalF_T(Location_ESM(1:num_Loc_ESM)) = 
     &               globalF_T(Location_ESM(1:num_Loc_ESM))  + 
     &               c_Inter_Force(1:num_Loc_ESM)*
     &               c_thick*weight(i_G)*detJ
              end do
          enddo  
          globalF = globalF + globalF_T
      endif
      
      if(Key_Analysis_Type==17 .and. Yes_Biot)then  
          globalF = globalF - MATMUL(Biot_c_MAT,Porous_P)  
      endif

      if(Num_Boux_nonzero/=0  .or.  Num_Bouy_nonzero/=0) then
          penalty_k =  penalty_k_bou_nonzero 
          do i = 1,Num_Boux_nonzero
              c_node = int(Bou_x_nonzero(i,1))
              c_dis  = Bou_x_nonzero(i,2)
              c_node_DOF = 2*c_node-1
              globalF(c_node_DOF)= globalF(c_node_DOF) + penalty_k*c_dis
            end do
          do i = 1,Num_Bouy_nonzero
              c_node = int(Bou_y_nonzero(i,1))
              c_dis  = Bou_y_nonzero(i,2)
              c_node_DOF = 2*c_node
              globalF(c_node_DOF)= globalF(c_node_DOF) + penalty_k*c_dis
            end do          
      end if

  199 continue  
      
      RETURN
      END SUBROUTINE Force_Vector_no_Killed_Ele
