 
      SUBROUTINE Force_Vector(c_Total_FD,isub,Yes_Biot,Lambda,globalF)

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
      integer local(8),i_row,i_col,G_Counter
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
      real(kind=FT) x1,y1,x2,y2,ori_n(2)
      real(kind=FT) ori_CT_elem
      real(kind=FT)Coor_Tip(2)
      real(kind=FT) r,Coor_AB(2,2),Orient_Tip,ele_Sign,c_Sign
      real(kind=FT) kesi_Enr_64(64),
     &              yita_Enr_64(64),
     &              weight_Enr_64(64)
      real(kind=FT) c_N(2,8),c_G_x,c_G_y
      logical Yes_Gauss_in_Incl
      integer c_Incl_Num,c_MatNum
      real(kind=FT) c_T_Alpha,c_TStress(3)
      real(kind=FT) c_v
      integer c_node_DOF
      real(kind=FT) penalty_k,c_dis
      integer num_Ele_Killed
      integer Local_W(60),num_Local_W,Num_Div_Points  
      real(kind=FT) CT_Length,tem_n(4),N_W(2,60),N(2,8)
      real(kind=FT) c_PC_Gauss(2)
      integer i_CT_GAUSS
      integer Num_Enriched_Node
      integer ref_elem,i_tem,find_element(8)
      real(kind=FT) tem_val1,tem_val2,tem_val3,tem_val4
      integer j_E,c_Adj_Ele,c_Junc_Ele
      integer j_C
      integer num_CT_Gauss,i_CT_Elem
      real(kind=FT) Contact_Force_x,Contact_Force_y
      real(kind=FT) kesi_CT(2),Weight_CT(2)
      real(kind=FT) N_CT(2)
      real(kind=FT) CT_GAUSS_x,CT_GAUSS_y
      integer CT_GAUSS_Elem
      integer c_Ele
      real(kind=FT) c_InSitu_x,c_InSitu_y
      real(kind=FT) ori_n_l,ori_n_m
      real(kind=FT) tem_F(Num_Node*2)
                  
      
      print *,'    Constructing global force vector...'
      
      if(Key_EKILL==1)then
          num_Ele_Killed = count(Ele_Killed_Each_Load_Step(isub,:)>0)
      endif
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
              
              if(Key_EKILL==1 .and. num_Ele_Killed>=1) then
                if(any(Ele_Killed_Each_Load_Step
     &                          (1:isub,1:num_Ele_Killed)==i_E))then
                    c_density = ZR
                endif
              endif  
              
              
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
      if(Key_Thermal_Stress==1 .and. Yes_XFEM .eqv. .True.)then
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
                      c_TStress = 
     &                     c_T_Alpha*Thermal_Str_Temper(Elem_Mat(i_E))*
     &                         MATMUL(c_D,[ONE,ONE,ZR]) 
                  elseif(Key_Type_2D==2)then
                      c_TStress = (ONE+c_v)*
     &                      c_T_Alpha*Thermal_Str_Temper(Elem_Mat(i_E))*
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

      if((Key_Thermal_Stress==1) .and. (Yes_XFEM .eqv. .False.))then
          call Cal_Gauss_Points_QUAD(Num_Gauss_P_FEM,kesi,yita,weight)
          G_Counter = 0
          do i_E = 1,Num_Elem
              c_thick   = thick(Elem_Mat(i_E))
              c_D       = D(Elem_Mat(i_E),:,:) 
              
              if(Flag_Weibull_E)then
                  if (Key_Weibull_E(Elem_Mat(i_E)) ==1)then
                      c_D = Weibull_Elements_D_Matrix(i_E,1:3,1:3)
                  endif
              endif
          
              c_T_Alpha = T_Alpha(Elem_Mat(i_E))
              c_v       = v(Elem_Mat(i_E),1)
                
              c_NN    = G_NN(:,i_E)
              c_X_NODES = G_X_NODES(:,i_E)
              c_Y_NODES = G_Y_NODES(:,i_E)                 
              local=[c_NN(1)*2-1,c_NN(1)*2,c_NN(2)*2-1,c_NN(2)*2,
     &               c_NN(3)*2-1,c_NN(3)*2,c_NN(4)*2-1,c_NN(4)*2]   
        
              do i_G = 1,Num_Gauss_P_FEM  
                  G_Counter = G_Counter +1
                  c_kesi = kesi(i_G)
                  c_yita = yita(i_G)
                  call Cal_detJ(c_kesi,c_yita,c_X_NODES,c_Y_NODES,detJ)  
                  call Cal_Ele_B_N4(i_E,c_X_NODES,c_Y_NODES,
     &                              c_kesi,c_yita,c_B)
                  c_TStress(1:3) =ZR
                  if(Key_Type_2D==1)then
                      c_TStress = 
     &                     c_T_Alpha*Thermal_Str_Temper(Elem_Mat(i_E))*
     &                         MATMUL(c_D,[ONE,ONE,ZR]) 
                  elseif(Key_Type_2D==2)then
                      c_TStress = (ONE+c_v)*
     &                      c_T_Alpha*Thermal_Str_Temper(Elem_Mat(i_E))*
     &                         MATMUL(c_D,[ONE,ONE,ZR]) 
                  endif
                  c_temp(1:8) =MATMUL(transpose(c_B(1:3,1:8)),c_TStress)
          
          
                  
                  globalF(local) =  globalF(local)
     &                       +c_temp(1:8)*c_thick*weight(i_G)*detJ
              end do            
          enddo
      endif


      if(Key_PoreP==1)then
        if(Key_Analysis_Type==1)then  
            Porous_P(1:Num_Node) = -Initial_PoreP
            tem_F(1:Num_Node*2) = MATMUL(Biot_c_MAT,Porous_P)
            if(Key_EKILL==1)then
              do i_E = 1,Num_Elem
                if(any(Ele_Killed_Each_Load_Step(1:isub,
     &                                   1:num_Ele_Killed) ==i_E))then
                    
                    tem_F(G_NN(1:4,i_E)*2)   = ZR
                    tem_F(G_NN(1:4,i_E)*2-1) = ZR
                endif  
              enddo
            endif
            globalF = globalF - tem_F          
        endif
      endif
      
      if(Key_FS_Seque_Coupling==1)then
        Porous_P(1:Num_Node) = -Fd_Value(1:Num_Node)
        tem_F(1:Num_Node*2) = MATMUL(Biot_c_MAT,Porous_P)
        if(Key_EKILL==1)then
          do i_E = 1,Num_Elem
            if(any(Ele_Killed_Each_Load_Step(1:isub,1:num_Ele_Killed)
     &                                                      ==i_E))then
                
                tem_F(G_NN(1:4,i_E)*2)   = ZR
                tem_F(G_NN(1:4,i_E)*2-1) = ZR
            endif  
          enddo
        endif
        globalF = globalF - tem_F
      endif  
      
      if(Key_Analysis_Type==17 .and. Yes_Biot)then  
          globalF = globalF - MATMUL(Biot_c_MAT,Porous_P)  
      endif
        
      if(Key_InSitu_Strategy==2)then
          if(abs(InSitu_x) <0.0001D6 .and. abs(InSitu_y) < 0.0001D6)then
              print *,'    #--#--#--#--#--#--#--#--#--#--#--#--#--#'
              print *,'    Warning :: cannot find InSitu_x and '
     &                   // 'InSitu_y (warning in Force_Vector.f)!'
              print *,'    #--#--#--#--#--#--#--#--#--#--#--#--#--#'
              goto 199
          endif
          call Cal_Gauss_Points_QUAD(Num_Gauss_P_FEM,kesi,yita,weight)
          G_Counter = 0
          do i_E = 1,Num_Elem
              c_thick = thick(Elem_Mat(i_E))
              c_D     = D(Elem_Mat(i_E),:,:)  
              
              if(Flag_Weibull_E)then
                  if (Key_Weibull_E(Elem_Mat(i_E)) ==1)then
                      c_D = Weibull_Elements_D_Matrix(i_E,1:3,1:3)
                  endif
              endif
          
              if(Key_EKILL==1 .and. num_Ele_Killed>=1) then
                if(any(
     &             Ele_Killed_Each_Load_Step(1:isub,1:num_Ele_Killed)
     &                                                    ==i_E))then
                    cycle
                endif
              endif  
              
              c_NN    = G_NN(:,i_E)
              c_X_NODES = G_X_NODES(:,i_E)
              c_Y_NODES = G_Y_NODES(:,i_E)   
              local=[c_NN(1)*2-1,c_NN(1)*2,c_NN(2)*2-1,c_NN(2)*2,
     &               c_NN(3)*2-1,c_NN(3)*2,c_NN(4)*2-1,c_NN(4)*2]       
              do i_G = 1,Num_Gauss_P_FEM  
                  G_Counter = G_Counter +1
                  c_kesi = kesi(i_G); c_yita = yita(i_G)
                  call Cal_detJ(c_kesi,c_yita,c_X_NODES,c_Y_NODES,detJ)  
                  call Cal_Ele_B_N4(i_E,c_X_NODES,c_Y_NODES,
     &                              c_kesi,c_yita,c_B) 
                  c_Stress = [InSitu_Strs_Gaus_xx(i_E,i_G),
     &                        InSitu_Strs_Gaus_yy(i_E,i_G),
     &                        InSitu_Strs_Gaus_xy(i_E,i_G)]
                  c_temp(1:8) = MATMUL(transpose(c_B(1:3,1:8)),c_Stress)
                  globalF(local) = globalF(local)  - 
     &               c_temp(1:8)*c_thick*weight(i_G)*detJ
              end do            
          enddo
              
          if(Key_Analysis_Type==1 .and. num_Crack/=0 .and.
     &       Key_Force_Control /= 5)then 
            num_CT_Gauss = 2
            kesi_CT(1)  = -0.577350269189626D0
            kesi_CT(2)  =  0.577350269189626D0
            Weight_CT(1)=  ONE
            Weight_CT(2)=  ONE            
            do i_C=1,num_Crack
              if(Key_Crack_Inner_Pressure==1  .and. 
     &          (Crack_Pressure(i_C)>Tol_10)) then
                  CONTINUE
              else
                  cycle
              endif
              Num_Div_Points = Cracks_CalP_Num(i_C)
              do i_CT_Elem = 1,Num_Div_Points - 1
                  N_W(1:2,1:60)  = ZR   
                  x1= Cracks_CalP_Coors(i_C,i_CT_Elem,1)
                  y1= Cracks_CalP_Coors(i_C,i_CT_Elem,2)
                  x2= Cracks_CalP_Coors(i_C,i_CT_Elem+1,1)
                  y2= Cracks_CalP_Coors(i_C,i_CT_Elem+1,2)
                  ori_CT_elem = atan2(y2-y1,x2-x1)
                  CT_Length = Cracks_HF_Ele_L(i_C,i_CT_Elem)
                  detJ   = CT_Length/TWO
                  ori_n = [cos(pi/TWO + ori_CT_elem),
     &                     sin(pi/TWO + ori_CT_elem)]
                  
                  c_Ele = Cracks_CalP_Elem(i_C,i_CT_Elem) 
                  c_InSitu_x = 
     &               sum(InSitu_Strs_Gaus_xx(c_Ele,1:Num_Gauss_P_FEM))/
     &                                    Num_Gauss_P_FEM
                  c_InSitu_y = 
     &               sum(InSitu_Strs_Gaus_yy(c_Ele,1:Num_Gauss_P_FEM))/
     &                                    Num_Gauss_P_FEM   
                  ori_n_l = ori_n(1)
                  ori_n_m = ori_n(2)
                  Contact_Force_x =c_InSitu_x * ori_n_l
                  Contact_Force_y =c_InSitu_y * ori_n_m
                  
                  do i_CT_GAUSS  = 1,num_CT_Gauss
                      call Cal_HF_Coor_by_Kesi(kesi_CT(i_CT_GAUSS),
     &                                         x1,y1,x2,y2,
     &                                         CT_GAUSS_x,CT_GAUSS_y)
                      call Cal_Ele_Num_by_Coors(CT_GAUSS_x,CT_GAUSS_y,
     &                                          CT_GAUSS_Elem)
                      call Cal_KesiYita_by_Coor([CT_GAUSS_x,CT_GAUSS_y],
     &                                    CT_GAUSS_Elem,c_Kesi,c_Yita)
                      c_NN    = G_NN(1:4,CT_GAUSS_Elem)     
                      call Cal_N(c_Kesi,c_Yita,N)      
                      tem_n(1) = N(1,1);    tem_n(2) = N(1,3)
                      tem_n(3) = N(1,5);    tem_n(4) = N(1,7)      
                      N_CT(1) = (ONE-kesi_CT(i_CT_GAUSS))/TWO      
                      N_CT(2) = (ONE+kesi_CT(i_CT_GAUSS))/TWO   
                      num_Local_W = 0
                      N_W(1:2,1:60) =ZR
                      do i_N = 1,4
                          if (Enriched_Node_Type(c_NN(i_N),i_C) ==2)then
                              num_Local_W = num_Local_W+1
                              Local_W(2*num_Local_W-1)=
     &                                          2*c_POS(c_NN(i_N),i_C)-1
                              Local_W(2*num_Local_W)  =
     &                                          2*c_POS(c_NN(i_N),i_C)
                              N_W(1,2*num_Local_W-1) = TWO*tem_n(i_N)
                              N_W(2,2*num_Local_W)   = TWO*tem_n(i_N)
                          elseif(Enriched_Node_Type(c_NN(i_N),i_C) ==3)
     &                                                              then
                              c_Junc_Ele = Node_Jun_elem(c_NN(i_N),i_C)
                              Coor_Tip=
     &                                Coors_Junction(c_Junc_Ele,i_C,3:4) 
                              Orient_Tip = atan2(Coor_Tip(2)-CT_GAUSS_y,
     &                                           Coor_Tip(1)-CT_GAUSS_x)
                              ele_Sign =  sign(ONE,cos(ori_CT_elem))
                              c_Sign   =  sign(ONE,cos(Orient_Tip))
                              tem_val1 = abs(ori_CT_elem - pi/TWO)
                              tem_val2 = abs(ori_CT_elem + pi/TWO)
                              if (tem_val1 <= Tol_10) then
                                  ele_Sign = ONE
                                  tem_val3 = abs(Orient_Tip - pi/TWO)
                                  tem_val4 = abs(Orient_Tip + pi/TWO)
                                  if(tem_val4 <= Tol_10) then
                                      c_Sign   = -ONE
                                  elseif (tem_val3 <= Tol_10) then
                                      c_Sign   =  ONE
                                  endif
                              elseif (tem_val2 <= Tol_10) then
                                  ele_Sign = ONE
                                  tem_val3 = abs(Orient_Tip - pi/TWO)
                                  tem_val4 = abs(Orient_Tip + pi/TWO)  
                                  if(tem_val4 <= Tol_10) then
                                      c_Sign   =  ONE
                                  elseif (tem_val3 <= Tol_10) then
                                      c_Sign   = -ONE
                                  endif                          
                              endif
                              num_Local_W = num_Local_W+1
                              Local_W(2*num_Local_W-1)=
     &                                          2*c_POS(c_NN(i_N),i_C)-1
                              Local_W(2*num_Local_W)  =
     &                                          2*c_POS(c_NN(i_N),i_C)
                              N_W(1,2*num_Local_W-1)=c_Sign*ele_Sign*
     &                                               TWO*tem_n(i_N)
                              N_W(2,2*num_Local_W)  =c_Sign*ele_Sign*
     &                                               TWO*tem_n(i_N)
                          elseif(Enriched_Node_Type(c_NN(i_N),i_C) ==1)
     &                                                             then
                              Num_Enriched_Node = c_POS(c_NN(i_N),i_C)
                              if(Elem_Type(CT_GAUSS_Elem,i_C).eq.1)then
                                  ref_elem = CT_GAUSS_Elem
                              else
                                  find_element = 
     &                                       Node_Elements(c_NN(i_N),:)
                                  do i_tem=
     &                                    1,num_Node_Elements(c_NN(i_N))
                                      if (Elem_Type(find_element(i_tem),
     &                                                i_C) .eq.1)  then
                                          ref_elem = find_element(i_tem)
                                      end if
                                  end do
                              end if
                              Coor_AB(1,:)=
     &                             [Coors_Element_Crack(ref_elem,i_C,1),
     &                              Coors_Element_Crack(ref_elem,i_C,2)]
                              Coor_AB(2,:)=
     &                             [Coors_Element_Crack(ref_elem,i_C,3),
     &                              Coors_Element_Crack(ref_elem,i_C,4)]
                              Coor_Tip = [Coor_AB(2,1),Coor_AB(2,2)]
                              r=sqrt((CT_GAUSS_x-Coor_Tip(1))**2 + 
     &                               (CT_GAUSS_y-Coor_Tip(2))**2)
                              Orient_Tip = atan2(Coor_Tip(2)-CT_GAUSS_y,
     &                                           Coor_Tip(1)-CT_GAUSS_x)
                              num_Local_W = num_Local_W+1
                              Local_W(2*num_Local_W-1)=
     &                                 2*(c_POS(c_NN(i_N),i_C)) - 1
                              Local_W(2*num_Local_W)  =
     &                                 2*(c_POS(c_NN(i_N),i_C))
                              ele_Sign =  sign(ONE,cos(ori_CT_elem))
                              c_Sign   =  sign(ONE,cos(Orient_Tip))
                              tem_val1 = abs(ori_CT_elem - pi/TWO)
                              tem_val2 = abs(ori_CT_elem + pi/TWO)
                              if (tem_val1 <= Tol_10) then
                                  ele_Sign = ONE
                                  tem_val3 = abs(Orient_Tip - pi/TWO)
                                  tem_val4 = abs(Orient_Tip + pi/TWO)
                                  if(tem_val4 <= Tol_10) then
                                      c_Sign   = -ONE
                                  elseif (tem_val3 <= Tol_10) then
                                      c_Sign   =  ONE
                                  endif
                              elseif (tem_val2 <= Tol_10) then
                                  ele_Sign = ONE
                                  tem_val3 = abs(Orient_Tip - pi/TWO)
                                  tem_val4 = abs(Orient_Tip + pi/TWO)  
                                  if(tem_val4 <= Tol_10) then
                                      c_Sign   =  ONE
                                  elseif (tem_val3 <= Tol_10) then
                                      c_Sign   = -ONE
                                  endif                          
                              endif
                              N_W(1,2*num_Local_W-1) = 
     &                            ele_Sign*c_Sign*TWO*sqrt(r)*tem_n(i_N)
                              N_W(2,2*num_Local_W)   = 
     &                            ele_Sign*c_Sign*TWO*sqrt(r)*tem_n(i_N)
                          endif
                          if(Enriched_Node_Type(c_NN(i_N),i_C)==2)then
                            do j_C = 1,num_Crack
                              if(Enriched_Node_Type(c_NN(i_N),j_C)==3)
     &                                                             then
                                  do j_E =1,num_Node_Elements(c_NN(i_N))
                                      c_Adj_Ele=
     &                                      Node_Elements(c_NN(i_N),j_E) 
                                      if(Elem_Type(c_Adj_Ele,j_C)==4)
     &                                                             then
                                          c_Junc_Ele = c_Adj_Ele 
                                          exit
                                      endif
                                  enddo
                                  Coor_Tip =Coors_Junction(
     &                                             c_Junc_Ele,i_C,3:4) 
                                  Orient_Tip = atan2(
     &                                          Coor_Tip(2)-CT_GAUSS_y,
     &                                          Coor_Tip(1)-CT_GAUSS_x)
                                  ele_Sign =  sign(ONE,cos(ori_CT_elem))
                                  c_Sign   =  sign(ONE,cos(Orient_Tip))
                                  tem_val1 = abs(ori_CT_elem - pi/TWO)
                                  tem_val2 = abs(ori_CT_elem + pi/TWO)
                                  if (tem_val1 <= Tol_10) then
                                      ele_Sign = ONE
                                      tem_val3 =abs(Orient_Tip - pi/TWO)
                                      tem_val4 =abs(Orient_Tip + pi/TWO)
                                      if(tem_val4 <= Tol_10) then
                                          c_Sign   = -ONE
                                      elseif (tem_val3 <= Tol_10) then
                                          c_Sign   =  ONE
                                      endif
                                  elseif (tem_val2 <= Tol_10) then
                                      ele_Sign = ONE
                                      tem_val3 =abs(Orient_Tip - pi/TWO)
                                      tem_val4 =abs(Orient_Tip + pi/TWO)  
                                      if(tem_val4 <= Tol_10) then
                                          c_Sign   =  ONE
                                      elseif (tem_val3 <= Tol_10) then
                                          c_Sign   = -ONE
                                      endif                          
                                  endif
                                  num_Local_W = num_Local_W+1
                                  Local_W(2*num_Local_W-1)=
     &                                         2*c_POS(c_NN(i_N),j_C)-1
                                  Local_W(2*num_Local_W)  =
     &                                         2*c_POS(c_NN(i_N),j_C)
                                  N_W(1,2*num_Local_W-1)= 
     &                                  -c_Sign*ele_Sign*ONE*tem_n(i_N)
                                  N_W(2,2*num_Local_W)  = 
     &                                  -c_Sign*ele_Sign*ONE*tem_n(i_N)
                              endif
                            enddo
                          endif
                      enddo
                      c_PC_Gauss = [Contact_Force_x,Contact_Force_y]
                      globalF(Local_W(1:2*num_Local_W)) = 
     &                       globalF(Local_W(1:2*num_Local_W))  + 
     &                 MATMUL(transpose(N_W(1:2,1:2*num_Local_W)),
     &                    c_PC_Gauss)* Weight_CT(i_CT_GAUSS)*detJ
                  enddo
              enddo
            enddo
          endif
      endif

      if(Key_InStress_for_Mat==1)then
          call Cal_Gauss_Points_QUAD(Num_Gauss_P_FEM,kesi,yita,weight)
          G_Counter = 0
          do i_E = 1,Num_Elem
              c_thick = thick(Elem_Mat(i_E))
              c_D     = D(Elem_Mat(i_E),:,:)  
              
              if(Flag_Weibull_E)then
                  if (Key_Weibull_E(Elem_Mat(i_E)) ==1)then
                      c_D = Weibull_Elements_D_Matrix(i_E,1:3,1:3)
                  endif
              endif
          
              if (any(Mat_Number_of_InStress == Elem_Mat(i_E)))then     
                  c_NN    = G_NN(:,i_E)
                  c_X_NODES = G_X_NODES(:,i_E)
                  c_Y_NODES = G_Y_NODES(:,i_E)   
                  local=[c_NN(1)*2-1,c_NN(1)*2,c_NN(2)*2-1,c_NN(2)*2,
     &                   c_NN(3)*2-1,c_NN(3)*2,c_NN(4)*2-1,c_NN(4)*2]       
                  do i_G = 1,Num_Gauss_P_FEM  
                      G_Counter = G_Counter +1
                      c_kesi = kesi(i_G); c_yita = yita(i_G)
                      call Cal_detJ(c_kesi,c_yita,c_X_NODES,
     &                                            c_Y_NODES,detJ)  
                      call Cal_Ele_B_N4(i_E,c_X_NODES,c_Y_NODES,
     &                                  c_kesi,c_yita,c_B) 
                      c_Stress = [Mat_InStress_x,Mat_InStress_y,ZR]
                      c_temp(1:8) = MATMUL(transpose(c_B(1:3,1:8)),
     &                                     c_Stress)
                      globalF(local) = 
     &                   globalF(local)  - 
     &                   c_temp(1:8)*c_thick*weight(i_G)*detJ
                  end do   
              end if
          enddo 
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
      END SUBROUTINE Force_Vector
