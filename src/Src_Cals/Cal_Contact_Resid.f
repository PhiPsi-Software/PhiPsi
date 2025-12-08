 
      subroutine Cal_Contact_Resid(
     &       iter,ifra,Counter_Iter,i_NR_P,
     &       c_Total_FD,c_num_freeDOF,F_U,
     &       c_DISP,freeDOF,PC_Gauss_x,PC_Gauss_y,R_PSI)
      

      use Global_Float_Type
      use Global_Common
      use Global_Model
      use Global_Crack
      use Global_Crack_Common
      use Global_HF
      use Global_Material
      
      implicit none
      integer, intent(in)::iter,ifra,Counter_Iter,i_NR_P
      integer, intent(in)::c_Total_FD,c_num_freeDOF
      real(kind=FT),intent(in)::F_U(c_Total_FD),c_DISP(c_Total_FD)
      integer, intent(in)::freeDOF(c_Total_FD)
      real(kind=FT),intent(in)::
     &              PC_Gauss_x(num_Crack,Max_Num_Cr_CalP-1,2),
     &              PC_Gauss_y(num_Crack,Max_Num_Cr_CalP-1,2)
      real(kind=FT),intent(out)::R_PSI(c_Total_FD)
      integer i_E,i_N,i_C,i_G
      real(kind=FT) c_thick,c_D(3,3),U(80)
      integer c_NN(4),c_Num_Gauss_Point
      real(kind=FT) c_X_NODES(4),c_Y_NODES(4),
     &                 c_kesi,c_yita,c_Stress(3),B(3,80),
     &                 tem_B(3,80)
      integer num_B,num_tem_B
      integer:: Location_ESM(200)
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
      real(kind=FT) detJ,c_Inter_Force(80)
      real(kind=FT) kesi(Num_Gauss_Points),yita(Num_Gauss_Points),
     &                 weight(Num_Gauss_Points)  
      real(kind=FT) kesi_CT(2),Weight_CT(2)
      real(kind=FT) N_CT(2)
      real(kind=FT) CT_GAUSS_x,CT_GAUSS_y
      integer CT_GAUSS_Elem
      integer num_CT_Gauss,i_CT_Elem
      integer Local_W(60),num_Local_W,Num_Div_Points  
      real(kind=FT) x1,y1,x2,y2,ori_n(2),ori_CT_elem
      real(kind=FT) CT_Length,tem_n(4),N_W(2,60),N(2,8)
      real(kind=FT) c_PC_Gauss(2)
      integer i_CT_GAUSS
      integer Num_Enriched_Node
      integer ref_elem,i_tem,find_element(8)
      real(kind=FT) Coor_Tip(2),
     &                 r,Coor_AB(2,2),Orient_Tip,ele_Sign,c_Sign
      real(kind=FT) tem_val1,tem_val2,tem_val3,tem_val4
      integer j_E,c_Adj_Ele,c_Junc_Ele
      integer j_C
      
      R_PSI(1:c_Total_FD)= ZR
      if (Key_Integral_Sol  == 2)then
          call Cal_Gauss_Points_QUAD(Num_Gauss_Points,
     &                     kesi_Enr,yita_Enr,
     &                     weight_Enr)
      elseif (Key_Integral_Sol  == 3)then
          call Cal_Gauss_Points_QUAD_for_SUBQUAD(Num_Sub_Quads,
     &                     kesi_Enr,yita_Enr,
     &                     weight_Enr)
          Num_Gauss_Points = Num_Sub_Quads*4
      endif
      call Cal_Gauss_Points_QUAD(4,kesi_N_Enr,yita_N_Enr,weight_N_Enr)
      
      if(Conta_Integ_Point==2)then
          num_CT_Gauss = 2
          kesi_CT(1)  = -0.577350269189626D0
          kesi_CT(2)  =  0.577350269189626D0
          Weight_CT(1)=  ONE
          Weight_CT(2)=  ONE
      elseif(Conta_Integ_Point==1)then
          num_CT_Gauss = 1
          kesi_CT(1)  =  0.0D0
          Weight_CT(1)=  TWO
      endif
      
      EleGaus_yes_FEM_asemd(1:Num_Elem,1:Num_Gauss_Points)= .false.
      do i_E = 1,Num_Elem
          c_thick = thick(Elem_Mat(i_E))
          
          c_D     = D(Elem_Mat(i_E),:,:)  
          
          if(Flag_Weibull_E)then
              if (Key_Weibull_E(Elem_Mat(i_E)) ==1)then
                  c_D = Weibull_Elements_D_Matrix(i_E,1:3,1:3)
              endif
          endif
          
          c_NN    = G_NN(:,i_E)
          c_X_NODES = G_X_NODES(:,i_E)
          c_Y_NODES = G_Y_NODES(:,i_E)   
          Location_ESM(1:200)  = 0; num_Loc_ESM = 0                  
          do i_C =1,num_Crack 
            call Location_Element_Stiff_Matrix(i_E,i_C,
     &                 c_POS(:,i_C),Location_ESM_C_Crack,
     &                 num_Loc_ESM_C_Crack,Location_ESM_C_Cr_NoFEM,
     &                 num_Loc_ESM_C_Cr_NoFEM)
            Location_ESM(num_Loc_ESM+1:num_Loc_ESM+num_Loc_ESM_C_Crack) 
     &                   =  Location_ESM_C_Crack(1:num_Loc_ESM_C_Crack)
            num_Loc_ESM  =  num_Loc_ESM + num_Loc_ESM_C_Crack                     
          end do
          U(1:num_Loc_ESM) = c_DISP(Location_ESM(1:num_Loc_ESM))
          if(Key_Integral_Sol.eq.1)then
          elseif(Key_Integral_Sol.eq.2 .or. Key_Integral_Sol.eq.3)then
              if(maxval(Enriched_Node_Type(c_NN,1:num_Crack)).ne.0)then
                  kesi(1:Num_Gauss_Points)   = kesi_Enr
                  yita(1:Num_Gauss_Points)   = yita_Enr
                  weight(1:Num_Gauss_Points)  = weight_Enr
                  c_Num_Gauss_Point = Num_Gauss_Points
              else 
                  kesi(1:4)   = kesi_N_Enr;  yita(1:4)   = yita_N_Enr
                  weight(1:4)  = weight_N_Enr
                  c_Num_Gauss_Point = 4
              end if             
          end if            
          do i_G = 1, c_Num_Gauss_Point  
            c_kesi = kesi(i_G); c_yita = yita(i_G)
              call Cal_detJ(c_kesi,c_yita,c_X_NODES,c_Y_NODES,detJ)    
              B(1:3,1:80) = ZR
              num_B = 0 
              do i_C =1,num_Crack 
                  call Cal_B_Matrix_Crack(c_kesi,c_yita,i_C,i_E,i_G,
     &                              c_NN,c_X_NODES,c_Y_NODES,
     &                              tem_B,num_tem_B)                  
                  B(1:3,num_B+1:num_B+num_tem_B)=tem_B(1:3,1:num_tem_B)
                  num_B = num_B + num_tem_B                      
              end do
              c_Stress = MATMUL(MATMUL(c_D,B(1:3,1:num_Loc_ESM)),
     &                          U(1:num_Loc_ESM))
              c_Inter_Force(1:num_Loc_ESM) = 
     &               MATMUL(transpose(B(1:3,1:num_Loc_ESM)),c_Stress)
              R_PSI(Location_ESM(1:num_Loc_ESM)) = 
     &           R_PSI(Location_ESM(1:num_Loc_ESM))  + 
     &           c_Inter_Force(1:num_Loc_ESM)*c_thick*weight(i_G)*detJ
          end do
      enddo

      R_PSI =  R_PSI - F_U
      
      do i_C=1,num_Crack
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
     &                 sin(pi/TWO + ori_CT_elem)]
              do i_CT_GAUSS  = 1,num_CT_Gauss
                  c_PC_Gauss(1) =  PC_Gauss_x(i_C,i_CT_Elem,i_CT_GAUSS)
                  c_PC_Gauss(2) =  PC_Gauss_y(i_C,i_CT_Elem,i_CT_GAUSS)
                  call Cal_HF_Coor_by_Kesi(kesi_CT(i_CT_GAUSS),
     &                                     x1,y1,x2,y2,
     &                                     CT_GAUSS_x,CT_GAUSS_y)
                  call Cal_Ele_Num_by_Coors(CT_GAUSS_x,CT_GAUSS_y,
     &                                      CT_GAUSS_Elem)
                  call Cal_KesiYita_by_Coor([CT_GAUSS_x,CT_GAUSS_y],
     &                                CT_GAUSS_Elem,c_Kesi,c_Yita)
                  c_NN    = G_NN(1:4,CT_GAUSS_Elem)     
                  call Cal_N(c_Kesi,c_Yita,N)      
                  tem_n(1) = N(1,1);    tem_n(2) = N(1,3)
                  tem_n(3) = N(1,5);    tem_n(4) = N(1,7)      
                  N_CT(1) = (ONE-kesi_CT(i_CT_GAUSS))/TWO      
                  N_CT(2) = (ONE+kesi_CT(i_CT_GAUSS))/TWO   
                  num_Local_W = 0
                  N_W(1:2,1:60) =ZR
                  do i_N = 1,4
                      if (Enriched_Node_Type(c_NN(i_N),i_C) ==2) then
                          num_Local_W = num_Local_W+1
                          Local_W(2*num_Local_W-1)=
     &                                      2*c_POS(c_NN(i_N),i_C)-1
                          Local_W(2*num_Local_W)  =
     &                                      2*c_POS(c_NN(i_N),i_C)
                          N_W(1,2*num_Local_W-1) = TWO*tem_n(i_N)
                          N_W(2,2*num_Local_W)   = TWO*tem_n(i_N)
                      elseif(Enriched_Node_Type(c_NN(i_N),i_C) ==3)then
                          c_Junc_Ele = Node_Jun_elem(c_NN(i_N),i_C)
                          Coor_Tip =Coors_Junction(c_Junc_Ele,i_C,3:4) 
                          Orient_Tip = atan2(Coor_Tip(2)-CT_GAUSS_y,
     &                                       Coor_Tip(1)-CT_GAUSS_x)
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
     &                                      2*c_POS(c_NN(i_N),i_C)-1
                          Local_W(2*num_Local_W)  =
     &                                      2*c_POS(c_NN(i_N),i_C)
                          N_W(1,2*num_Local_W-1)=c_Sign*ele_Sign*
     &                                           TWO*tem_n(i_N)
                          N_W(2,2*num_Local_W)  =c_Sign*ele_Sign*
     &                                           TWO*tem_n(i_N)
                      elseif(Enriched_Node_Type(c_NN(i_N),i_C) ==1)then
                          Num_Enriched_Node = c_POS(c_NN(i_N),i_C)
                          if (Elem_Type(CT_GAUSS_Elem,i_C) .eq. 1) then
                              ref_elem = CT_GAUSS_Elem
                          else
                              find_element = Node_Elements(c_NN(i_N),:)
                              do i_tem=1,num_Node_Elements(c_NN(i_N))
                                  if (Elem_Type(find_element(i_tem),i_C) 
     &                                                     .eq.1)  then
                                      ref_elem = find_element(i_tem)
                                  end if
                              end do
                          end if
                          Coor_AB(1,:)=
     &                            [Coors_Element_Crack(ref_elem,i_C,1),
     &                             Coors_Element_Crack(ref_elem,i_C,2)]
                          Coor_AB(2,:)=
     &                            [Coors_Element_Crack(ref_elem,i_C,3),
     &                             Coors_Element_Crack(ref_elem,i_C,4)]
                          Coor_Tip = [Coor_AB(2,1),Coor_AB(2,2)]
                          r=sqrt((CT_GAUSS_x-Coor_Tip(1))**2 + 
     &                           (CT_GAUSS_y-Coor_Tip(2))**2)
                          Orient_Tip = atan2(Coor_Tip(2)-CT_GAUSS_y,
     &                                       Coor_Tip(1)-CT_GAUSS_x)
                          num_Local_W = num_Local_W+1
                          Local_W(2*num_Local_W-1)=
     &                             2*(c_POS(c_NN(i_N),i_C)) - 1
                          Local_W(2*num_Local_W)  =
     &                             2*(c_POS(c_NN(i_N),i_C))
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
     &                        ele_Sign*c_Sign*TWO*sqrt(r)*tem_n(i_N)
                          N_W(2,2*num_Local_W)   = 
     &                        ele_Sign*c_Sign*TWO*sqrt(r)*tem_n(i_N)
                      endif
                      if(Enriched_Node_Type(c_NN(i_N),i_C)==2)then
                        do j_C = 1,num_Crack
                          if(Enriched_Node_Type(c_NN(i_N),j_C)==3)then
                              do j_E =1,num_Node_Elements(c_NN(i_N))
                                  c_Adj_Ele=Node_Elements(c_NN(i_N),j_E) 
                                  if(Elem_Type(c_Adj_Ele,j_C)==4)then
                                      c_Junc_Ele = c_Adj_Ele 
                                      exit
                                  endif
                              enddo
                              Coor_Tip =Coors_Junction(
     &                                              c_Junc_Ele,i_C,3:4) 
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
     &                                        2*c_POS(c_NN(i_N),j_C)-1
                              Local_W(2*num_Local_W)  =
     &                                        2*c_POS(c_NN(i_N),j_C)
                              N_W(1,2*num_Local_W-1)= -c_Sign*ele_Sign*
     &                                                ONE*tem_n(i_N)
                              N_W(2,2*num_Local_W)  = -c_Sign*ele_Sign*
     &                                                ONE*tem_n(i_N)
                          endif
                        enddo
                      endif
                  enddo
                  R_PSI(Local_W(1:2*num_Local_W)) = 
     &                   R_PSI(Local_W(1:2*num_Local_W))  + 
     &                   MATMUL(transpose(N_W(1:2,1:2*num_Local_W)),
     &                 c_PC_Gauss)*Weight_CT(i_CT_GAUSS)*detJ   
              enddo
          enddo
      enddo
      
      if(i_NR_P==1)then
          call Vector_Norm2(c_num_freeDOF,
     &                      R_PSI(freeDOF(1:c_num_freeDOF)),
     &                      Norm2_Contact_R_PSI_0)   
      endif

      return 
      end SUBROUTINE Cal_Contact_Resid         
