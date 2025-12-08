 
      subroutine Cal_HF_Matrix_Q_Linear(Counter,Coupled_Q,c_Total_FD)
      
      use Global_Float_Type
      use Global_Common
      use Global_Model
      use Global_Elem_Area_Vol
      use Global_Crack
      use Global_Crack_Common
      use Global_HF

      implicit none
      integer,intent(in)::Counter,c_Total_FD
      real(kind=FT),intent(out)::
     &                      Coupled_Q(c_Total_FD,num_Tol_CalP_Water)
      integer i_N,i_C,i_HF_Elem,i_HF_GAUSS    
      real(kind=FT) Kesi,Yita,Orient1,Orient2,x1,y1,x2,y2
      integer c_NN(4),Num_Div_Points 
      real(kind=FT) Div_Points(Max_Num_Cr_CalP,2),N(2,8),
     &                 tem_n(4),Coor_Tip(2),
     &                 r,Coor_AB(2,2),Orient_Tip
      integer ref_elem,i_tem,find_element(8)
      
      real(kind=FT) kesi_HF(4),Weight_HF(4)
      real(kind=FT) N_HF(2)
      real(kind=FT) HF_GAUSS_x,HF_GAUSS_y
      integer HF_GAUSS_Elem
      real(kind=FT) HF_Length
      integer Local_P(2),Local_W(60),num_Local_W                           
      real(kind=FT) N_W(2,60),detJ,ori_n(2)
      real(kind=FT) tem_Q(60,2),tem(60,2),c_Sign,ori_HF_elem,
     &                 ele_Sign,tem_val1,tem_val2,tem_val3,tem_val4
      integer num_HF_Gauss,j_C,c_Crac
      integer Num_Enriched_Node,j_E,c_Adj_Ele,c_Junc_Ele
      real(kind=FT) c_X_NODES(4),c_Y_NODES(4),detJ_N_W
      print *,'    Calculating coupled matrix Q......'
      
      num_HF_Gauss = 2
      kesi_HF(1)  = -0.577350269189626D0
      kesi_HF(2)  =  0.577350269189626D0
      Weight_HF(1)=  ONE
      Weight_HF(2)=  ONE
      
 

      Coupled_Q(1:c_Total_FD,1:num_Tol_CalP_Water) = ZR

      Local_P(1:2)=0
      do i_C=1,num_Crack
        if (Cracks_HF_State(i_C) == 1) then  
          Num_Div_Points = Cracks_CalP_Num(i_C)
          Div_Points(1:Num_Div_Points,1:2) = 
     &                     Cracks_CalP_Coors(i_C,1:Num_Div_Points,1:2)
          if(i_C>1)then
              Local_P(1) =  Local_P(1)+1       
          endif
          do i_HF_Elem = 1,Num_Div_Points - 1
              tem(1:60,1:2)  = ZR
              tem_Q(1:60,1:2)= ZR   
              Local_P(1) =  Local_P(1)+1        
              Local_P(2) =  Local_P(1)+1      
            x1= Cracks_CalP_Coors(i_C,i_HF_Elem,1)
            y1= Cracks_CalP_Coors(i_C,i_HF_Elem,2)
            x2= Cracks_CalP_Coors(i_C,i_HF_Elem+1,1)
            y2= Cracks_CalP_Coors(i_C,i_HF_Elem+1,2)
              ori_HF_elem = atan2(y2-y1,x2-x1)
              HF_Length = Cracks_HF_Ele_L(i_C,i_HF_Elem)
              detJ   = HF_Length/TWO
              ori_n = [cos(pi/TWO+ori_HF_elem),
     &                 sin(pi/TWO+ori_HF_elem)]
              do i_HF_GAUSS  = 1,num_HF_Gauss
                  call Cal_HF_Coor_by_Kesi(kesi_HF(i_HF_GAUSS),
     &                                     x1,y1,x2,y2,
     &                                     HF_GAUSS_x,HF_GAUSS_y)
                  call Cal_Ele_Num_by_Coors(HF_GAUSS_x,HF_GAUSS_y,
     &                                      HF_GAUSS_Elem)
                  call Cal_KesiYita_by_Coor([HF_GAUSS_x,HF_GAUSS_y],
     &                                       HF_GAUSS_Elem,Kesi,Yita)
                  c_NN    = G_NN(1:4,HF_GAUSS_Elem)     
                  call Cal_N(Kesi,Yita,N)      
                  c_X_NODES = G_X_NODES(:,HF_GAUSS_Elem)
                  c_Y_NODES = G_Y_NODES(:,HF_GAUSS_Elem)    
                  call Cal_detJ(Kesi,Yita,c_X_NODES,c_Y_NODES,detJ_N_W)                  
                  tem_n(1) = N(1,1);    tem_n(2) = N(1,3)
                  tem_n(3) = N(1,5);    tem_n(4) = N(1,7)      
                  N_HF(1) = (ONE-kesi_HF(i_HF_GAUSS))/TWO      
                  N_HF(2) = (ONE+kesi_HF(i_HF_GAUSS))/TWO   
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
                          Orient_Tip = atan2(Coor_Tip(2)-HF_GAUSS_y,
     &                                       Coor_Tip(1)-HF_GAUSS_x)
                          ele_Sign =  sign(ONE,cos(ori_HF_elem))
                          c_Sign   =  sign(ONE,cos(Orient_Tip))
                          tem_val1 = abs(ori_HF_elem - pi/TWO)
                          tem_val2 = abs(ori_HF_elem + pi/TWO)
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
                          if (Elem_Type(HF_GAUSS_Elem,i_C) .eq. 1) then
                              ref_elem = HF_GAUSS_Elem
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
     &                           [Coors_Element_Crack(ref_elem,i_C,1),
     &                            Coors_Element_Crack(ref_elem,i_C,2)]
                          Coor_AB(2,:)=
     &                           [Coors_Element_Crack(ref_elem,i_C,3),
     &                            Coors_Element_Crack(ref_elem,i_C,4)]
                          Coor_Tip = [Coor_AB(2,1),Coor_AB(2,2)]
                          r=sqrt((HF_GAUSS_x-Coor_Tip(1))**2 + 
     &                           (HF_GAUSS_y-Coor_Tip(2))**2)
                          Orient_Tip = atan2(Coor_Tip(2)-HF_GAUSS_y,
     &                                       Coor_Tip(1)-HF_GAUSS_x)
                          num_Local_W = num_Local_W+1
                          Local_W(2*num_Local_W-1)=
     &                             2*(c_POS(c_NN(i_N),i_C)) - 1
                          Local_W(2*num_Local_W)  =
     &                             2*(c_POS(c_NN(i_N),i_C))
                          ele_Sign =  sign(ONE,cos(ori_HF_elem))
                          c_Sign   =  sign(ONE,cos(Orient_Tip))
                          tem_val1 = abs(ori_HF_elem - pi/TWO)
                          tem_val2 = abs(ori_HF_elem + pi/TWO)
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
                      end if
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
     &                                       c_Junc_Ele,i_C,3:4) 
                              Orient_Tip = atan2(Coor_Tip(2)-HF_GAUSS_y,
     &                                           Coor_Tip(1)-HF_GAUSS_x)
                              ele_Sign =  sign(ONE,cos(ori_HF_elem))
                              c_Sign   =  sign(ONE,cos(Orient_Tip))
                              tem_val1 = abs(ori_HF_elem - pi/TWO)
                              tem_val2 = abs(ori_HF_elem + pi/TWO)
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
                  end do 
                  call Vectors_Multi
     &              (MATMUL(transpose(N_W(1:2,1:2*num_Local_W)),ori_n),
     &                  2*num_Local_W,N_HF,2,tem(1:2*num_Local_W,1:2))
     
                  tem_Q(1:2*num_Local_W,1:2) = 
     &                          tem_Q(1:2*num_Local_W,1:2)  +
     &                          tem(1:2*num_Local_W,1:2)*
     &                          Weight_HF(i_HF_GAUSS)*detJ  
              end do
              Coupled_Q(Local_W(1:2*num_Local_W),Local_P(1:2)) = 
     &             Coupled_Q(Local_W(1:2*num_Local_W),Local_P(1:2)) +
     &             tem_Q(1:2*num_Local_W,1:2)

          end do  
        end if      
      end do
      
      return 
      end SUBROUTINE Cal_HF_Matrix_Q_Linear           
