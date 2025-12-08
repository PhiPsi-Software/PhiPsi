 
      subroutine Cal_Contact_Jacobian(
     &       iter,ifra,Counter_Iter,i_NR_P,
     &       c_Total_FD,c_num_freeDOF,ori_globalK,
     &       freeDOF,Kn,Kn_Gauss,Kt_Gauss,CT_State_Gauss,CT_Jacobian)


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
      real(kind=FT),intent(in)::
     &              ori_globalK(c_Total_FD,c_Total_FD),Kn
      real(kind=FT),intent(inout)::
     &        Kt_Gauss(num_Crack,Max_Num_Cr_CalP-1,2),
     &        Kn_Gauss(num_Crack,Max_Num_Cr_CalP-1,2)
      integer, intent(in)::freeDOF(c_Total_FD)
      integer,intent(in)::
     &              CT_State_Gauss(num_Crack,Max_Num_Cr_CalP-1,2)
      real(kind=FT),intent(out)::
     &              CT_Jacobian(c_Total_FD,c_Total_FD)
      integer i_N,i_C
      integer c_NN(4)
      real(kind=FT) c_kesi,c_yita
      real(kind=FT) detJ
      real(kind=FT) kesi_CT(2),Weight_CT(2)
      real(kind=FT) N_CT(2)
      real(kind=FT) CT_GAUSS_x,CT_GAUSS_y
      integer CT_GAUSS_Elem
      integer num_CT_Gauss,i_CT_Elem
      integer Local_W(60),num_Local_W,Num_Div_Points  
      real(kind=FT) x1,y1,x2,y2,ori_n(2),ori_CT_elem
      real(kind=FT) CT_Length,tem_n(4),N_W(2,60),N(2,8)
      real(kind=FT) D_ep(2,2),D_ep_o(2,2),T(2,2)
      integer i_CT_GAUSS
      integer Num_Enriched_Node
      integer ref_elem,i_tem,find_element(8)
      real(kind=FT) Coor_Tip(2),
     &                 r,Coor_AB(2,2),Orient_Tip,ele_Sign,c_Sign
      real(kind=FT) tem_val1,tem_val2,tem_val3,tem_val4,c_kt,c_kn
      integer j_E,c_Adj_Ele,c_Junc_Ele
      integer j_C
      
      CT_Jacobian(1:c_Total_FD,1:c_Total_FD) = ZR
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
      
      CT_Jacobian = ori_globalK
      
      
      do i_C=1,num_Crack
          Num_Div_Points = Cracks_CalP_Num(i_C)
          do i_CT_Elem = 1,Num_Div_Points - 1
              x1= Cracks_CalP_Coors(i_C,i_CT_Elem,1)
              y1= Cracks_CalP_Coors(i_C,i_CT_Elem,2)
              x2= Cracks_CalP_Coors(i_C,i_CT_Elem+1,1)
              y2= Cracks_CalP_Coors(i_C,i_CT_Elem+1,2)
              ori_CT_elem = atan2(y2-y1,x2-x1)
              CT_Length = Cracks_HF_Ele_L(i_C,i_CT_Elem)
              
              detJ   = CT_Length/TWO
              
              ori_n = [cos(pi/TWO + ori_CT_elem),
     &                 sin(pi/TWO + ori_CT_elem)]
              T(1,1:2)   = [ cos(ori_CT_elem), sin(ori_CT_elem)]   
              T(2,1:2)   = [-sin(ori_CT_elem), cos(ori_CT_elem)]  
              do i_CT_GAUSS  = 1,num_CT_Gauss
                  if(CT_State_Gauss(i_C,i_CT_Elem,i_CT_GAUSS) /= 0)then
                      c_kt = Kt_Gauss(i_C,i_CT_Elem,i_CT_GAUSS)
                      c_kn = Kn_Gauss(i_C,i_CT_Elem,i_CT_GAUSS)
                      
                      D_ep_o(1,1:2) = [c_kt, ZR]
                      D_ep_o(2,1:2) = [ZR, c_kn]
                      
                      
                      call Cal_HF_Coor_by_Kesi(kesi_CT(i_CT_GAUSS),
     &                                         x1,y1,x2,y2,
     &                                         CT_GAUSS_x,CT_GAUSS_y)
                      call Cal_Ele_Num_by_Coors(CT_GAUSS_x,CT_GAUSS_y,
     &                                          CT_GAUSS_Elem)
                      call Cal_KesiYita_by_Coor([CT_GAUSS_x,CT_GAUSS_y],
     &                                CT_GAUSS_Elem,c_Kesi,c_Yita)
                      c_NN    = G_NN(1:4,CT_GAUSS_Elem)     
                      call Cal_N(c_Kesi,c_Yita,N)      
                      tem_n(1) = N(1,1);    tem_n(2) = N(1,3)
                      tem_n(3) = N(1,5);    tem_n(4) = N(1,7)      
                      N_CT(1) = (ONE-kesi_CT(i_CT_GAUSS))/TWO      
                      N_CT(2) = (ONE+kesi_CT(i_CT_GAUSS))/TWO   
                      D_ep = MATMUL(MATMUL(transpose(T),D_ep_o),T)
                      num_Local_W = 0
                      N_W(1:2,1:60) =ZR
                      do i_N = 1,4
                          if(Enriched_Node_Type(c_NN(i_N),i_C) ==2)then
                              num_Local_W = num_Local_W+1
                              Local_W(2*num_Local_W-1)=
     &                                        2*c_POS(c_NN(i_N),i_C)-1
                              Local_W(2*num_Local_W)  =
     &                                        2*c_POS(c_NN(i_N),i_C)
                              N_W(1,2*num_Local_W-1) = TWO*tem_n(i_N)
                              N_W(2,2*num_Local_W)   = TWO*tem_n(i_N)
                          elseif(Enriched_Node_Type(c_NN(i_N),i_C) ==3)
     &                             then

                              c_Junc_Ele = Node_Jun_elem(c_NN(i_N),i_C)
                              Coor_Tip =Coors_Junction(
     &                                  c_Junc_Ele,i_C,3:4) 
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
                          elseif(Enriched_Node_Type(c_NN(i_N),i_C)==1)
     &                                                       then
                              Num_Enriched_Node = c_POS(c_NN(i_N),i_C)
                              if(Elem_Type(CT_GAUSS_Elem,i_C).eq.1)then
                                  ref_elem = CT_GAUSS_Elem
                              else
                                find_element=Node_Elements(c_NN(i_N),:)
                                do i_tem=1,num_Node_Elements(c_NN(i_N))
                                  if (Elem_Type(find_element(i_tem),i_C) 
     &                                                    .eq.1)  then
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
     &                         ele_Sign*c_Sign*TWO*sqrt(r)*tem_n(i_N)
                              N_W(2,2*num_Local_W)   = 
     &                         ele_Sign*c_Sign*TWO*sqrt(r)*tem_n(i_N)
                          endif
                          
                        if(Enriched_Node_Type(c_NN(i_N),i_C)==2)then
                          do j_C = 1,num_Crack
                           if(j_C==i_C) then
                               exit
                           endif
                           if(Enriched_Node_Type(c_NN(i_N),j_C)==3)then
                              c_Junc_Ele = Node_Jun_elem(c_NN(i_N),j_C)
                              Coor_Tip=Coors_Junction(
     &                                             c_Junc_Ele,i_C,3:4) 
                              Orient_Tip =atan2(Coor_Tip(2)-CT_GAUSS_y,
     &                                          Coor_Tip(1)-CT_GAUSS_x)
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
                      CT_Jacobian(Local_W(1:2*num_Local_W),
     &                            Local_W(1:2*num_Local_W))
     &                    = CT_Jacobian(Local_W(1:2*num_Local_W),
     &                                  Local_W(1:2*num_Local_W)) 
     &                       + 
     &                MATMUL(MATMUL(transpose(N_W(1:2,1:2*num_Local_W)),
     &                              D_ep),                            
     &                        N_W(1:2,1:2*num_Local_W))
     &                        * Weight_CT(i_CT_GAUSS)*detJ 
      
                      
                  endif
              enddo
          enddo
      enddo
      
      
      return 
      end SUBROUTINE Cal_Contact_Jacobian         
