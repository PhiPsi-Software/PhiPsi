 
      subroutine Cal_HF_HFNode_Location(A,B,
     &                            Div_Points,Num_Div_Points)
      
      use Global_Float_Type
      use Global_Crack
      use Global_Model
      use Global_Filename
      use Global_Common
      use Global_Elem_Area_Vol
      
      implicit none
      
      real(kind=FT),intent(in)::A(2),B(2)
      real(kind=FT),intent(out)::Div_Points(Max_Num_Seg_CalP,2)
      integer,intent(out)::Num_Div_Points
      
      integer i,i_Side,cc_Location
      real(kind=FT) a_x,a_y,b_x,b_y,delta_L,c_L,c_x_r,c_y_r
      real(kind=FT) theta,c_x_l,c_y_l
      integer Ele_l,Ele_r,num_Inter,N1,N2,N3,N4,NN_L(5),n_node_1,
     &        n_node_2
      real(kind=FT) c_Intersec(5,2),point_side_1(2),point_side_2(2)
      logical c_Yes_Cross,Yes_Already_has
      real(kind=FT) c_Inter_x,c_Inter_y,L_Seg,tem1,tem2
      integer number_Check
      
      
      a_x = A(1)
      a_y = A(2)
      b_x = B(1)
      b_y = B(2)
      
      L_Seg = sqrt(( b_x- a_x)**2+( b_y- a_y)**2) 
     
      delta_L = 0.009999D0*Ave_Elem_L_Enrich
      
      Num_Div_Points = 0
      
      Div_Points(1,:)= A
      Num_Div_Points = 1
      
      
      
      tem1 = b_x-a_x 
      tem2 = b_y-a_y
      
      if (b_x ==a_x) then
          tem1 = Tol_20
      endif
      if (b_y ==a_y) then
          tem2 = Tol_20
      endif
      
      theta = atan2(tem2,tem1)

      
      number_Check = int(L_Seg/delta_L*1.1)
      do i = 1,number_Check
          c_Intersec(1:5,1:2) = ZR
          c_x_l = a_x + (i-1)*delta_L*cos(theta)
          c_y_l = a_y + (i-1)*delta_L*sin(theta)
          c_x_r = a_x + i*delta_L*cos(theta)
          c_y_r = a_y + i*delta_L*sin(theta)
          c_L = i*delta_L
          call Cal_Ele_Num_by_Coors(c_x_l,c_y_l,Ele_l)
          call Cal_Ele_Num_by_Coors(c_x_r,c_y_r,Ele_r)
          c_Intersec(1:5,1:2) = ZR
          if(Ele_l /= Ele_r .and. c_L<=L_Seg)then
              num_Inter = 0
              N1  = Elem_Node(Ele_l,1)                            
              N2  = Elem_Node(Ele_l,2)                             
              N3  = Elem_Node(Ele_l,3)     
              N4  = Elem_Node(Ele_l,4)      
              NN_L= [N1,N2,N3,N4,N1]
              do i_Side = 1,4
                  n_node_1 = NN_L(i_Side)
                  n_node_2 = NN_L(i_Side+1)
                  point_side_1 = [Coor(n_node_1,1),Coor(n_node_1,2)]
                  point_side_2 = [Coor(n_node_2,1),Coor(n_node_2,2)]
                  call Tool_Intersection(
     &                 [c_x_l,c_y_l],[c_x_r,c_y_r],
     &                 point_side_1,point_side_2,
     &                 c_Inter_x,c_Inter_y,c_Yes_Cross) 
                  if (c_Yes_Cross.eqv..True.) then
                      call Vector_belongs_Matrix_Is_Dou(5,2,
     &                                c_Intersec,[c_Inter_x,c_Inter_y],
     &                                cc_Location,Yes_Already_has)   
                      if(Yes_Already_has.eqv..False.) then
                        num_Inter=num_Inter+1
                        c_Intersec(num_Inter,1:2)=[c_Inter_x,c_Inter_y]
                      endif
                  end if
              end do  
              N1  = Elem_Node(Ele_r,1)                                         
              N2  = Elem_Node(Ele_r,2)                                             
              N3  = Elem_Node(Ele_r,3)                                             
              N4  = Elem_Node(Ele_r,4)                                            
              NN_L= [N1,N2,N3,N4,N1]
              do i_Side = 1,4
                  n_node_1 = NN_L(i_Side)
                  n_node_2 = NN_L(i_Side+1)
                  point_side_1 = [Coor(n_node_1,1),Coor(n_node_1,2)]
                  point_side_2 = [Coor(n_node_2,1),Coor(n_node_2,2)]
                  call Tool_Intersection(
     &                 [c_x_l,c_y_l],[c_x_r,c_y_r],
     &                 point_side_1,point_side_2,
     &                 c_Inter_x,c_Inter_y,c_Yes_Cross) 
                  if (c_Yes_Cross.eqv..True.) then
                      call Vector_belongs_Matrix_Is_Dou(5,2,
     &                                c_Intersec,[c_Inter_x,c_Inter_y],
     &                                cc_Location,Yes_Already_has)   
                      if(Yes_Already_has.eqv..False.) then
                        num_Inter=num_Inter+1
                        c_Intersec(num_Inter,1:2)=[c_Inter_x,c_Inter_y]
                      endif
                  end if
              end do
              if(num_Inter/=0)then
                  Div_Points(Num_Div_Points+1:
     &                       Num_Div_Points+num_Inter,:) = 
     &                                   c_Intersec(1:num_Inter,:)
                  Num_Div_Points = Num_Div_Points +num_Inter
              end if
          end if
          
          if(i*delta_L>=L_Seg)then
              exit
          end if
      end do 
              
      Num_Div_Points  = Num_Div_Points + 1
      Div_Points(Num_Div_Points,:) = B
      
      return 
      end SUBROUTINE Cal_HF_HFNode_Location
