 
      subroutine Cal_and_Check_Junction_Point(ifra,iter,Point_B,
     &                                Point_C,Point_D,
     &                                Inter_x,Inter_y)
     
     
      use Global_Float_Type
      use Global_Model 
      use Global_Elem_Area_Vol
      use Global_Common
      implicit none
      integer,intent(in)::ifra,iter
      real(kind=FT),intent(in)::Point_B(2),Point_C(2),Point_D(2)
      real(kind=FT),intent(inout)::Inter_x,Inter_y
      integer N1,N2,N3,N4,NN_L(5),i_Side
      real(kind=FT) point_side_1(2),point_side_2(2),c_Y,c_X
      logical c_Yes_Cross_1,c_Yes_Cross_2,Yes_Accept
      integer num_Cross_Side(4),c_Elem,n_node_1,n_node_2
      real(kind=FT) Point_J_L(2),Point_J_R(2),delta_L
      integer i_Try
      real(kind=FT) line_J_C(2,2),line_J_D(2,2),Point_J(2),incre_L
      real(kind=FT) shorted_Line_AB(2,2)

      Yes_Accept = .False.
      
      Point_J(1) = Inter_x; Point_J(2) = Inter_y
      call Cal_Ele_Num_by_Coors(Inter_x,Inter_y,c_Elem)
      N1  = Elem_Node(c_Elem,1); N2  = Elem_Node(c_Elem,2)                               
      N3  = Elem_Node(c_Elem,3); N4  = Elem_Node(c_Elem,4)                                   
      NN_L= [N1,N2,N3,N4,N1]
      
      num_Cross_Side(1:4) = 0
      do i_Side = 1,4
          n_node_1 = NN_L(i_Side); n_node_2 = NN_L(i_Side+1)
        point_side_1 = [Coor(n_node_1,1),Coor(n_node_1,2)]
        point_side_2 = [Coor(n_node_2,1),Coor(n_node_2,2)]
          call Tool_Intersection(Point_C,Point_D,
     &              point_side_1,point_side_2,c_X,c_Y,c_Yes_Cross_1) 
          
          if (c_Yes_Cross_1.eqv..True.) then
              num_Cross_Side(i_Side)= num_Cross_Side(i_Side)+1
          end if
           call Tool_Intersection(point_side_1,point_side_2,
     &                  Point_B,Point_J,c_X,c_Y,c_Yes_Cross_2) 
            if (c_Yes_Cross_2.eqv..True.) then
                num_Cross_Side(i_Side)= num_Cross_Side(i_Side)+1
            end if
      end do  
      if (maxval(num_Cross_Side) < 2) then
          Yes_Accept =.True.
          goto 99
      endif
      
      delta_L = Ave_Elem_L_Enrich*0.05D0
      line_J_C(1,1:2) =  Point_J; line_J_C(2,1:2) =  Point_C
      line_J_D(1,1:2) =  Point_J; line_J_D(2,1:2) =  Point_D
      do i_Try = 1,50
          print *, '    Junction corection process:',i_Try
          incre_L = -delta_L* i_Try
          call Tool_Shorten_or_Extend_Line(line_J_C,incre_L,'A',
     &                               shorted_Line_AB,Point_J_L)
          call Tool_Shorten_or_Extend_Line(line_J_D,incre_L,'A',
     &                               shorted_Line_AB,Point_J_R)
          call Cal_Ele_Num_by_Coors(Point_J_L(1),Point_J_L(2),c_Elem)
          N1  = Elem_Node(c_Elem,1); N2  = Elem_Node(c_Elem,2)                               
          N3  = Elem_Node(c_Elem,3); N4  = Elem_Node(c_Elem,4)                                   
          NN_L= [N1,N2,N3,N4,N1]
          num_Cross_Side(1:4) = 0
          do i_Side = 1,4
              n_node_1 = NN_L(i_Side); n_node_2 = NN_L(i_Side+1)
            point_side_1 = [Coor(n_node_1,1),Coor(n_node_1,2)]
            point_side_2 = [Coor(n_node_2,1),Coor(n_node_2,2)]
              call Tool_Intersection(Point_C,Point_D,
     &                  point_side_1,point_side_2,
     &                  c_X,c_Y,c_Yes_Cross_1) 
              if (c_Yes_Cross_1.eqv..True.) then
                  num_Cross_Side(i_Side)= num_Cross_Side(i_Side)+1
              end if
              call Tool_Intersection(point_side_1,point_side_2,
     &                  Point_B,Point_J_L,c_X,c_Y,c_Yes_Cross_2) 
              if (c_Yes_Cross_2.eqv..True.) then
                  num_Cross_Side(i_Side)= num_Cross_Side(i_Side)+1
              end if
          end do  
      
          if (maxval(num_Cross_Side) < 2) then
              Inter_x = Point_J_L(1); Inter_y = Point_J_L(2)
              Yes_Accept =.True.
              goto 99
          endif
          call Cal_Ele_Num_by_Coors(Point_J_R(1),Point_J_R(2),c_Elem)
          N1  = Elem_Node(c_Elem,1); N2  = Elem_Node(c_Elem,2)                               
          N3  = Elem_Node(c_Elem,3); N4  = Elem_Node(c_Elem,4)                                   
          NN_L= [N1,N2,N3,N4,N1]
          num_Cross_Side(1:4) = 0
          do i_Side = 1,4
              n_node_1 = NN_L(i_Side); n_node_2 = NN_L(i_Side+1)
            point_side_1 = [Coor(n_node_1,1),Coor(n_node_1,2)]
            point_side_2 = [Coor(n_node_2,1),Coor(n_node_2,2)]
              call Tool_Intersection(Point_C,Point_D,
     &                  point_side_1,point_side_2,
     &                  c_X,c_Y,c_Yes_Cross_1) 
              if (c_Yes_Cross_1.eqv..True.) then
                  num_Cross_Side(i_Side)= num_Cross_Side(i_Side)+1
              end if
              call Tool_Intersection(
     &                  point_side_1,point_side_2,
     &                  Point_B,Point_J_R,
     &                  c_X,c_Y,c_Yes_Cross_2) 
              if (c_Yes_Cross_2.eqv..True.) then
                  num_Cross_Side(i_Side)= num_Cross_Side(i_Side)+1
              end if
          end do  
          if (maxval(num_Cross_Side) < 2) then
              Inter_x = Point_J_R(1); Inter_y = Point_J_R(2)
              Yes_Accept =.True.
              goto 99
          endif      
      enddo
      
   99 continue
      
      if(Yes_Accept .eqv. .False.)then
          print *,'    Error :: Cannot find avaliable Junction point!'
          print *,'             Error in Cal_and_Check_Junction_Point.f'
          call Warning_Message('S',Keywords_Blank) 
      endif
      
      return 
      end SUBROUTINE Cal_and_Check_Junction_Point                  
