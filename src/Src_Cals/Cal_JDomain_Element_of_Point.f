 
      subroutine Cal_JDomain_Element_of_Point(x,y,r_JDomain,
     &                                        J_Elem,num_J_Elem,
     &                                        Q_Elem_Nodes)
     
      use Global_Float_Type
      use Global_Model
      use Global_Common
      
      implicit none
      
      real(kind=FT),intent(in)::x,y,r_JDomain
      integer,intent(out)::J_Elem(300),num_J_Elem          
      real(kind=FT),intent(out)::Q_Elem_Nodes(num_Elem,4)
      
      integer i_Elem,i_Node
      integer all_J_Elem(300),all_num_J_Elem 
      real(kind=FT) distance,c_Node_x,c_Node_y
      integer i_Try,Total_Try 
      real(kind=FT) c_r_JDomain
      
      
      num_J_Elem = 0
      Total_Try  = 10
      
      do i_Try = 1,Total_Try 
          c_r_JDomain = r_JDomain * dble(i_Try) 
          J_Elem(1:300) = 0
          all_J_Elem(300) = 0
          Q_Elem_Nodes(1:num_Elem,1:4) =ZR
          all_num_J_Elem = 0
          do i_Elem=1,Num_Elem
              do i_Node=1,4
                  c_Node_x = G_X_NODES(i_Node,i_Elem)
                  c_Node_y = G_Y_NODES(i_Node,i_Elem)
                  distance = sqrt((x-c_Node_x)**2 + (y-c_Node_y)**2)
                  if (distance <= c_r_JDomain) then
                      all_num_J_Elem = all_num_J_Elem + 1
                      all_J_Elem(all_num_J_Elem) = i_Elem
                      Q_Elem_Nodes(i_Elem,i_Node) = 1
                  end if
              end do
          end do
          if (all_num_J_Elem>=4) then 
              goto 1000
          else
              if(i_Try>=Total_Try)then
                print *, '    Error :: Cannot find J domain elements!'
                print *, '             Error in Cal_JDomain_Element_of'
     &                            //'_Point.f!'
                call Warning_Message('S',Keywords_Blank) 
              endif
          endif
          
      enddo
      
 1000 continue
 
      call Vector_Unique_Int(300,all_num_J_Elem,all_J_Elem,
     &                             J_Elem,num_J_Elem)   


      return 
      end SUBROUTINE Cal_JDomain_Element_of_Point                
