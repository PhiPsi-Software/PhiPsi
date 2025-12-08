 
      SUBROUTINE Force_Vector_3D_no_Killed_Ele(c_Total_FD,
     &                                         isub,Lambda,globalF)

      use Global_Float_Type
      use Global_Common
      use Global_Model
      use Global_Elem_Area_Vol
      
      implicit none

      integer,intent(in)::isub,c_Total_FD
      real(kind=FT),intent(in)::Lambda
      real(kind=FT),intent(out)::globalF(c_Total_FD)
      
      real(kind=FT) h,density
      integer i,i_E,i_N,cur_Node,mat_num,c_mat_type,c_node
      
      print *,'    Constructing global force vector...'        
      
      globalF(1:c_Total_FD) = ZR
      
      do i = 1,Num_Foc_x
            cur_Node = int(Foc_x(i,1))
          globalF(3*cur_Node-2) = Lambda*Foc_x(i,2)                   
      end do

      do i = 1,Num_Foc_y
            cur_Node = int(Foc_y(i,1))
            globalF(3*cur_Node-1) =   Lambda*Foc_y(i,2)                  
      end do

      do i = 1,Num_Foc_z
            cur_Node = int(Foc_z(i,1))
            globalF(3*cur_Node) =   Lambda*Foc_z(i,2)                  
      end do     
      
      if(Key_Gravity==1) then
          do i_E = 1,Num_Elem
                mat_num    = Elem_Mat(i_E)
                c_mat_type = Material_Type(mat_num)
              if (c_mat_type ==1) then
                  density = Material_Para(mat_num,3)
              else
            density = Material_Para(mat_num,3); 
              end if
                do i_N = 1,8
                  c_node = G_NN(i_N,i_E)  
                  globalF(3*c_node-2)  = globalF(3*c_node-2)-
     &                          g_X_Y_Z(1)*density*Elem_Vol(i_E)/EIG
              end do              
                do i_N = 1,8
                  c_node = G_NN(i_N,i_E)  
                  globalF(3*c_node-1)  = globalF(3*c_node-1)-
     &                          g_X_Y_Z(2)*density*Elem_Vol(i_E)/EIG
              end do
                do i_N = 1,8
                  c_node = G_NN(i_N,i_E)  
                  globalF(3*c_node)  = globalF(3*c_node)-
     &                          g_X_Y_Z(3)*density*Elem_Vol(i_E)/EIG
              end do
            end do
      end if
      

      
      
  199 continue 
  
  
      RETURN
      END SUBROUTINE Force_Vector_3D_no_Killed_Ele
