 
      subroutine Cal_Ele_Num_by_Coors_3D_Try(x,y,z,offset_dis,
     &                            x_out,y_out,z_out,OUT_Elem)

      use Global_Float_Type
      use Global_Model
      use Global_INTERFACE_Tool_Yes_Point_in_3D_Hexahedron

      implicit none

      real(kind=FT),intent(in):: x,y,z,offset_dis
      real(kind=FT),intent(out):: x_out,y_out,z_out
      integer,intent(out)::OUT_Elem
      real(kind=FT) c_x_max,c_x_min,c_y_max,c_y_min,c_z_max,c_z_min

      integer c_count,Potent_Elem(500)
      integer i
      integer c_NN(8)
      logical c_Yes_in,c_Yes_on
      integer i_Try
      real(kind=FT) x_new,y_new,z_new
      real(kind=FT) offset_factor_1,offset_factor_2,offset_factor_3
      real(kind=FT) offset_factor_4
      
      OUT_Elem= 0
      x_out = x
      y_out = y
      z_out = z
      
      offset_factor_1  = ONE
      offset_factor_2  = ONE/TWO
      offset_factor_3  = ONE/THR
      offset_factor_4  = ONEP5
          
      do i_Try = 1,80
          select case(i_Try)
          case(1)
              x_new = x + offset_dis*offset_factor_1
              y_new = y
              z_new = z
          case(2)
              x_new = x - offset_dis*offset_factor_1
              y_new = y
              z_new = z   
          case(3)
              x_new = x 
              y_new = y + offset_dis*offset_factor_1
              z_new = z     
          case(4)
              x_new = x 
              y_new = y - offset_dis*offset_factor_1
              z_new = z     
          case(5)
              x_new = x 
              y_new = y 
              z_new = z + offset_dis*offset_factor_1
          case(6)
              x_new = x 
              y_new = y 
              z_new = z - offset_dis*offset_factor_1       
          case(7)
              x_new = x + offset_dis*offset_factor_1
              y_new = y + offset_dis*offset_factor_1
              z_new = z
          case(8)
              x_new = x
              y_new = y + offset_dis*offset_factor_1
              z_new = z + offset_dis*offset_factor_1   
          case(9)
              x_new = x + offset_dis*offset_factor_1
              y_new = y
              z_new = z + offset_dis*offset_factor_1   
          case(10)
              x_new = x - offset_dis*offset_factor_1
              y_new = y - offset_dis*offset_factor_1
              z_new = z     
          case(11)
              x_new = x 
              y_new = y - offset_dis*offset_factor_1
              z_new = z - offset_dis*offset_factor_1  
          case(12)
              x_new = x - offset_dis*offset_factor_1
              y_new = y 
              z_new = z - offset_dis*offset_factor_1        
          case(13)
              x_new = x + offset_dis*offset_factor_1
              y_new = y + offset_dis*offset_factor_1
              z_new = z + offset_dis*offset_factor_1   
          case(14)
              x_new = x - offset_dis*offset_factor_1
              y_new = y - offset_dis*offset_factor_1
              z_new = z - offset_dis*offset_factor_1      
          case(15)
              x_new = x + offset_dis*offset_factor_1
              y_new = y - offset_dis*offset_factor_1
              z_new = z         
          case(16)
              x_new = x - offset_dis*offset_factor_1
              y_new = y + offset_dis*offset_factor_1
              z_new = z   
          case(17)
              x_new = x + offset_dis*offset_factor_1
              y_new = y 
              z_new = z - offset_dis*offset_factor_1        
          case(18)
              x_new = x - offset_dis*offset_factor_1
              y_new = y 
              z_new = z + offset_dis*offset_factor_1    
          case(19)
              x_new = x 
              y_new = y + offset_dis*offset_factor_1
              z_new = z - offset_dis*offset_factor_1    
          case(20)
              x_new = x 
              y_new = y - offset_dis*offset_factor_1
              z_new = z + offset_dis*offset_factor_1
          case(21)
              x_new = x + offset_dis*offset_factor_2
              y_new = y
              z_new = z
          case(22)
              x_new = x - offset_dis*offset_factor_2
              y_new = y
              z_new = z   
          case(23)
              x_new = x 
              y_new = y + offset_dis*offset_factor_2
              z_new = z     
          case(24)
              x_new = x 
              y_new = y - offset_dis*offset_factor_2
              z_new = z     
          case(25)
              x_new = x 
              y_new = y 
              z_new = z + offset_dis*offset_factor_2
          case(26)
              x_new = x 
              y_new = y 
              z_new = z - offset_dis*offset_factor_2       
          case(27)
              x_new = x + offset_dis*offset_factor_2
              y_new = y + offset_dis*offset_factor_2
              z_new = z
          case(28)
              x_new = x
              y_new = y + offset_dis*offset_factor_2
              z_new = z + offset_dis*offset_factor_2   
          case(29)
              x_new = x + offset_dis*offset_factor_2
              y_new = y
              z_new = z + offset_dis*offset_factor_2   
          case(30)
              x_new = x - offset_dis*offset_factor_2
              y_new = y - offset_dis*offset_factor_2
              z_new = z     
          case(31)
              x_new = x 
              y_new = y - offset_dis*offset_factor_2
              z_new = z - offset_dis*offset_factor_2  
          case(32)
              x_new = x - offset_dis*offset_factor_2
              y_new = y 
              z_new = z - offset_dis*offset_factor_2        
          case(33)
              x_new = x + offset_dis*offset_factor_2
              y_new = y + offset_dis*offset_factor_2
              z_new = z + offset_dis*offset_factor_2   
          case(34)
              x_new = x - offset_dis*offset_factor_2
              y_new = y - offset_dis*offset_factor_2
              z_new = z - offset_dis*offset_factor_2      
          case(35)
              x_new = x + offset_dis*offset_factor_2
              y_new = y - offset_dis*offset_factor_2
              z_new = z         
          case(36)
              x_new = x - offset_dis*offset_factor_2
              y_new = y + offset_dis*offset_factor_2
              z_new = z   
          case(37)
              x_new = x + offset_dis*offset_factor_2
              y_new = y 
              z_new = z - offset_dis*offset_factor_2        
          case(38)
              x_new = x - offset_dis*offset_factor_2
              y_new = y 
              z_new = z + offset_dis*offset_factor_2    
          case(39)
              x_new = x 
              y_new = y + offset_dis*offset_factor_2
              z_new = z - offset_dis*offset_factor_2    
          case(40)
              x_new = x 
              y_new = y - offset_dis*offset_factor_2
              z_new = z + offset_dis*offset_factor_2   
          case(41)
              x_new = x + offset_dis*offset_factor_3
              y_new = y
              z_new = z
          case(42)
              x_new = x - offset_dis*offset_factor_3
              y_new = y
              z_new = z   
          case(43)
              x_new = x 
              y_new = y + offset_dis*offset_factor_3
              z_new = z     
          case(44)
              x_new = x 
              y_new = y - offset_dis*offset_factor_3
              z_new = z     
          case(45)
              x_new = x 
              y_new = y 
              z_new = z + offset_dis*offset_factor_3
          case(46)
              x_new = x 
              y_new = y 
              z_new = z - offset_dis*offset_factor_3       
          case(47)
              x_new = x + offset_dis*offset_factor_3
              y_new = y + offset_dis*offset_factor_3
              z_new = z
          case(48)
              x_new = x
              y_new = y + offset_dis*offset_factor_3
              z_new = z + offset_dis*offset_factor_3   
          case(49)
              x_new = x + offset_dis*offset_factor_3
              y_new = y
              z_new = z + offset_dis*offset_factor_3   
          case(50)
              x_new = x - offset_dis*offset_factor_3
              y_new = y - offset_dis*offset_factor_3
              z_new = z     
          case(51)
              x_new = x 
              y_new = y - offset_dis*offset_factor_3
              z_new = z - offset_dis*offset_factor_3  
          case(52)
              x_new = x - offset_dis*offset_factor_3
              y_new = y 
              z_new = z - offset_dis*offset_factor_3        
          case(53)
              x_new = x + offset_dis*offset_factor_3
              y_new = y + offset_dis*offset_factor_3
              z_new = z + offset_dis*offset_factor_3   
          case(54)
              x_new = x - offset_dis*offset_factor_3
              y_new = y - offset_dis*offset_factor_3
              z_new = z - offset_dis*offset_factor_3      
          case(55)
              x_new = x + offset_dis*offset_factor_3
              y_new = y - offset_dis*offset_factor_3
              z_new = z         
          case(56)
              x_new = x - offset_dis*offset_factor_3
              y_new = y + offset_dis*offset_factor_3
              z_new = z   
          case(57)
              x_new = x + offset_dis*offset_factor_3
              y_new = y 
              z_new = z - offset_dis*offset_factor_3        
          case(58)
              x_new = x - offset_dis*offset_factor_3
              y_new = y 
              z_new = z + offset_dis*offset_factor_3    
          case(59)
              x_new = x 
              y_new = y + offset_dis*offset_factor_3
              z_new = z - offset_dis*offset_factor_3    
          case(60)
              x_new = x 
              y_new = y - offset_dis*offset_factor_3
              z_new = z + offset_dis*offset_factor_3           
          case(61)
              x_new = x + offset_dis*offset_factor_4
              y_new = y
              z_new = z
          case(62)
              x_new = x - offset_dis*offset_factor_4
              y_new = y
              z_new = z   
          case(63)
              x_new = x 
              y_new = y + offset_dis*offset_factor_4
              z_new = z     
          case(64)
              x_new = x 
              y_new = y - offset_dis*offset_factor_4
              z_new = z     
          case(65)
              x_new = x 
              y_new = y 
              z_new = z + offset_dis*offset_factor_4
          case(66)
              x_new = x 
              y_new = y 
              z_new = z - offset_dis*offset_factor_4       
          case(67)
              x_new = x + offset_dis*offset_factor_4
              y_new = y + offset_dis*offset_factor_4
              z_new = z
          case(68)
              x_new = x
              y_new = y + offset_dis*offset_factor_4
              z_new = z + offset_dis*offset_factor_4   
          case(69)
              x_new = x + offset_dis*offset_factor_4
              y_new = y
              z_new = z + offset_dis*offset_factor_4   
          case(70)
              x_new = x - offset_dis*offset_factor_4
              y_new = y - offset_dis*offset_factor_4
              z_new = z     
          case(71)
              x_new = x 
              y_new = y - offset_dis*offset_factor_4
              z_new = z - offset_dis*offset_factor_4  
          case(72)
              x_new = x - offset_dis*offset_factor_4
              y_new = y 
              z_new = z - offset_dis*offset_factor_4        
          case(73)
              x_new = x + offset_dis*offset_factor_4
              y_new = y + offset_dis*offset_factor_4
              z_new = z + offset_dis*offset_factor_4   
          case(74)
              x_new = x - offset_dis*offset_factor_4
              y_new = y - offset_dis*offset_factor_4
              z_new = z - offset_dis*offset_factor_4     
          case(75)
              x_new = x + offset_dis*offset_factor_4
              y_new = y - offset_dis*offset_factor_4
              z_new = z         
          case(76)
              x_new = x - offset_dis*offset_factor_4
              y_new = y + offset_dis*offset_factor_4
              z_new = z   
          case(77)
              x_new = x + offset_dis*offset_factor_4
              y_new = y 
              z_new = z - offset_dis*offset_factor_4        
          case(78)
              x_new = x - offset_dis*offset_factor_4
              y_new = y 
              z_new = z + offset_dis*offset_factor_4    
          case(79)
              x_new = x 
              y_new = y + offset_dis*offset_factor_4
              z_new = z - offset_dis*offset_factor_4    
          case(80)
              x_new = x 
              y_new = y - offset_dis*offset_factor_4
              z_new = z + offset_dis*offset_factor_4                  
          end select
          
          
          c_count = 0
          Potent_Elem  = 0

          if (x_new > Max_X_Coor .or. x_new < Min_X_Coor)then
              continue
          endif
          if (y_new > Max_Y_Coor .or. y_new < Min_Y_Coor)then
              continue
          endif
          if (z_new > Max_Z_Coor .or. z_new < Min_Z_Coor)then
              continue
          endif

          do i=1,Num_Elem
              c_x_max = x_max_Elements(i)
              c_x_min = x_min_Elements(i)
              c_y_max = y_max_Elements(i)
              c_y_min = y_min_Elements(i)
              c_z_max = z_max_Elements(i)
              c_z_min = z_min_Elements(i)
              if   ((x_new.le.c_x_max) .and. (x_new.ge.c_x_min).
     &         and. (y_new.le.c_y_max) .and. (y_new.ge.c_y_min).
     &         and. (z_new.le.c_z_max) .and. (z_new.ge.c_z_min))then
                  c_count = c_count +1
                  Potent_Elem(c_count) = i
              end if
          end do

          do i =1,c_count
              c_NN  = G_NN(1:8,Potent_Elem(i))
              call Tool_Yes_Point_in_3D_Hexahedron([x_new,y_new,z_new],
     &                   Coor(c_NN(1),1:3),Coor(c_NN(2),1:3),
     &                   Coor(c_NN(3),1:3),Coor(c_NN(4),1:3),
     &                   Coor(c_NN(5),1:3),Coor(c_NN(6),1:3),
     &                   Coor(c_NN(7),1:3),Coor(c_NN(8),1:3),
     &                   c_Yes_in,c_Yes_on)
             if(c_Yes_in)then
                 OUT_Elem = Potent_Elem(i)
                 x_out = x_new
                 y_out = y_new
                 z_out = z_new
                 return
             endif
          end do
      enddo
      
      return
      end SUBROUTINE Cal_Ele_Num_by_Coors_3D_Try
