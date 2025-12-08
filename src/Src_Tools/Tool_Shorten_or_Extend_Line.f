 
      subroutine Tool_Shorten_or_Extend_Line(Line_AB,delta_L,
     &                                      Point_String,
     &                                      new_Line_AB,new_Point)

      use Global_Float_Type      
      implicit none
      character*1,intent(in) :: Point_String
      real(kind=FT),intent(in) ::Line_AB(2,2),delta_L
      real(kind=FT),intent(out)::new_Line_AB(2,2),new_Point(2)
      
      real(kind=FT) a_x,a_y,b_x,b_y,theta
      
      a_x = Line_AB(1,1)
      a_y = Line_AB(1,2)
      b_x = Line_AB(2,1)
      b_y = Line_AB(2,2)
      

      theta = atan2(b_y-a_y,b_x-a_x)
      
      
      select case(Point_String)
      case('A')
          new_Line_AB(1,1) = Line_AB(1,1)-delta_L*cos(theta)
          new_Line_AB(1,2) = Line_AB(1,2)-delta_L*sin(theta)
          new_Line_AB(2,:) = Line_AB(2,:)
          new_Point = [new_Line_AB(1,1),new_Line_AB(1,2)]
      case('B')  
          new_Line_AB(2,1) = Line_AB(2,1)+delta_L*cos(theta)
          new_Line_AB(2,2) = Line_AB(2,2)+delta_L*sin(theta)
          new_Line_AB(1,:) = Line_AB(1,:)
          new_Point = [new_Line_AB(2,1),new_Line_AB(2,2)]
      end select
      
      return 
      end SUBROUTINE Tool_Shorten_or_Extend_Line                          
