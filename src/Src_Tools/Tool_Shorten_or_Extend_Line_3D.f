 
      subroutine Tool_Shorten_or_Extend_Line_3D(Line_AB,delta_L,
     &                                      Point_String,
     &                                      new_Line_AB,new_Point)
      
      use Global_Float_Type      
      implicit none
      character*1,intent(in) :: Point_String
      real(kind=FT),intent(in) ::Line_AB(2,3),delta_L
      real(kind=FT),intent(out)::new_Line_AB(2,3),new_Point(3)
      real(kind=FT) Vector_AB(3),norm_Vector_AB
      real(kind=FT) theta1,theta2,theta3
      
      Vector_AB = Line_AB(2,1:3) - Line_AB(1,1:3) 
      
      norm_Vector_AB = sqrt(Vector_AB(1)**2 + 
     &                      Vector_AB(2)**2 + Vector_AB(3)**2)
      theta1 = acos(Vector_AB(1)/norm_Vector_AB)
      theta2 = acos(Vector_AB(2)/norm_Vector_AB) 
      theta3 = acos(Vector_AB(3)/norm_Vector_AB) 
      
      select case(Point_String)
      case('A')
          new_Line_AB(1,1) = Line_AB(1,1)-delta_L*cos(theta1)
          new_Line_AB(1,2) = Line_AB(1,2)-delta_L*cos(theta2)
          new_Line_AB(1,3) = Line_AB(1,3)-delta_L*cos(theta3)
          new_Line_AB(2,1:3) = Line_AB(2,1:3)
          new_Point = new_Line_AB(1,1:3)
          
      case('B')  
          new_Line_AB(2,1) = Line_AB(2,1)+delta_L*cos(theta1)
          new_Line_AB(2,2) = Line_AB(2,2)+delta_L*cos(theta2)
          new_Line_AB(2,3) = Line_AB(2,3)+delta_L*cos(theta3)
          new_Line_AB(1,1:3) = Line_AB(1,1:3)
          new_Point = new_Line_AB(2,1:3)
      end select
      
      return 
      end SUBROUTINE Tool_Shorten_or_Extend_Line_3D                          
