 
      subroutine Tool_Point_Given_Point_Normal_and_Distance_3D(
     &                                      Point,
     &                                      Normal_Vector,Distance,
     &                                      New_Point)
      
      use Global_Float_Type      
      implicit none
      real(kind=FT),intent(in) ::Point(3),Normal_Vector(3),Distance
      real(kind=FT),intent(out)::New_Point(3)
      real(kind=FT) theta1,theta2,theta3,norm_Normal_Vector
      
      norm_Normal_Vector = sqrt(Normal_Vector(1)**2 + 
     &                      Normal_Vector(2)**2 + Normal_Vector(3)**2)
      theta1 = acos(Normal_Vector(1)/norm_Normal_Vector)
      theta2 = acos(Normal_Vector(2)/norm_Normal_Vector) 
      theta3 = acos(Normal_Vector(3)/norm_Normal_Vector) 
      
      New_Point(1) = Point(1) + Distance*cos(theta1)
      New_Point(2) = Point(2) + Distance*cos(theta2)
      New_Point(3) = Point(3) + Distance*cos(theta3)
      
      return 
      end SUBROUTINE Tool_Point_Given_Point_Normal_and_Distance_3D                      
