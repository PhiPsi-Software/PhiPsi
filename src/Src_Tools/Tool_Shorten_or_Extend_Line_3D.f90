!-----------------------------------------------------------
! Brief: Shorten or extend a 3D line segment at endpoint A or B.
!
! Parameters:
!   Input:  Line_AB     - Endpoints of the original segment (2x3)
!   Input:  delta_L     - Signed length change (may be negative)
!   Input:  Point_String- 'A' modifies endpoint A; 'B' modifies B
!   Output: new_Line_AB - Updated segment endpoints (2x3)
!   Output: new_Point   - Coordinates of the moved endpoint
!
! Notes:   Direction cosines of AB drive the displacement of
!   the chosen endpoint along the line direction.
!-----------------------------------------------------------

subroutine Tool_Shorten_or_Extend_Line_3D(Line_AB,delta_L, Point_String, new_Line_AB,new_Point)
      ! Shorten or extend line_AB at point a or b by the increment of offset_L.
      ! delta_L can be negative.
      ! Point_String ='A' or 'B'.
      !                A                       B
      !
      !                ---------------------
      !
      ! Point_String ='A',delta_L < 0:
      !
      !                |<---delta_L--->|
      !                                ------
      !              
      !              
      ! Point_String ='B',delta_L > 0:
      !
      !                                        |<---delta_L--->|
      !                --------------------------------------
      ! Reference: http://wenku.baidu.com/view/2613ac1a227916888486d7cb.html, Vectors in Three-Dimensional
      ! Space 2011
      
      use Global_Float_Type      
      implicit none
      character*1,intent(in) :: Point_String
      real(kind=FT),intent(in) ::Line_AB(2,3),delta_L
      real(kind=FT),intent(out)::new_Line_AB(2,3),new_Point(3)
      real(kind=FT) Vector_AB(3),norm_Vector_AB
      real(kind=FT) theta1,theta2,theta3
      
      Vector_AB = Line_AB(2,1:3) - Line_AB(1,1:3) 
      
norm_Vector_AB = sqrt(Vector_AB(1)**2 + Vector_AB(2)**2 + Vector_AB(3)**2)
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
      end subroutine Tool_Shorten_or_Extend_Line_3D                          
