 
      subroutine Tool_Angle_of_Vectors_a_and_b_3D(a,b,angle,Key_Method)

      use Global_Float_Type
      use Global_Common
      
      implicit none
      integer,intent(in)::Key_Method
      real(kind=FT),intent(in)::a(3),b(3)
      real(kind=FT),intent(out)::angle
      real(kind=FT) x1,y1,z1,x2,y2,z2
      real(kind=FT) a_plus_b,norm_a,norm_b
      real(kind=FT) Direction(3),tem_angle
      real(kind=FT) tem_check
      real(kind=FT) angle_rad_3d
      real(kind=FT) Cross_Vector_a_b(3),norm_Cross_Vector_a_b
      real(kind=FT) p1(3), p2(3), p3(3)
      real(kind=FT) Direction_Flag(3)
      real(kind=FT) a_Normalized(3),b_Normalized(3)
      
      a_Normalized(1:3) = a
      b_Normalized(1:3) = b
      call Vector_Normalize(3,a_Normalized)   
      call Vector_Normalize(3,b_Normalized)   
      
      x1 = a_Normalized(1)
      y1 = a_Normalized(2)
      z1 = a_Normalized(3)
      
      x2 = b_Normalized(1)
      y2 = b_Normalized(2)
      z2 = b_Normalized(3)
      
      if(abs(x1+x2)< Tol_8  .and.
     &   abs(y1+y2)< Tol_8  .and.
     &   abs(z1+z2)< Tol_8) then
         angle = pi
         return
      endif
      
      if(sum(abs(a-b))<=Tol_10)then
          angle = ZR
          return
      endif
          
      a_plus_b = x1*x2 + y1*y2 + z1*z2
      norm_a   = sqrt(x1**2 + y1**2 + z1**2)
      norm_b   = sqrt(x2**2 + y2**2 + z2**2)
      
      
      if(Key_Method==1) then
          tem_angle    = acos(a_plus_b/norm_a/norm_b)
          
          
          call Vector_Cross_Product_3(a_Normalized,
     &                                b_Normalized,Direction)
          
          
          
          Direction_Flag(1) = 0.52132278048555858D0      
          Direction_Flag(2) = 0.57277733891948324D0
          Direction_Flag(3) = 0.63257306184117679D0
          
          tem_check = Direction(1)*Direction_Flag(1) + 
     &                Direction(2)*Direction_Flag(2) + 
     &                Direction(3)*Direction_Flag(3)
          
          
          if (tem_check >= ZR)then
              angle = tem_angle
          else
              angle = TWO*pi-tem_angle
          endif
      endif



      if(Key_Method==2) then
          
          call Vector_Cross_Product_3(a_Normalized,b_Normalized,
     &                                Cross_Vector_a_b)  
          norm_Cross_Vector_a_b = sqrt(Cross_Vector_a_b(1)**2+
     &                                 Cross_Vector_a_b(2)**2+
     &                                 Cross_Vector_a_b(3)**2)
          angle = atan2(norm_Cross_Vector_a_b,a_plus_b)
      endif
      


           
      if(Key_Method==3) then
          p2 = ZR
          p1 = a_Normalized
          p3 = b_Normalized
          angle = angle_rad_3d(p1, p2, p3)
      endif    
      
      
#ifndef Silverfrost
      if (isnan(angle)) then
          print *, '    Error :: NAN in Tool_Angle_of'
     &                  //'_Vectors_a_and_b_3D.f'
      endif
#endif  
      return 
      
      end SUBROUTINE Tool_Angle_of_Vectors_a_and_b_3D          
