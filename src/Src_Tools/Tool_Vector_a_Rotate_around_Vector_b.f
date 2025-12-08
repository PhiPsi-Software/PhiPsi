 
      SUBROUTINE Tool_Vector_a_Rotate_around_Vector_b(a,b,theta,output)

      use Global_Float_Type 
      
      implicit none
      real(kind=FT),intent(in)::a(3),b(3),theta
      real(kind=FT),intent(out)::output(3)
      real(kind=FT) tem(3),ax,ay,az,bx,by,bz
      
      tem(1:3) = b(1:3)
      call Vector_Normalize(3,tem(1:3))   
      
      ax = a(1)
      ay = a(2)
      az = a(3)
      bx = tem(1)
      by = tem(2)
      bz = tem(3)
      
      output(1)= ax*cos(theta)+(by*az-bz*ay)*sin(theta) + 
     &           bx*(bx*ax + ay*ay + az*az)*(ONE - cos(theta))
      output(2)= ay*cos(theta)+(bx*az-bz*ax)*sin(theta) + 
     &           by*(bx*ax + ay*ay + az*az)*(ONE  - cos(theta))
      output(3)= az*cos(theta)+(bx*ay-by*ax)*sin(theta) + 
     &           bz*(bx*ax + ay*ay + az*az)*(ONE  - cos(theta))
      
      return
      END SUBROUTINE Tool_Vector_a_Rotate_around_Vector_b