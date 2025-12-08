 
      subroutine Tool_von_Mises(stress_x,stress_y,stress_xy,
     &                          c_v,stress_VM)

      use Global_Float_Type
      use Global_Common

      implicit none
      real(kind=FT),intent(in)::stress_x,stress_y,stress_xy,c_v
      real(kind=FT),intent(out)::stress_VM
      real(kind=FT) Stress_z,I_1,I_2

      if(Key_Type_2D==1)then
          Stress_z = ZR
      elseif(Key_Type_2D==2)then
          Stress_z = c_v*(stress_x + stress_y)
      endif
      
      stress_VM = 
     &      sqrt(((stress_x-stress_y)**2 + 
     &            (stress_y-stress_z)**2 +
     &            (stress_x-stress_z)**2 +
     &        SIX*(stress_xy**2))/TWO)      

      return
      end SUBROUTINE Tool_von_Mises
