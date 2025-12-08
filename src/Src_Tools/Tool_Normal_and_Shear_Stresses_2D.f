 
      subroutine Tool_Normal_and_Shear_Stresses_2D(S_xx,S_yy,S_xy,
     &                           theta_stress,S_normal,S_shear)
      use Global_Float_Type
      use Global_Common
      
      implicit none
      real(kind=FT),intent(in)::S_xx,S_yy,S_xy,theta_stress
      real(kind=FT),intent(out)::S_normal,S_shear


      S_normal=(S_xx+S_yy)/TWO + (S_xx-S_yy)/TWO*cos(TWO*theta_stress)-
     &          S_xy*sin(TWO*theta_stress)
      S_shear =(S_xx-S_yy)/TWO*sin(TWO*theta_stress)+
     &          S_xy*cos(TWO*theta_stress)      
      RETURN
      END SUBROUTINE Tool_Normal_and_Shear_Stresses_2D