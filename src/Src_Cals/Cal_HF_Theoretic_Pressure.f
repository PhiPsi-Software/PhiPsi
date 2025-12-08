 
      subroutine Cal_HF_Theoretic_Pressure(iter,Counter_Iter,
     &            Inject_Crack_L,Theor_Time,Theor_Press,Theor_Aper)
      
      use Global_Float_Type
      use Global_Common
      use Global_Crack
      use Global_Material
      use Global_HF

      implicit none
      integer, intent(in)::iter,Counter_Iter
      real(kind=FT),intent(in)::Inject_Crack_L
      real(kind=FT),intent(out)::Theor_Time,Theor_Press,Theor_Aper
      real(kind=FT)  A0,A1,B1,E_1,u_1
      real(kind=FT) Q_0 
      real(kind=FT) tem,tem1,tem2,tem3,ttt1,ttt2,ttt3,Kesi
      
      A0 = sqrt(THR)
      A1 = -0.156D0
      B1 = 0.0663D0
      E_1 = Material_Para(1,1)/(ONE-Material_Para(1,2)**2)
      u_1 = 12.0D0*Viscosity
      Q_0 = Inject_Q_Val(1)
      
      Theor_Time=((Inject_Crack_L/0.616D0)**SIX*u_1/E_1/(Q_0**THR))
     &           **(ONE/FOU)
      
      Kesi = ZR
      call Tool_Function_BETA(ONE/TWO,TWO/THR,ttt1)
      tem1 = ttt1/(THR*pi)
      ttt2 = ONE
      ttt3 = ONE
      tem2 =  A0*ttt2+TEN/SEV*A1*ttt3
      tem3 =  B1*(TWO-pi*abs(Kesi))
      tem  = tem1*tem2+tem3
      Theor_Press = ((u_1/E_1/Theor_Time)**(ONE/THR))*E_1*tem
      print *,'    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      print *,'    Theoretic pressure(MPa):  ',Theor_Press/1.0D6
      print *,'    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      tem1 = A0+A1
        tem2 = FOU
      tem = tem1 + B1*(tem2+tem3)
      Theor_Aper = ((u_1/E_1/Theor_Time)**(ONE/THR))*
     &              Inject_Crack_L*tem
      
      return 
      end subroutine Cal_HF_Theoretic_Pressure