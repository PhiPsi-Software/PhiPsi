 
      subroutine Tool_Point_Giv_2P_and_L(A,B,L,Out_C,Statu,L_BC)




      use Global_Float_Type     
      implicit none
      real(kind=FT),intent(in)::A(2),B(2),L
      real(kind=FT),intent(out)::Out_C(2),L_BC
      integer,intent(out)::Statu
      real(kind=FT) L_AB,theta_AB
     
      L_AB = sqrt((A(1)-B(1))**2+(A(2)-B(2))**2)
      if(L_AB > L)then
          Statu = -1
      elseif(L_AB < L)then
          Statu = 1
      else
          Statu = 0
      endif
      theta_AB = atan2(B(2)-A(2),B(1)-A(1))

      Out_C(1) = A(1)+L*cos(theta_AB)
      Out_C(2) = A(2)+L*sin(theta_AB)
      
      L_BC = sqrt((Out_C(1)-B(1))**2+(Out_C(2)-B(2))**2)
      
      return 
      end SUBROUTINE Tool_Point_Giv_2P_and_L                          
