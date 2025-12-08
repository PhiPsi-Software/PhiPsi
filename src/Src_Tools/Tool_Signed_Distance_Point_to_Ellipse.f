 
      subroutine Tool_Signed_Distance_Point_to_Ellipse(
     &                      x0,y0,a,b,theta,Point_C,S_Distance)

      use Global_Float_Type      
      implicit none
      real(kind=FT),intent(in)::x0,y0,a,b,theta,Point_C(2)
      real(kind=FT),intent(out)::S_Distance
      
      real(kind=FT) c_Dis,Tol
      integer Return_Statu
      
      Tol = 1.0D-10
      
      call Tool_Yes_Point_in_Oblique_Ellipse(Point_C,
     &                  x0,y0,a,b,theta,
     &                  Return_Statu,Tol)
      if(Return_Statu==1)then
          S_Distance =-ONE
      elseif(Return_Statu==2)then
          S_Distance =ONE
      else
          S_Distance =ZR
      endif
      return 
      end SUBROUTINE Tool_Signed_Distance_Point_to_Ellipse                        
