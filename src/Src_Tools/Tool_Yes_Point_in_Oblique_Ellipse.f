 
      subroutine Tool_Yes_Point_in_Oblique_Ellipse(point,
     &                  x0,y0,a,b,theta,
     &                  Return_Statu,Tol)

      use Global_Float_Type
      use Global_Inclusion
      use Global_Common
      
      implicit none
      real(kind=FT),intent(in)::point(2),x0,y0,a,b,theta
      real(kind=FT),intent(in):: Tol
      integer,intent(out):: Return_Statu
      real(kind=FT) x,y,temp1,temp2,c,c_theta
      
      c_theta = theta*pi/180.0D0
      c = sqrt(a**2-b**2)
      x = point(1)
      y = point(2)
      temp1 =(a**2-c**2*cos(c_theta)**2)*(x-x0)**2 +
     &       (a**2-c**2*sin(c_theta)**2)*(y-y0)**2 -
     &       c**2*sin(TWO*c_theta)*(x-x0)*(y-y0)
      temp2 = (a**2)*(b**2)
      
      if((temp1-temp2)>=Tol)then
          Return_Statu =2
      elseif((temp1-temp2)<=-Tol)then
          Return_Statu =1
      else 
          Return_Statu =0
      endif
      
      return 
      end SUBROUTINE Tool_Yes_Point_in_Oblique_Ellipse                          
