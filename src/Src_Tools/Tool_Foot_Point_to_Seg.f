 
      subroutine Tool_Foot_Point_to_Seg(a,b,c,Foot)

      use Global_Float_Type 
      implicit none
      real(kind=FT), intent(in) :: a(2), b(2), c(2)
      real(kind=FT), intent(out) :: Foot(2)
      real(kind=FT):: x1,y1,x2,y2,x3,y3,x,y
      real(kind=FT):: px,py,dAB,u
      
      x1=a(1)
      y1=a(2)
      x2=b(1)
      y2=b(2)
      x3=c(1)
      y3=c(2)
      px = x2-x1
      py = y2-y1
      dAB = px*px + py*py
      u = ((x3 - x1) * px + (y3 - y1) * py) / dAB
      x = x1 + u * px
      y = y1 + u * py
      Foot(1:2) = [x,y]
      
      return 
      end SUBROUTINE Tool_Foot_Point_to_Seg                       
