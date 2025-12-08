 
      subroutine Tool_Yes_Two_Circles_Intersects(x0,y0,R0,
     &                                           x1,y1,R1,
     &                                           Yes_Intersect)
      use Global_Float_Type
      implicit none
      real(kind=FT),intent(in)::x0,y0,R0,x1,y1,R1
      logical,intent(out)::Yes_Intersect
      real(kind=FT) tem
      
      Yes_Intersect = .False.
      
      tem = SQRT((x0-x1)**2+(y0-y1)**2) 
      
      if (tem <= (R0+R1))then
          Yes_Intersect = .True.
      endif
      
      return 
      end SUBROUTINE Tool_Yes_Two_Circles_Intersects                          
