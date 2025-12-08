 
      subroutine Tool_Yes_On_Line(x,y,A,B,Yes_ON)
      use Global_Float_Type
      use Global_Elem_Area_Vol
      
      implicit none
      real(kind=FT),intent(in):: x,y,A(2),B(2)
      logical,intent(out):: Yes_ON

      real(kind=FT) L_AB,L_AP,L_BP
      
      Yes_ON = .False.

      L_AB = sqrt((A(1)-B(1))**2+(A(2)-B(2))**2)
      L_AP = sqrt((A(1)-x)**2+(A(2)-y)**2)
      L_BP = sqrt((B(1)-x)**2+(B(2)-y)**2)
      
      if ((abs((L_AP + L_BP)-L_AB)) .le. 
     &               1.0D-12*sqrt(Ave_Elem_Area)) then
          Yes_ON = .True. 
      end if
      
      return 
      end SUBROUTINE Tool_Yes_On_Line                          
