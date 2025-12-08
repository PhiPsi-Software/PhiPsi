 
      subroutine Tool_Yes_Point_on_Line_Segment_3D(A,B,Point,Yes)

      use Global_Float_Type 
      implicit none
      real(kind=FT), intent(in) :: A(3),B(3),Point(3)
      logical, intent(out) :: Yes
      real(kind=FT) Tool_Function_2Point_Dis_3D
      real(kind=FT) DIS_AB,DIS_PA,DIS_PB 
      Yes = .False.
      DIS_AB = Tool_Function_2Point_Dis_3D(A,B)
      DIS_PA = Tool_Function_2Point_Dis_3D(Point,A)
      DIS_PB = Tool_Function_2Point_Dis_3D(Point,B)
      if(abs((DIS_PA + DIS_PB)-DIS_AB) <= Tol_11)then
          Yes  = .True.
      endif
      
      return 
      end SUBROUTINE Tool_Yes_Point_on_Line_Segment_3D                      
