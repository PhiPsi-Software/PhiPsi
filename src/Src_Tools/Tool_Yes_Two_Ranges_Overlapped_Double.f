 
      subroutine Tool_Yes_Two_Ranges_Overlapped_Double(Range_1,Range_2,
     &                            Logical_Yes)   
     

      use Global_Float_Type
      
      implicit none
      real(kind=FT),intent(in)::Range_1(2),Range_2(2)
      logical,intent(out)::Logical_Yes
      real(kind=FT) max_x1_y1,min_x2_y2
      
      Logical_Yes =.False.
      
      max_x1_y1 = max(Range_1(1),Range_2(1))
      min_x2_y2 = min(Range_1(2),Range_2(2))
      
      if(max_x1_y1 <= min_x2_y2)then
          Logical_Yes =.True.
      endif
      
      if(abs(max_x1_y1-min_x2_y2) <=Tol_11) then
          Logical_Yes =.True.
      endif
      
      return 
      end SUBROUTINE Tool_Yes_Two_Ranges_Overlapped_Double          
