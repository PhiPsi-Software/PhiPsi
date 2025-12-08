 
      function Tool_Function_2Point_Dis(A,B)
      
      use Global_Float_Type
      implicit none
      real(kind=FT),intent(in)::A(2),B(2)
      real(kind=FT) :: Tool_Function_2Point_Dis
      
      Tool_Function_2Point_Dis = sqrt((A(1)-B(1))**2+(A(2)-B(2))**2)
      
      end function Tool_Function_2Point_Dis                         
