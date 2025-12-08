 
      function Tool_Function_2Point_Dis_3D(A,B)
      
      use Global_Float_Type
      implicit none
      real(kind=FT),intent(in)::A(3),B(3)
      real(kind=FT) :: Tool_Function_2Point_Dis_3D
      
      Tool_Function_2Point_Dis_3D = sqrt((A(1)-B(1))**2
     &                                  +(A(2)-B(2))**2
     &                                  +(A(3)-B(3))**2)
      
      end function Tool_Function_2Point_Dis_3D                        
