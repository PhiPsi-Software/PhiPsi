 
      SUBROUTINE Vector_Cross_Product_3(a,b,Out_Vec)   
      use Global_Float_Type          
      implicit none
      real(kind=FT),intent(in)::a(3),b(3)
      real(kind=FT),intent(out)::Out_Vec(3)
      
      Out_Vec(1) = a(2) * b(3) - a(3) * b(2)
      Out_Vec(2) = a(3) * b(1) - a(1) * b(3)
      Out_Vec(3) = a(1) * b(2) - a(2) * b(1)
      
      return
      END SUBROUTINE Vector_Cross_Product_3
    


