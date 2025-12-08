 
      subroutine Vector_ZeroOut_Neg_value(Vector,n)   
      use Global_Float_Type     
      implicit none
      integer,intent(in):: n
      real(kind=FT),intent(inout)::Vector(n)
      
      where (Vector<ZR)
          Vector=ZR
      end where
      
      return 
      end subroutine Vector_ZeroOut_Neg_value
    


