 
      subroutine Tool_Volume_Hexahedron(A,B,C,D,E,F,G,H,Volume)
       
      use Global_Float_Type
      IMPLICIT NONE

      real(kind=FT),intent(in)::A(3),B(3),C(3),D(3),
     &                             E(3),F(3),G(3),H(3)
      real(kind=FT),intent(out)::Volume
      real(kind=FT) Vol_A_CDHG,Vol_A_BDGF,Vol_A_EFGH
      
      call Tool_Volume_Pyramid(A,C,D,H,G,Vol_A_CDHG)
      call Tool_Volume_Pyramid(A,B,C,G,F,Vol_A_BDGF)
      call Tool_Volume_Pyramid(A,E,F,G,H,Vol_A_EFGH)
      
      Volume = Vol_A_CDHG + Vol_A_BDGF + Vol_A_EFGH
      
      
      RETURN
      END subroutine Tool_Volume_Hexahedron
