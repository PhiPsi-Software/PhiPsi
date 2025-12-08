 
      subroutine Tool_Centroid_Quad(x,y,Centroid)
       
      
      use Global_Float_Type
      IMPLICIT NONE
      real(kind=FT),intent(in):: x(4)
      real(kind=FT),intent(in):: y(4)
      real(kind=FT),intent(out)::Centroid(2)
      Centroid(1) = sum(x(1:4))/FOU
      Centroid(2) = sum(y(1:4))/FOU
      

      RETURN
      END subroutine Tool_Centroid_Quad
