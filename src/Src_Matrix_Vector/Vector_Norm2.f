 
      SUBROUTINE Vector_Norm2(n,Vector,Norm_2)   
      use Global_Float_Type
      implicit none
      integer,intent(in):: n
      real(kind=FT),intent(in):: Vector(n)
      real(kind=FT),intent(out)::Norm_2
     
      real(kind=FT) tem
     
      tem = sum(Vector(1:n)**2)
      Norm_2  = sqrt(tem)
      
      return
      END SUBROUTINE Vector_Norm2
    


