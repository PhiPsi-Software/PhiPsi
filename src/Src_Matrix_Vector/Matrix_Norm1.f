 
      SUBROUTINE Matrix_Norm1(n,Matrix,Norm_1)   

      use Global_Float_Type
      implicit none
      integer,intent(in):: n
      real(kind=FT),intent(in):: Matrix(n,n)
      real(kind=FT),intent(out)::Norm_1
     
      Norm_1  = sum(abs(Matrix(1:n,1:n)))   

      
      return
      END SUBROUTINE Matrix_Norm1
    


