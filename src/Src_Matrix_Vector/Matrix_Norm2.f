 
      SUBROUTINE Matrix_Norm2(n,Matrix,Norm_2)   

      use Global_Float_Type
      implicit none
      integer,intent(in):: n
      real(kind=FT),intent(in):: Matrix(n,n)
      real(kind=FT),intent(out)::Norm_2
     
      Norm_2  = sqrt(sum(Matrix(1:n,1:n)**2))   
      
      return
      END SUBROUTINE Matrix_Norm2
    


