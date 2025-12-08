 
      SUBROUTINE Matrix_Get_Diagonal_Elements(n,Matrix,vector_D)   

      use Global_Float_Type       
      implicit none
      integer,intent(in)::n
      real(kind=FT),intent(in)::Matrix(n,n)
      real(kind=FT),intent(out)::vector_D(n)
      INTEGER :: i
      
      DO i = 1, n
          vector_D(i) = Matrix(i,i)
      END DO
      
      return
      END SUBROUTINE Matrix_Get_Diagonal_Elements
    


