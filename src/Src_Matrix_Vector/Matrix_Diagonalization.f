 
      SUBROUTINE Matrix_Diagonalization(n,Matrix)   

      use Global_Float_Type       
      implicit none
      integer,intent(in)::n
      real(kind=FT),intent(inout)::Matrix(n,n)
      real(kind=FT) ::Matrix_in(n,n),Matrix_out(n,n)
      
      INTEGER :: i
      
      Matrix_in = Matrix
      Matrix_out(1:n,1:n) = ZR
      
      DO i = 1, n
          Matrix_out(i,i) = sum(Matrix_in(i,1:n))
          
      END DO
      
      Matrix = Matrix_out
      return
      END SUBROUTINE Matrix_Diagonalization
    


