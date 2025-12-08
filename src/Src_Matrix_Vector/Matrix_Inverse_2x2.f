 
      SUBROUTINE Matrix_Inverse_2x2(Matrix_A,Matrix_invA)   
      
      use Global_Float_Type  
      implicit none
      real(kind=FT),intent(in)::Matrix_A(2,2)
      real(kind=FT),intent(out)::Matrix_invA(2,2)
      
      real(kind=FT) tem
      
      tem = Matrix_A(1,2)*Matrix_A(2,1) - Matrix_A(1,1)*Matrix_A(2,2) 
      
      Matrix_invA(1,1) = -Matrix_A(2,2) /tem
      Matrix_invA(1,2) =  Matrix_A(1,2)/tem
      Matrix_invA(2,1) =  Matrix_A(2,1)/tem
      Matrix_invA(2,2) = -Matrix_A(1,1)/tem      
      
      return
      END SUBROUTINE Matrix_Inverse_2x2
    


