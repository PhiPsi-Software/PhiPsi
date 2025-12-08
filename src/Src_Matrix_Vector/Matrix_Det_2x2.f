 
      SUBROUTINE Matrix_Det_2x2(Matrix,det)   
      use Global_Float_Type       
      implicit none
      real(kind=FT),intent(in)::Matrix(2,2)
      real(kind=FT),intent(out)::det
      
      det = Matrix(1,1)*Matrix(2,2) - Matrix(1,2)*Matrix(2,1)
        
      return
      END SUBROUTINE Matrix_Det_2x2
    


