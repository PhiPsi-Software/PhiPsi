 
SUBROUTINE Matrix_Inverse_3x3(Matrix_A,Matrix_invA)   

use Global_Float_Type
implicit none
real(kind=FT),intent(in)::Matrix_A(3,3)
real(kind=FT),intent(out)::Matrix_invA(3,3)

real(kind=FT) tem

tem =  Matrix_A(1,1)*(Matrix_A(2,2)*Matrix_A(3,3) -Matrix_A(2,3)*Matrix_A(3,2))- &
       Matrix_A(1,2)*(Matrix_A(2,1)*Matrix_A(3,3) -Matrix_A(2,3)*Matrix_A(3,1))+ &
       Matrix_A(1,3)*(Matrix_A(2,1)*Matrix_A(3,2) -Matrix_A(2,2)*Matrix_A(3,1))

Matrix_invA(1,1) = Matrix_A(2,2)*Matrix_A(3,3) -Matrix_A(2,3)*Matrix_A(3,2)
Matrix_invA(2,1) = Matrix_A(2,3)*Matrix_A(3,1) -Matrix_A(2,1)*Matrix_A(3,3)
Matrix_invA(3,1) = Matrix_A(2,1)*Matrix_A(3,2) -Matrix_A(2,2)*Matrix_A(3,1)
Matrix_invA(1,2) = Matrix_A(1,3)*Matrix_A(3,2) -Matrix_A(1,2)*Matrix_A(3,3)
Matrix_invA(2,2) = Matrix_A(1,1)*Matrix_A(3,3) -Matrix_A(1,3)*Matrix_A(3,1)
Matrix_invA(3,2) = Matrix_A(1,2)*Matrix_A(3,1) -Matrix_A(1,1)*Matrix_A(3,2)
Matrix_invA(1,3) = Matrix_A(1,2)*Matrix_A(2,3) -Matrix_A(1,3)*Matrix_A(2,2)
Matrix_invA(2,3) = Matrix_A(1,3)*Matrix_A(2,1) -Matrix_A(1,1)*Matrix_A(2,3)
Matrix_invA(3,3) = Matrix_A(1,1)*Matrix_A(2,2) -Matrix_A(1,2)*Matrix_A(2,1)
Matrix_invA = Matrix_invA/tem

return
END SUBROUTINE Matrix_Inverse_3x3
    


