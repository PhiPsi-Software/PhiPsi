!-----------------------------------------------------------
! Brief: Compute the 1- or 2-norm condition number of a dense matrix.
!
! Parameters:
!   Input:  norm_Type         - which norm to use (1 or 2)
!           n                 - order of the square matrix
!           Matrix            - the n by n matrix
!   Output: Condition_Num_Norm - the resulting condition number
!
! Notes:   Uses an explicit inverse to obtain Norm(K)*Norm(K^-1).
!          If the matrix is singular the output may be NaN.
!-----------------------------------------------------------

subroutine Matrix_Condition_Number_v2(norm_Type,n,Matrix, Condition_Num_Norm)
!     Calculate the condition number of a matrix. Suitable for calculating the condition number of smaller matrices.
!     Fist written on 2022-07-10 by Fang Shi.
!     Test passed.
!     Note: If the matrix does not have an inverse, the output may be NAN.

use Global_Float_Type

implicit none
integer,intent(in)::n,norm_Type
real(kind=FT),intent(in)::Matrix(n,n)
real(kind=FT),intent(out)::Condition_Num_Norm

real(kind=FT) Norm_2_K,Norm_2_K_inv
real(kind=FT) Norm_1_K,Norm_1_K_inv
real(kind=FT) Matrix_inv(n,n)


call Matrix_Inverse(Matrix,Matrix_inv,n)   

if (norm_Type==2) then
    call Matrix_Norm2(n,Matrix,Norm_2_K) 

    call Matrix_Norm2(n,Matrix_inv,Norm_2_K_inv) 

    Condition_Num_Norm =Norm_2_K* Norm_2_K_inv
endif

if (norm_Type==1) then
    call Matrix_Norm1(n,Matrix,Norm_1_K) 

    call Matrix_Norm1(n,Matrix_inv,Norm_1_K_inv) 

    Condition_Num_Norm =Norm_1_K* Norm_1_K_inv
endif      

return
END subroutine Matrix_Condition_Number_v2

