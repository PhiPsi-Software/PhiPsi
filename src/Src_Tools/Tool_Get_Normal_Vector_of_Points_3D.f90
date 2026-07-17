!-----------------------------------------------------------
! Brief: Estimate the plane normal of a set of 3D points via least squares.
!
! Parameters:
!   Input:  num_Point - number of input points
!   Input:  In_Points - 3D coordinates of the points
!   Output: n_Vector(3) - unit normal vector of the fitted plane
!
! Notes:   Solves (M^T M) n = M^T 1 by inverting the 3x3 normal
!   equations matrix and normalizing the resulting vector.
!-----------------------------------------------------------

subroutine Tool_Get_Normal_Vector_of_Points_3D(In_Points,num_Point,n_Vector)
! Calculate the outward normal vector of the plane on which the 3D point In_Points lies.
! For information on obtaining an out-of-plane normal vector, refer to the following materials
!Theory Ref-https://www.freesion.com/article/7157945105/
!            or
! \theory_documents\026 Optimal Spatial Circle Fitting for 3D Discrete Points and
! Implementation-2021-10-30.pdf
!2021-10-30.

use Global_Float_Type

implicit none
integer,intent(in)::num_Point
real(kind=FT),intent(in)::In_Points(num_Point,3)
real(kind=FT),intent(out)::n_Vector(3)
real(kind=FT) M(num_Point,3),MTM_Inv(3,3),MTM(3,3),L(num_Point)
real(kind=FT) norm_n_Vector

M = In_Points
MTM = matmul(TRANSPOSE(M),M)
call Matrix_Inverse_3x3(MTM,MTM_Inv)  
L(1:num_Point) = ONE
n_Vector =  matmul(matmul(MTM_Inv,TRANSPOSE(M)),L)
call Vector_Norm2(3,n_Vector,norm_n_Vector) 
n_Vector =  n_Vector /norm_n_Vector


return 
end SUBROUTINE Tool_Get_Normal_Vector_of_Points_3D                        
