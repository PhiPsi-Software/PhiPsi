!-----------------------------------------------------------
! Brief: Compute 3D rotation matrix from original to crack coordinate frame
!
! Parameters:
!   Input:  i_V              - unused dimension flag
!   Input:  Vector_Origi_1..3 - 3-component vectors in original frame
!   Input:  Vector_Crack_1..3 - 3-component vectors in target crack frame
!   Output: ThetaX,ThetaY,ThetaZ - rotation angles about each axis
!   Output: T_Matrix         - 3x3 rotation matrix mapping original to crack frame
!
! Notes:   Builds T_Matrix row-by-row via dot products between paired vectors.
!-----------------------------------------------------------

subroutine Tool_ThetaX_ThetaY_ThetaZ_3D_rotation(i_V, Vector_Origi_1, Vector_Origi_2, Vector_Origi_3, Vector_Crack_1, &
Vector_Crack_2, Vector_Crack_3, ThetaX,ThetaY,ThetaZ,T_Matrix)
!     After the 3D spatial coordinate system is rotated, its vectors also transform accordingly. 
!     For example: from Vector_Origi_1(1:3) to Vector_Crack_1(1:3)
!     Vector_Origi_2(1:3) to Vector_Crack_2(1:3)
!     From Vector_Origi_3(1:3) to Vector_Crack_3(1:3)
!     The purpose of this subroutine is to determine the corresponding rotation matrix 
!     T_Matrix based on known transformations.

!......................
! Variable Declaration
!......................
use Global_Float_Type   
use Global_Crack_3D

implicit none
integer,intent(in)::i_V
real(kind=FT),intent(in):: Vector_Origi_1(3),Vector_Origi_2(3),Vector_Origi_3(3), &
Vector_Crack_1(3),Vector_Crack_2(3),Vector_Crack_3(3)
real(kind=FT),intent(out):: ThetaX,ThetaY,ThetaZ,T_Matrix(3,3)
real(kind=FT) c_T(3,3)      
T_Matrix(1:3,1:3) = ZR
ThetaX = ZR
ThetaY = ZR
ThetaZ = ZR


call Tool_Cos_Value_of_Vectors_a_and_b_3D(Vector_Origi_1, Vector_Crack_1,c_T(1,1))
call Tool_Cos_Value_of_Vectors_a_and_b_3D(Vector_Origi_1, Vector_Crack_2,c_T(1,2))
call Tool_Cos_Value_of_Vectors_a_and_b_3D(Vector_Origi_1, Vector_Crack_3,c_T(1,3))
call Tool_Cos_Value_of_Vectors_a_and_b_3D(Vector_Origi_2, Vector_Crack_1,c_T(2,1))
call Tool_Cos_Value_of_Vectors_a_and_b_3D(Vector_Origi_2, Vector_Crack_2,c_T(2,2))
call Tool_Cos_Value_of_Vectors_a_and_b_3D(Vector_Origi_2, Vector_Crack_3,c_T(2,3))
call Tool_Cos_Value_of_Vectors_a_and_b_3D(Vector_Origi_3, Vector_Crack_1,c_T(3,1))
call Tool_Cos_Value_of_Vectors_a_and_b_3D(Vector_Origi_3, Vector_Crack_2,c_T(3,2))
call Tool_Cos_Value_of_Vectors_a_and_b_3D(Vector_Origi_3, Vector_Crack_3,c_T(3,3))
T_Matrix = c_T


return 
end subroutine Tool_ThetaX_ThetaY_ThetaZ_3D_rotation                   
