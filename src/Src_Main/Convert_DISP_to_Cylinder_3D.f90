!-----------------------------------------------------------
! Brief: Transform nodal displacements to cylindrical CS.
!
! Parameters:
!   Input:  n_FD - length of the displacement vector
!           in_Disp - Cartesian-system displacement vector
!   Output: out_Disp - cylindrical-system displacement vector
!
! Notes:   Uses per-node rotation angle from
!          Theta_Cartesian_to_Cylinder_Node for 3D conversion.
!-----------------------------------------------------------

subroutine Convert_DISP_to_Cylinder_3D(n_FD,in_Disp,out_Disp)
!     Transform node displacements from the Cartesian coordinate system to the cylindrical coordinate system. 3D.
!     Ref: \theory_documents\024.2 Stress Tensor and Strain Tensor Transformation from Cartesian Coordinate System to Cylindrical Coordinate System --- 2.4-2021-09-10.html
!     2021-09-11.

!     ----------------------------
!     Read public variable module
!     ----------------------------
use Global_Float_Type
use Global_Common
use Global_Model

!     --------------------------
!     Variable Type Declaration
!     --------------------------
implicit none
integer,intent(in)::n_FD
real(kind=FT),intent(in)::  in_Disp(n_FD)
real(kind=FT),intent(out):: out_Disp(n_FD)
integer i_E,i_N
real(kind=FT) c_Theta,c_Disp(1:3),sTheta,cTheta
real(kind=FT) :: tem_Disp(1:3)

print *,'    Converting disp from Cartesian CS cylindical CS...'
do i_N = 1,Num_Node
    c_Disp(1)  = in_Disp(3*i_N-2)
    c_Disp(2)  = in_Disp(3*i_N-1)
    c_Disp(3)  = in_Disp(3*i_N-0)

    c_Theta = Theta_Cartesian_to_Cylinder_Node(i_N)

    sTheta      =  sin(c_Theta)
    cTheta      =  cos(c_Theta)
    tem_Disp(1) =  cTheta*c_Disp(1) + sTheta*c_Disp(2)
    tem_Disp(2) = -sTheta*c_Disp(1) + cTheta*c_Disp(2)
    tem_Disp(3) = c_Disp(3)

    out_Disp(3*i_N-2) = tem_Disp(1)
    out_Disp(3*i_N-1) = tem_Disp(2)
    out_Disp(3*i_N-0) = tem_Disp(3)

end do
return
END subroutine Convert_DISP_to_Cylinder_3D
