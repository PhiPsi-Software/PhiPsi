!-----------------------------------------------------------
! Brief: Compute the Jacobian and its determinant for a HEX8.
!
! Parameters:
!   Input:  X(8), Y(8), Z(8)     - HEX8 nodal coordinates
!           dN_dxi(8), dN_deta(8), dN_dzeta(8) - shape derivs
!   Output: J_mat(3,3)           - Jacobian matrix
!           det_J                 - Jacobian determinant
!
! Notes:   Standard 3x3 Jacobian assembly via the chain rule and
!          cofactor-expansion determinant for HEX8 finite
!          elements.
!-----------------------------------------------------------

subroutine Cal_Jacobian_3D_Hex8(X, Y, Z, dN_dxi, dN_deta, dN_dzeta, J_mat, det_J)
!2026-01-28.

use Global_Float_Type      
implicit none

real(kind=FT), intent(in) :: X(8), Y(8), Z(8)
real(kind=FT), intent(in) :: dN_dxi(8), dN_deta(8), dN_dzeta(8)
real(kind=FT), intent(out) :: J_mat(3,3), det_J

J_mat(1,1) = sum(dN_dxi * X)
J_mat(1,2) = sum(dN_dxi * Y)
J_mat(1,3) = sum(dN_dxi * Z)
J_mat(2,1) = sum(dN_deta * X)
J_mat(2,2) = sum(dN_deta * Y)
J_mat(2,3) = sum(dN_deta * Z)
J_mat(3,1) = sum(dN_dzeta * X)
J_mat(3,2) = sum(dN_dzeta * Y)
J_mat(3,3) = sum(dN_dzeta * Z)

det_J = J_mat(1,1)*(J_mat(2,2)*J_mat(3,3) - J_mat(2,3)*J_mat(3,2)) - &
J_mat(1,2)*(J_mat(2,1)*J_mat(3,3) - J_mat(2,3)*J_mat(3,1)) + J_mat(1,3)*(J_mat(2,1)*J_mat(3,2) - J_mat(2,2)*J_mat(3,1))
    
end subroutine Cal_Jacobian_3D_Hex8
