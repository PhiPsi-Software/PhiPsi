 

subroutine Cal_Transform_Derivatives_Hex8(dN_dxi, dN_deta, dN_dzeta, inv_J, dN_dx, dN_dy, dN_dz)

use Global_Float_Type      
implicit none
real(kind=FT), intent(in) :: dN_dxi(8), dN_deta(8), dN_dzeta(8), inv_J(3,3)
real(kind=FT), intent(out) :: dN_dx(8), dN_dy(8), dN_dz(8)
integer :: i

do i = 1, 8
    dN_dx(i) = inv_J(1,1)*dN_dxi(i) + inv_J(1,2)*dN_deta(i) + inv_J(1,3)*dN_dzeta(i)
    dN_dy(i) = inv_J(2,1)*dN_dxi(i) + inv_J(2,2)*dN_deta(i) + inv_J(2,3)*dN_dzeta(i)
    dN_dz(i) = inv_J(3,1)*dN_dxi(i) + inv_J(3,2)*dN_deta(i) + inv_J(3,3)*dN_dzeta(i)
enddo
    
end subroutine Cal_Transform_Derivatives_Hex8
