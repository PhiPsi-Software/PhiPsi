!-----------------------------------------------------------
! Brief: Compute the in-situ normal stress on a crack surface from
!        the Gauss-point stress tensor at a given element.
!
! Parameters:
!   Input:  c_ELE        - element index
!           ori_n        - crack surface normal direction (3-vector)
!   Output: Normal_Stress - normal traction on the surface
!
! Notes:   Uses the in-situ stress tensor stored at the first
!          Gauss point and rotates it via Cauchy formula.
!-----------------------------------------------------------

subroutine Get_Normal_Insitu_Stress(c_ELE,ori_n,Normal_Stress)
! Obtain the ground stress in the normal direction of the crack surface. 2026-01-23.

use Global_Float_Type
use Global_Stress
use Global_Common, only: Key_InSitu_Strategy

implicit none
integer, intent(in) :: c_ELE
real(kind=FT), intent(in) :: ori_n(3)
real(kind=FT), intent(out) :: Normal_Stress
real(kind=FT) :: c_Insitu_Stress(6)       

if(Key_InSitu_Strategy == 4 )then
    ! Obtained in-situ stress tensor c_Insitu_Stress
    c_Insitu_Stress(1) =InSitu_Strs_Gaus_xx(c_ELE,1)
    c_Insitu_Stress(2) =InSitu_Strs_Gaus_yy(c_ELE,1)                    
    c_Insitu_Stress(3) =InSitu_Strs_Gaus_zz(c_ELE,1)                    
    c_Insitu_Stress(4) =InSitu_Strs_Gaus_xy(c_ELE,1)                    
    c_Insitu_Stress(5) =InSitu_Strs_Gaus_yz(c_ELE,1)                    
    c_Insitu_Stress(6) =InSitu_Strs_Gaus_xz(c_ELE,1)                    
    ! Calculate the normal stress on the fracture surface
    call Tool_Get_Normal_Stress_on_Plane_by_Stress_Tensor(c_Insitu_Stress,ori_n,Normal_Stress)
endif
end subroutine Get_Normal_Insitu_Stress
