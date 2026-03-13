!     ================================================= !
!             ____  _       _   ____  _____   _         !
!            |  _ \| |     |_| |  _ \|  ___| |_|        !
!            | |_) | |___   _  | |_) | |___   _         !
!            |  _ /|  _  | | | |  _ /|___  | | |        !
!            | |   | | | | | | | |    ___| | | |        !
!            |_|   |_| |_| |_| |_|   |_____| |_|        !
!     ================================================= !
!     PhiPsi:     a general-purpose computational       !
!                 mechanics program written in Fortran. !
!     Website:    http://phipsi.top                     !
!     Author:     Shi Fang, Huaiyin Institute of        !
!                 Technology, Huaian, JiangSu, China    !
!     Email:      shifang@hyit.edu.cn                   !
!     ------------------------------------------------- !
!     Please cite the following papers:                 !
!     (1)Shi F., Lin C. Modeling fluid-driven           !
!        propagation of 3D complex crossing fractures   !
!        with the extended finite element method.       !
!        Computers and Geotechnics, 2024, 172, 106482.  !
!     (2)Shi F., Wang D., Li H. An XFEM-based approach  !
!        for 3D hydraulic fracturing simulation         !
!        considering crack front segmentation. Journal  !
!        of Petroleum Science and Engineering, 2022,    !
!        214, 110518.                                   !
!     (3)Shi F., Wang D., Yang Q. An XFEM-based         !
!        numerical strategy to model three-dimensional  !
!        fracture propagation regarding crack front     !
!        segmentation. Theoretical and Applied Fracture !
!        Mechanics, 2022, 118, 103250.                  !
!     (4)Shi F., Liu J. A fully coupled hydromechanical !
!        XFEM model for the simulation of 3D non-planar !
!        fluid-driven fracture propagation. Computers   !
!        and Geotechnics, 2021, 132: 103971.            !
!     (5)Shi F., Wang X.L., Liu C., Liu H., Wu H.A. An  !
!        XFEM-based method with reduction technique     !
!        for modeling hydraulic fracture propagation    !
!        in formations containing frictional natural    !
!        fractures. Engineering Fracture Mechanics,     !
!        2017, 173: 64-90.                              !
!     ------------------------------------------------- !
 

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
