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
 
subroutine Tool_Get_3D_Rect_by_n_Vector_Center_LVector_LW_Ratio(n_Vector,&
           c_Crack_Center,L,W,LongSide_Vecotr,Out_Coors)
! Obtain a 3D rectangle based on the normal vector, center coordinates, and rectangle side lengths.
! LongSide_Vector(3) is the vector pointing in the direction of the long side.
! L and W are the lengths of the long side and the short side.
! 2023-02-28. NEWFTU2023022802.
!   
!
!     P1-------------------------------P2
!      |                               |
!      |                               |
!      A------------Center-------------B
!      |                               |
!      |                               |
!     P4-------------------------------P3
!
!
! Long side vector P2P1
! Short-side vector P1P4
! The normal vector points from the center towards itself
!

use Global_Float_Type 
implicit none
real(kind=FT),intent(in) :: n_Vector(3),c_Crack_Center(3),L,W
real(kind=FT),intent(in) :: LongSide_Vecotr(3)
real(kind=FT),intent(out) :: Out_Coors(4,3)

real(kind=FT) W_Vector(3),P_A(3),P_B(3)
real(kind=FT) P1(3),P2(3),P3(3),P4(3)
real(kind=FT) LongSide_Vecotr_Norm(3)

LongSide_Vecotr_Norm = LongSide_Vecotr
call Vector_Normalize(3,LongSide_Vecotr_Norm)     

call Vector_Cross_Product_3(n_Vector,LongSide_Vecotr,W_Vector)   
call Vector_Normalize(3,W_Vector)    

P_A = c_Crack_Center + LongSide_Vecotr_Norm*L/TWO

P_B = c_Crack_Center - LongSide_Vecotr_Norm*L/TWO

P1 = P_A - W_Vector*W/TWO

P2 = P1 - LongSide_Vecotr_Norm*L

P3 = P2 + W_Vector*W

P4 = P_A + W_Vector*W/TWO

Out_Coors(1,1:3) = P1
Out_Coors(2,1:3) = P2
Out_Coors(3,1:3) = P3
Out_Coors(4,1:3) = P4

return 
end SUBROUTINE Tool_Get_3D_Rect_by_n_Vector_Center_LVector_LW_Ratio                     
