!-----------------------------------------------------------
! Brief: Determine average outward normal of a crack outline polygon.
!
! Parameters:
!   Input:  i_C                 - crack index
!           Outline_num         - number of outline edges
!           c_Outline           - outline edge list (Outline_num x 2)
!   Output: Cros_Product_Vector - averaged cross-product (outward normal)
!
! Notes:   Sums edge-to-center cross products and normalises the result.
!-----------------------------------------------------------

subroutine D3_Get_Crack_Mesh_Outline_Clockwise(i_C, Outline_num,c_Outline,Cros_Product_Vector)
!     Determine the clockwise or counterclockwise direction of the outline
!     using the cross product. Ref: My PhiPsi Development Notebook V1-P41.
!     2022-07-13.

!     ----------------------------
!     Read public variable module
!     ----------------------------
use Global_Float_Type
use Global_Common
use Global_Model
use Global_Elem_Area_Vol
use Global_Crack_3D

!     --------------------------
!     Variable Type Declaration
!     --------------------------
implicit none
integer,intent(in)::i_C,Outline_num
integer,intent(in)::c_Outline(Outline_num,2)
real(kind=FT),intent(out)::Cros_Product_Vector(3)
integer :: i_outline,c_P1,c_P2,c_P3
real(kind=FT) c_Vector_1(3),c_Vector_2(3)
real(kind=FT) c_P1_Coor(3),c_P2_Coor(3),c_P3_Coor(3)
real(kind=FT) c_Cros_Product(3)
real(kind=FT) Sum_Cros_Product(3)
real(kind=FT) c_Center(3)

!     ---------------------
!     Main program section
!     ---------------------
!Cros_Product_Vector = ZR
Sum_Cros_Product(1:3) =  ZR

! Obtain the centroid of the boundary discrete points.
c_Center(1) = sum(Crack3D_Meshed_Node(i_C)%row( c_Outline(1:Outline_num,1),1))/Outline_num
c_Center(2) = sum(Crack3D_Meshed_Node(i_C)%row( c_Outline(1:Outline_num,1),2))/Outline_num
c_Center(3) = sum(Crack3D_Meshed_Node(i_C)%row( c_Outline(1:Outline_num,1),3))/Outline_num

do i_outline =1,Outline_num
    c_P1 = c_Outline(i_outline,1)
    c_P2 = c_Outline(i_outline,2)
    !c_P3 = c_Outline(i_outline+1,2)
    c_P1_Coor = Crack3D_Meshed_Node(i_C)%row(c_P1,1:3)
    c_P2_Coor = Crack3D_Meshed_Node(i_C)%row(c_P2,1:3)
    !c_P3_Coor = Crack3D_Meshed_Node(i_C,c_P3,1:3)
    c_Vector_1 = c_P2_Coor - c_P1_Coor
    !c_Vector_2 = c_P3_Coor - c_P2_Coor
    c_Vector_2 = c_Center - c_P1_Coor
    call Vector_Normalize(3,c_Vector_1)
    call Vector_Normalize(3,c_Vector_2)
call Vector_Cross_Product_3(c_Vector_1,c_Vector_2, c_Cros_Product(1:3))
    call Vector_Normalize(3,c_Cros_Product)
    Sum_Cros_Product(1:3) = Sum_Cros_Product(1:3) +c_Cros_Product
enddo
Cros_Product_Vector = Sum_Cros_Product(1:3)/Outline_num

call Vector_Normalize(3,Cros_Product_Vector)

return
end subroutine D3_Get_Crack_Mesh_Outline_Clockwise
