!-----------------------------------------------------------
! Brief: Snap a crack-front extension point onto the crack plane.
!
! Parameters:
!   In/Out: point_inout - extension point corrected in place
!   Input:  i_C         - crack index
!
! Notes:   Replaces the point with its perpendicular foot on the
!          plane of the first meshed crack-surface triangle.
!-----------------------------------------------------------

subroutine D3_In_Plane_Growth(point_inout,i_C)
!     In-plane extension. Correct the extension point at the crack front
!     vertex to the foot of the perpendicular from that point to the plane of the crack.
!     Added date: 2022-04-18, NEWFTU2022041801


!     **************
!     Public Module
!     **************
use Global_Crack_3D
use Global_Model
use Global_Common

!     *********************
!     Variable Declaration
!     *********************
implicit none
real(kind=FT),intent(inout)::point_inout(3)
integer,intent(in)::i_C
real(kind=FT) c_Point(3),c_PER(3),c_Distance
logical c_Yes_PER_in,c_Yes_PER_on
real(kind=FT) c_Tri_P1(3),c_Tri_P2(3),c_Tri_P3(3)
integer c_CalP_1,c_CalP_2,c_CalP_3
integer Crack_Node1,Crack_Node2,Crack_Node3


c_Point = point_inout


!*****************************************
!       OPTION - 1
! There is a bug. For example, when three
! fluid nodes are on a straight line.
!*****************************************
!     Select the three nodes of the first fluid element of the crack as the projection surface.
!     c_CalP_1 = Cracks_FluidEle_CalP_3D(i_C,1,1) 
!     c_CalP_2 = Cracks_FluidEle_CalP_3D(i_C,1,2) 
!     c_CalP_3 = Cracks_FluidEle_CalP_3D(i_C,1,3) 
!     Coordinates of 3 fluid nodes
!     c_Tri_P1 = Cracks_CalP_Coors_3D(i_C,c_CalP_1,1:3) 
!     c_Tri_P2 = Cracks_CalP_Coors_3D(i_C,c_CalP_2,1:3) 
!     c_Tri_P3 = Cracks_CalP_Coors_3D(i_C,c_CalP_3,1:3) 

!*****************************************
!       OPTION - 2
! OPTION - 1 has a bug. BUGFIX2022070801.
!*****************************************
! Select the three nodes of the first discrete crack surface element
! of the crack as the projection plane. 2022-07-08.
Crack_Node1 = Crack3D_Meshed_Ele(i_C)%row(1,1)
Crack_Node2 = Crack3D_Meshed_Ele(i_C)%row(1,2)
Crack_Node3 = Crack3D_Meshed_Ele(i_C)%row(1,3)
c_Tri_P1 = Crack3D_Meshed_Node(i_C)%row(Crack_Node1,1:3)
c_Tri_P2 = Crack3D_Meshed_Node(i_C)%row(Crack_Node2,1:3)
c_Tri_P3 = Crack3D_Meshed_Node(i_C)%row(Crack_Node3,1:3)


! Calculate the foot of the perpendicular projection
call Tool_Dis_Point_to_3D_Tri(c_Point, c_Tri_P1,c_Tri_P2,c_Tri_P3, c_Distance,c_PER,c_Yes_PER_in,c_Yes_PER_on)
! Change the Final_Point coordinates to the foot of the perpendicular
point_inout = c_PER

return
end subroutine D3_In_Plane_Growth
