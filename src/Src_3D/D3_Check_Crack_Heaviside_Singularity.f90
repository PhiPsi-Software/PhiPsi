!-----------------------------------------------------------
! Brief: Detect Heaviside enrichment singularity for a single 3D crack.
!
! Parameters:
!   Input:  i_C              - crack index
!           Singularity_Factor - threshold ratio for min/max signed distance
!   Output: Logical_Yes_Singu - true if any element has a degenerate
!                                |min/max| signed-distance ratio
!
! Notes:   Walks every node, computes signed distance via
!          D3_Get_Signed_Dis_to_Crack_Mesh, then flags elements whose
!          node signed distances straddle zero with an extreme ratio.
!-----------------------------------------------------------

SUBROUTINE D3_Check_Crack_Heaviside_Singularity(i_C,Singularity_Factor,Logical_Yes_Singu)
! Determine whether the Heaviside enrichment function distance for a certain crack is singular.
!2022-08-01.


!//////////////////////////////////
! Read the public variable module.
!//////////////////////////////////
use Global_Float_Type
use Global_Common
use Global_Model
use Global_Elem_Area_Vol
use Global_Crack_Common
use Global_Crack_3D
!
use Global_INTERFACE_D3_Get_Signed_Dis_to_Crack_Mesh



!////////////////////////////
! Variable type declaration.
!////////////////////////////
implicit none
integer,intent(in)::i_C
real(kind=FT),intent(in)::Singularity_Factor
logical,intent(out)::Logical_Yes_Singu
integer i_Node,i_E
real(kind=FT) c_Node_Coor(3),c_Min_Signed_Dis,Check_Ball_R
logical Yes_Found_Min_Signed_Dis
real(kind=FT) c_n_Vector(3)
real(kind=FT) tem_PER_Node_to_FS(3)
logical tem_Yes_Node_PER_in_FS(num_Node)
real(kind=FT) tem_Dis_Node_to_FS(num_Node)
integer c_NN(8)
real(kind=FT) c_max,c_min
real(kind=FT) c_Singularity
real(kind=FT) c_Signed_Dis_v2

!/////////////////
! Initialization.
!/////////////////
Logical_Yes_Singu = .False.

!////////////////////////////////////
! The signed distance from the 
! computational node to the fracture 
! surface.
!////////////////////////////////////
! Loop through all nodes
do i_Node = 1,num_Node
    c_Node_Coor = Coor(i_Node,1:3)
    c_Min_Signed_Dis  = Con_Big_20
    ! Calculate the signed distance from the current node to the fracture plane and determine the foot
    ! position
    Check_Ball_R  = 3.0D0*Ave_Elem_L
    call D3_Get_Signed_Dis_to_Crack_Mesh(c_Node_Coor,i_C,Check_Ball_R,tem_Dis_Node_to_FS(i_Node),&
    c_Signed_Dis_v2,tem_Yes_Node_PER_in_FS(i_Node),tem_PER_Node_to_FS(1:3),          &
    Yes_Found_Min_Signed_Dis,c_n_Vector)
enddo

!///////////////////////////////////
! Looking for elements with unusual 
! symbol distances.
!///////////////////////////////////
do i_E = 1,Num_Elem
    c_NN  = G_NN(1:8,i_E)
    c_max = maxval(tem_Dis_Node_to_FS(c_NN))
    c_min = minval(tem_Dis_Node_to_FS(c_NN))
    ! If the element's nodes include both positive and negative nodes
    if( c_max > Tol_11 .and. c_min < -Tol_11)then
        ! Strange Factor
        c_Singularity = abs(c_min/c_max)
        if(c_Singularity>=Singularity_Factor .or. c_Singularity<= ONE/Singularity_Factor) then
            Logical_Yes_Singu = .True.
            ! exit    !Exit the loop
            return
        endif
    endif
enddo

RETURN
END SUBROUTINE D3_Check_Crack_Heaviside_Singularity
