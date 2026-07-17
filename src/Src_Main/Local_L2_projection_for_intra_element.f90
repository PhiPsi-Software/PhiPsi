!-----------------------------------------------------------
! Brief: 3D history-preserving local L2 projection that transfers
!        crack-tip enriched DOFs from the old tip basis to the
!        updated tip basis when the crack tip moves inside the same
!        element.
!
! Parameters:
!   Input:  isub                   - sub-step index
!           crack_number           - crack index
!           num_gp                 - Gauss points per direction
!           c_DOFs_Value           - global DOF vector
!           number_nodes_projected - number of tip-enriched nodes
!           nodes_to_be_projected  - node IDs
!   Output: nodes_new_enriched_value - new 3-component tip DOFs
!
! Notes:   Solves M_tt * b_new = P_tt * b_old; only the first tip
!          branch function is used.
!-----------------------------------------------------------

subroutine Local_L2_projection_for_intra_element( isub, crack_number, num_gp, c_DOFs_Value, number_nodes_projected, &
nodes_to_be_projected, nodes_new_enriched_value)
!-----------------------------------------------------------------------
! Subroutine Local_L2_projection_for_intra_element
!
! Purpose:
!   Perform a history-preserving local L2 projection of crack-tip enriched
!   DOFs during intra-element crack-tip motion.
!
!   Unlike the inter-element case, the enriched node set remains unchanged.
!   However, because the crack-tip enrichment function is defined with
!   respect to the instantaneous crack-tip position, the old crack-tip
!   enriched basis and the updated crack-tip enriched basis are different.
!   Therefore, the old enriched DOFs cannot be directly reused after the
!   crack tip moves inside the same element.
!
!   This routine implements the local transfer described in the paper:
!
!       M_tt * b_new = P_tt * b_old
!
!   where
!     - b_old : vector of old crack-tip enriched displacement/velocity/
!               acceleration DOFs (collected from c_DOFs_Value)
!     - b_new : resulting crack-tip enriched DOFs in the updated tip basis
!
!   The matrices M_tt and P_tt are assembled by numerical integration over
!   the support domain (union of elements containing the projected nodes).
!   The same operator is applied to displacement, velocity and acceleration
!   to preserve kinematic consistency.
!
! Assumptions:
!   - Only a single crack-tip enrichment function F (the first branch
!     function, typically sqrt(r)*sin(theta/2)) is used.
!   - The crack is represented in 3D; nodes are enriched with a single
!     tip enrichment DOF set per node.
!   - The updated crack-tip position and orientation are stored in global
!     data structures (Solid_El_Tip_BaseLine_*).
!   - The old crack-tip local system is obtained from global OLD arrays.
!
! Usage example:
!  crack_number   = 1
!  number_nodes_projected   = 8
!  if(allocated(nodes_to_be_projected)) deallocate(nodes_to_be_projected)
!  allocate(nodes_to_be_projected(number_nodes_projected)) 
!  if(allocated(nodes_new_enriched_value)) deallocate(nodes_new_enriched_value)
!  allocate(nodes_new_enriched_value(number_nodes_projected*3))
!  nodes_to_be_projected(1:number_nodes_projected) = [265,625,626,266,279,280,639,640]
!  nodes_new_enriched_value(1:number_nodes_projected*3) = ZR
!  num_gp = 14   !2,3,4,5,6,7,8,9,10,11,12,13,14,15
!  call Local_L2_projection_for_intra_element(isub,crack_number,num_gp,&
!                                             DISP(1:Total_FD),number_nodes_projected, &
!                                             nodes_to_be_projected(1:number_nodes_projected), &
!                                             nodes_new_enriched_value(1:number_nodes_projected*3))
!
! Arguments:
!   isub                       (in)  integer : sub-stepping index
!   crack_number               (in)  integer : which crack is propagating
!   num_gp                     (in)  integer : number of Gauss points per direction
!   c_DOFs_Value               (in)  real(FT) : global displacement (or velocity/acceleration)
!                                              vector; contains both standard and enriched DOFs.
!   number_nodes_projected     (in)  integer : number of crack-tip enriched nodes whose
!                                              enriched DOFs are transferred from the old
!                                              tip basis to the updated tip basis.
!   nodes_to_be_projected      (in)  integer(:) : list of crack-tip enriched node IDs.
!   nodes_new_enriched_value   (out) real(FT) : output array of length
!                                              number_nodes_projected * 3,
!                                              containing the new crack-tip enriched
!                                              DOFs (x, y, z components interleaved).
!
! Global data dependencies (modules):
!   Global_Float_Type, Global_Common, Global_Model,
!   Global_Elem_Area_Vol, Global_Crack_Common, Global_Crack_3D,
!   Global_Inter_Cal_N_dNdkesi_J_detJ_3D
!
! Notes:
!   - The projection domain is the union of all elements that contain any
!     of the nodes in nodes_to_be_projected.
!   - The routine checks that each projected node is currently a crack-tip
!     enriched node.
!   - Numerical integration uses Gauss-Legendre quadrature.
!   - The updated crack-tip local system is read from the current global
!     crack-tip data structure.
!   - The old crack-tip local system is read from the global OLD arrays.
!
! Implementation references:
!   Section 2.7 and Section 2.9 of
!   "A history-preserving enrichment-transition strategy for explicit
!    dynamic XFEM simulation of crack propagation"
!
! Author:  Fang Shi.
! Date:    2026-05-14.
! Website: phipsi.top
!-----------------------------------------------------------------------

use Global_Float_Type
use Global_Common
use Global_Model
use Global_Elem_Area_Vol
use Global_Crack_Common
use Global_Crack_3D
use Global_Inter_Cal_N_dNdkesi_J_detJ_3D
implicit none

integer, intent(in) :: isub, crack_number, num_gp
integer, intent(in) :: number_nodes_projected
integer, intent(in) :: nodes_to_be_projected(number_nodes_projected)
real(kind=FT), intent(in) :: c_DOFs_Value(Total_FD)
real(kind=FT), intent(out) :: nodes_new_enriched_value(number_nodes_projected*3)

integer :: i, j, k
integer :: c_node, node_j
integer :: g1, g2, g3
integer :: loc_i, loc_j
integer :: ref_elem, c_Cr_Location
integer :: mat_num, c_mat_type
integer :: c_NN(8)
integer :: temp_vector(Solid_El_Max_num_Crs)
real(kind=FT), allocatable :: b_old(:), b_new(:), rhs_vec(:)
real(kind=FT), allocatable :: M_tt(:,:), P_tt(:,:), M_tt_inv(:,:)
real(kind=FT) :: c_X_NODES(8), c_Y_NODES(8), c_Z_NODES(8)
real(kind=FT), allocatable  :: gp_loc(:), gp_w(:)
real(kind=FT) :: kesi, yita, zeta, weight
real(kind=FT) :: detJ
real(kind=FT) :: Jmat(3,3), Nmat(3,24), dNdkesi(8,3)
real(kind=FT) :: N8(8)
real(kind=FT) :: scalar_M, scalar_P
real(kind=FT) :: phi_tip_new_i, phi_tip_new_j, phi_tip_old_j
real(kind=FT) :: dFdx_tip(4)
real(kind=FT) :: dFdy_tip(4)
real(kind=FT) :: dFdz_tip(4)
real(kind=FT) :: F_Gauss_new(4), F_Node_new_i(4), F_Node_new_j(4)
real(kind=FT) :: F_Gauss_old(4), F_Node_old_j(4)
real(kind=FT) :: BaseLine_A(3), BaseLine_B(3), BaseLine_Mid(3)
real(kind=FT) :: BaseLine_x_Vec(3), BaseLine_y_Vec(3), BaseLine_z_Vec(3)
real(kind=FT) :: c_T_Matrx(3,3)
real(kind=FT) :: BaseLine_A_OLD(3), BaseLine_B_OLD(3), BaseLine_Mid_OLD(3)
real(kind=FT) :: BaseLine_x_Vec_OLD(3), BaseLine_y_Vec_OLD(3), BaseLine_z_Vec_OLD(3)
real(kind=FT) :: c_T_Matrx_OLD(3,3)
real(kind=FT) :: Node_Point(3), Node_Coor_Local_New(3), Node_Coor_Local_Old(3)
real(kind=FT) :: Gauss_Coor_Local_New(3), Gauss_Coor_Local_Old(3)
real(kind=FT) :: Global_coor_Gauss(3)
real(kind=FT) :: Temp_Point(3)
real(kind=FT) :: r_Node_new, theta_Node_new, r_Gauss_new, theta_Gauss_new
real(kind=FT) :: r_Node_old, theta_Node_old, r_Gauss_old, theta_Gauss_old
real(kind=FT) :: omega
real(kind=FT) :: det_Mtt, Condition_Num_Norm
integer :: i_E, element_number
integer :: local_pos_of_proj_nodes(number_nodes_projected)
integer :: Temp_Domain_Elements(1000), Uniqued_Temp_Domain_Elements(1000), i_E_count
integer :: number_of_elements
integer, allocatable :: Domain_Elements_List(:)
integer i_N

print *,"    Perform L2 projection of crack-tip enriched DOFs onto updated crack-tip enriched DOFs..."

!****************************
! Check the projected nodes.
!****************************
do i = 1, number_nodes_projected
    c_node = nodes_to_be_projected(i)
    if (Enriched_Node_Type_3D(c_node, crack_number) /= 1) then
        print *, '    ERROR-2026051402 :: projected node is not tip enriched node.'
        print *, '    In Local_L2_projection_for_intra_element.'
        print *, '    Node number: ', c_node
        call Warning_Message('S', Keywords_Blank)
    end if
end do

!****************************
! Initialization.
!****************************
nodes_new_enriched_value(1:number_nodes_projected*3) = ZR
allocate(b_old(number_nodes_projected*3))
allocate(b_new(number_nodes_projected*3))
allocate(rhs_vec(number_nodes_projected*3))
allocate(M_tt(number_nodes_projected*3, number_nodes_projected*3))
allocate(P_tt(number_nodes_projected*3, number_nodes_projected*3))
allocate(M_tt_inv(number_nodes_projected*3, number_nodes_projected*3))
b_old    = ZR
b_new    = ZR
rhs_vec  = ZR
M_tt     = ZR
P_tt     = ZR
M_tt_inv = ZR

!*****************************************************
! Get the support domain C of all tip-enriched nodes.
!*****************************************************
Temp_Domain_Elements(1:1000) = 0
i_E_count = 0
do i = 1, number_nodes_projected
    c_node = nodes_to_be_projected(i)
Temp_Domain_Elements(i_E_count+1:i_E_count + num_Node_Elements(c_node)) = &
Node_Elements_3D(c_node)%row(1:num_Node_Elements(c_node))
    i_E_count = i_E_count + num_Node_Elements(c_node)
end do

call Vector_Unique_Int(1000, i_E_count, Temp_Domain_Elements(1:1000), &
Uniqued_Temp_Domain_Elements(1:1000), number_of_elements)

allocate(Domain_Elements_List(number_of_elements))
Domain_Elements_List(1:number_of_elements) = Uniqued_Temp_Domain_Elements(1:number_of_elements)

!*****************************************
! Get old tip-enriched nodal DOFs values.
!*****************************************
do i = 1, number_nodes_projected
    c_node = nodes_to_be_projected(i)
    b_old((i-1)*3 + 1) = c_DOFs_Value(c_POS_3D(c_node, crack_number)*3 - 2)
    b_old((i-1)*3 + 2) = c_DOFs_Value(c_POS_3D(c_node, crack_number)*3 - 1)
    b_old((i-1)*3 + 3) = c_DOFs_Value(c_POS_3D(c_node, crack_number)*3)
end do

print *,"    -------------------------------------"
print *,"    DOFs of old tip-enirched nodes:      "
print *,"    -------------------------------------"
do i = 1, number_nodes_projected
write(*,'("     Node=", I6, ", X=", F12.6, ", Y=", F12.6, ", Z=", F12.6)') nodes_to_be_projected(i), &
b_old((i-1)*3 + 1), b_old((i-1)*3 + 2), b_old((i-1)*3 + 3)
end do

!***********************************
! Prepare the gauss point.
!***********************************
allocate(gp_loc(num_gp), gp_w(num_gp))
select case(num_gp)
case (2)
    gp_loc(1) = -0.5773502691896258D0
    gp_loc(2) =  0.5773502691896258D0
    gp_w(1) = 1.0000000000000000D0
    gp_w(2) = 1.0000000000000000D0
case (3)
    gp_loc(1) = -0.7745966692414834D0
    gp_loc(2) =  0.0000000000000000D0
    gp_loc(3) =  0.7745966692414834D0
    gp_w(1) = 0.5555555555555556D0
    gp_w(2) = 0.8888888888888889D0
    gp_w(3) = 0.5555555555555556D0
case (4)
    gp_loc(1) = -0.8611363115940526D0
    gp_loc(2) = -0.3399810435848563D0
    gp_loc(3) =  0.3399810435848563D0
    gp_loc(4) =  0.8611363115940526D0
    gp_w(1) = 0.3478548451374538D0
    gp_w(2) = 0.6521451548625461D0
    gp_w(3) = 0.6521451548625461D0
    gp_w(4) = 0.3478548451374538D0
case (5)
    gp_loc(1) = -0.9061798459386640D0
    gp_loc(2) = -0.5384693101056831D0
    gp_loc(3) =  0.0000000000000000D0
    gp_loc(4) =  0.5384693101056831D0
    gp_loc(5) =  0.9061798459386640D0
    gp_w(1) = 0.2369268850561891D0
    gp_w(2) = 0.4786286704993665D0
    gp_w(3) = 0.5688888888888889D0
    gp_w(4) = 0.4786286704993665D0
    gp_w(5) = 0.2369268850561891D0
case (6)
    gp_loc(1) = -0.9324695142031521D0
    gp_loc(2) = -0.6612093864662645D0
    gp_loc(3) = -0.2386191860831969D0
    gp_loc(4) =  0.2386191860831969D0
    gp_loc(5) =  0.6612093864662645D0
    gp_loc(6) =  0.9324695142031521D0
    gp_w(1) = 0.1713244923791704D0
    gp_w(2) = 0.3607615730481386D0
    gp_w(3) = 0.4679139345726910D0
    gp_w(4) = 0.4679139345726910D0
    gp_w(5) = 0.3607615730481386D0
    gp_w(6) = 0.1713244923791704D0
case (7)
    gp_loc(1) = -0.9491079123427585D0
    gp_loc(2) = -0.7415311855993945D0
    gp_loc(3) = -0.4058451513773972D0
    gp_loc(4) =  0.0000000000000000D0
    gp_loc(5) =  0.4058451513773972D0
    gp_loc(6) =  0.7415311855993945D0
    gp_loc(7) =  0.9491079123427585D0
    gp_w(1) = 0.1294849661688697D0
    gp_w(2) = 0.2797053914892766D0
    gp_w(3) = 0.3818300505051189D0
    gp_w(4) = 0.4179591836734694D0
    gp_w(5) = 0.3818300505051189D0
    gp_w(6) = 0.2797053914892766D0
    gp_w(7) = 0.1294849661688697D0
case (8)
    gp_loc(1) = -0.9602898564975363D0
    gp_loc(2) = -0.7966664774136267D0
    gp_loc(3) = -0.5255324099163290D0
    gp_loc(4) = -0.1834346424956498D0
    gp_loc(5) =  0.1834346424956498D0
    gp_loc(6) =  0.5255324099163290D0
    gp_loc(7) =  0.7966664774136267D0
    gp_loc(8) =  0.9602898564975363D0
    gp_w(1) = 0.1012285362903763D0
    gp_w(2) = 0.2223810344533745D0
    gp_w(3) = 0.3137066458778873D0
    gp_w(4) = 0.3626837833783620D0
    gp_w(5) = 0.3626837833783620D0
    gp_w(6) = 0.3137066458778873D0
    gp_w(7) = 0.2223810344533745D0
    gp_w(8) = 0.1012285362903763D0
case (9)
    gp_loc(1) = -0.9681602395076261D0
    gp_loc(2) = -0.8360311073266358D0
    gp_loc(3) = -0.6133714327005904D0
    gp_loc(4) = -0.3242534234038089D0
    gp_loc(5) =  0.0000000000000000D0
    gp_loc(6) =  0.3242534234038089D0
    gp_loc(7) =  0.6133714327005904D0
    gp_loc(8) =  0.8360311073266358D0
    gp_loc(9) =  0.9681602395076261D0
    gp_w(1) = 0.0812743883615744D0
    gp_w(2) = 0.1806481606948574D0
    gp_w(3) = 0.2606106964029354D0
    gp_w(4) = 0.3123470770400029D0
    gp_w(5) = 0.3302393550012598D0
    gp_w(6) = 0.3123470770400029D0
    gp_w(7) = 0.2606106964029354D0
    gp_w(8) = 0.1806481606948574D0
    gp_w(9) = 0.0812743883615744D0
case (10)
    gp_loc(1)  = -0.9739065285171717D0
    gp_loc(2)  = -0.8650633666889845D0
    gp_loc(3)  = -0.6794095682990244D0
    gp_loc(4)  = -0.4333953941292472D0
    gp_loc(5)  = -0.1488743389816312D0
    gp_loc(6)  =  0.1488743389816312D0
    gp_loc(7)  =  0.4333953941292472D0
    gp_loc(8)  =  0.6794095682990244D0
    gp_loc(9)  =  0.8650633666889845D0
    gp_loc(10) =  0.9739065285171717D0
    gp_w(1)  = 0.0666713443086881D0
    gp_w(2)  = 0.1494513491505806D0
    gp_w(3)  = 0.2190863625159820D0
    gp_w(4)  = 0.2692667193099963D0
    gp_w(5)  = 0.2955242247147529D0
    gp_w(6)  = 0.2955242247147529D0
    gp_w(7)  = 0.2692667193099963D0
    gp_w(8)  = 0.2190863625159820D0
    gp_w(9)  = 0.1494513491505806D0
    gp_w(10) = 0.0666713443086881D0
case (11)
    gp_loc(1)  = -0.9782286581460570D0
    gp_loc(2)  = -0.8870625997680953D0
    gp_loc(3)  = -0.7301520055740494D0
    gp_loc(4)  = -0.5190961292068118D0
    gp_loc(5)  = -0.2695431559523450D0
    gp_loc(6)  =  0.0000000000000000D0
    gp_loc(7)  =  0.2695431559523450D0
    gp_loc(8)  =  0.5190961292068118D0
    gp_loc(9)  =  0.7301520055740494D0
    gp_loc(10) =  0.8870625997680953D0
    gp_loc(11) =  0.9782286581460570D0
    gp_w(1)  = 0.0556685671161737D0
    gp_w(2)  = 0.1255803694649046D0
    gp_w(3)  = 0.1862902109277343D0
    gp_w(4)  = 0.2331937645919905D0
    gp_w(5)  = 0.2628045445102467D0
    gp_w(6)  = 0.2729250867779006D0
    gp_w(7)  = 0.2628045445102467D0
    gp_w(8)  = 0.2331937645919905D0
    gp_w(9)  = 0.1862902109277343D0
    gp_w(10) = 0.1255803694649046D0
    gp_w(11) = 0.0556685671161737D0
case (12)
    gp_loc(1)  = -0.9815606342467193D0
    gp_loc(2)  = -0.9041172563704749D0
    gp_loc(3)  = -0.7699026741943047D0
    gp_loc(4)  = -0.5873179542866175D0
    gp_loc(5)  = -0.3678314989981802D0
    gp_loc(6)  = -0.1252334085114689D0
    gp_loc(7)  =  0.1252334085114689D0
    gp_loc(8)  =  0.3678314989981802D0
    gp_loc(9)  =  0.5873179542866175D0
    gp_loc(10) =  0.7699026741943047D0
    gp_loc(11) =  0.9041172563704749D0
    gp_loc(12) =  0.9815606342467193D0
    gp_w(1)  = 0.0471753363865118D0
    gp_w(2)  = 0.1069393259953184D0
    gp_w(3)  = 0.1600783285433462D0
    gp_w(4)  = 0.2031674267230659D0
    gp_w(5)  = 0.2334925365383548D0
    gp_w(6)  = 0.2491470458134028D0
    gp_w(7)  = 0.2491470458134028D0
    gp_w(8)  = 0.2334925365383548D0
    gp_w(9)  = 0.2031674267230659D0
    gp_w(10) = 0.1600783285433462D0
    gp_w(11) = 0.1069393259953184D0
    gp_w(12) = 0.0471753363865118D0
case (13)
    gp_loc(1)  = -0.9841830547185881D0
    gp_loc(2)  = -0.9175983992229779D0
    gp_loc(3)  = -0.8015780907333099D0
    gp_loc(4)  = -0.6423493394403402D0
    gp_loc(5)  = -0.4484927510364469D0
    gp_loc(6)  = -0.2304583159551348D0
    gp_loc(7)  =  0.0000000000000000D0
    gp_loc(8)  =  0.2304583159551348D0
    gp_loc(9)  =  0.4484927510364469D0
    gp_loc(10) =  0.6423493394403402D0
    gp_loc(11) =  0.8015780907333099D0
    gp_loc(12) =  0.9175983992229779D0
    gp_loc(13) =  0.9841830547185881D0
    gp_w(1)  = 0.0404840047653159D0
    gp_w(2)  = 0.0921214998377285D0
    gp_w(3)  = 0.1388735102197872D0
    gp_w(4)  = 0.1781459807619457D0
    gp_w(5)  = 0.2078160475368885D0
    gp_w(6)  = 0.2262831802628972D0
    gp_w(7)  = 0.2325515532308739D0
    gp_w(8)  = 0.2262831802628972D0
    gp_w(9)  = 0.2078160475368885D0
    gp_w(10) = 0.1781459807619457D0
    gp_w(11) = 0.1388735102197872D0
    gp_w(12) = 0.0921214998377285D0
    gp_w(13) = 0.0404840047653159D0
case (14)
    gp_loc(1)  = -0.9862838086968123D0
    gp_loc(2)  = -0.9284348836635735D0
    gp_loc(3)  = -0.8272013150697650D0
    gp_loc(4)  = -0.6872929048116855D0
    gp_loc(5)  = -0.5152486363581541D0
    gp_loc(6)  = -0.3191123689278897D0
    gp_loc(7)  = -0.1080549487073437D0
    gp_loc(8)  =  0.1080549487073437D0
    gp_loc(9)  =  0.3191123689278897D0
    gp_loc(10) =  0.5152486363581541D0
    gp_loc(11) =  0.6872929048116855D0
    gp_loc(12) =  0.8272013150697650D0
    gp_loc(13) =  0.9284348836635735D0
    gp_loc(14) =  0.9862838086968123D0
    gp_w(1)  = 0.0351194603317519D0
    gp_w(2)  = 0.0801580871597602D0
    gp_w(3)  = 0.1215185706879032D0
    gp_w(4)  = 0.1572031671581935D0
    gp_w(5)  = 0.1855383974779378D0
    gp_w(6)  = 0.2051984637212956D0
    gp_w(7)  = 0.2152638534631578D0
    gp_w(8)  = 0.2152638534631578D0
    gp_w(9)  = 0.2051984637212956D0
    gp_w(10) = 0.1855383974779378D0
    gp_w(11) = 0.1572031671581935D0
    gp_w(12) = 0.1215185706879032D0
    gp_w(13) = 0.0801580871597602D0
    gp_w(14) = 0.0351194603317519D0
case (15)
    gp_loc(1)  = -0.9879925180204854D0
    gp_loc(2)  = -0.9372733924007060D0
    gp_loc(3)  = -0.8482065834104272D0
    gp_loc(4)  = -0.7244177313601701D0
    gp_loc(5)  = -0.5709721726085388D0
    gp_loc(6)  = -0.3941513470775634D0
    gp_loc(7)  = -0.2011940939974345D0
    gp_loc(8)  =  0.0000000000000000D0
    gp_loc(9)  =  0.2011940939974345D0
    gp_loc(10) =  0.3941513470775634D0
    gp_loc(11) =  0.5709721726085388D0
    gp_loc(12) =  0.7244177313601701D0
    gp_loc(13) =  0.8482065834104272D0
    gp_loc(14) =  0.9372733924007060D0
    gp_loc(15) =  0.9879925180204854D0
    gp_w(1)  = 0.0307532419961173D0
    gp_w(2)  = 0.0703660474881081D0
    gp_w(3)  = 0.1071592204671719D0
    gp_w(4)  = 0.1395706779261543D0
    gp_w(5)  = 0.1662692058169939D0
    gp_w(6)  = 0.1861610000155798D0
    gp_w(7)  = 0.1984314853271116D0
    gp_w(8)  = 0.2025782419255613D0
    gp_w(9)  = 0.1984314853271116D0
    gp_w(10) = 0.1861610000155798D0
    gp_w(11) = 0.1662692058169939D0
    gp_w(12) = 0.1395706779261543D0
    gp_w(13) = 0.1071592204671719D0
    gp_w(14) = 0.0703660474881081D0
    gp_w(15) = 0.0307532419961173D0
case default
    print *, '    ERROR-2026051403 :: Unsupported number of Gauss points: ', num_gp
    print *, '                        Supported values: 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15'
    call Warning_Message('S', Keywords_Blank)
end select

!******************************************************    
! Loop over Domain C to perform integral calculation.
!******************************************************
do i_E = 1, number_of_elements
    element_number = Domain_Elements_List(i_E)

    c_NN(1:8) = G_NN(1:8, element_number)
    omega = ZR

    local_pos_of_proj_nodes(1:number_nodes_projected) = 0
    do i = 1, number_nodes_projected
        c_node = nodes_to_be_projected(i)
        do j = 1, 8
            if (G_NN(j, element_number) == c_node) then
                local_pos_of_proj_nodes(i) = j
            endif
        enddo
    enddo

    c_X_NODES(1:8) = G_X_NODES(1:8, element_number)
    c_Y_NODES(1:8) = G_Y_NODES(1:8, element_number)
    c_Z_NODES(1:8) = G_Z_NODES(1:8, element_number)

    !----------------------------------------------------------
    ! Get the updated crack-tip local system from global arrays.
    ! Solid_El_Crs uses current-step data.
    !----------------------------------------------------------
    ref_elem = 0
    !ref_elem = element_number
    
    ! Obtain the reference element number.
    if ((Elem_Type_3D(element_number,crack_number).eq.1 )) then
        ref_elem = element_number
    else
        ! ref_elem = Ele_Num_Tip_Enriched_Node_3D(c_NN(i_N), i_C) !The enriched element number corresponding to the enriched node (reference element number)
        do i_N = 1,8
            if(allocated(Ele_Num_Tip_Enriched_Node_3D(c_NN(i_N))%row)) then
                if(Ele_Num_Tip_Enriched_Node_3D(c_NN(i_N))%row(crack_number) >= 0) then
                    ref_elem=Ele_Num_Tip_Enriched_Node_3D(c_NN(i_N))%row(crack_number)
                    exit
                endif
            endif
        enddo
    endif
    
    if(ref_elem==0) then
        cycle
    endif
          
    temp_vector = Solid_El_Crs(ref_elem, 1:Solid_El_Max_num_Crs)
call Vector_Location_Int_v2( Solid_El_Max_num_Crs, temp_vector, crack_number, c_Cr_Location)
    BaseLine_A   = Solid_El_Tip_BaseLine(ref_elem)%row(c_Cr_Location, 1, 1:3)
    BaseLine_B   = Solid_El_Tip_BaseLine(ref_elem)%row(c_Cr_Location, 2, 1:3)
    BaseLine_Mid = (BaseLine_A + BaseLine_B) / TWO
    BaseLine_x_Vec = Solid_El_Tip_BaseLine_x_Vec(ref_elem)%row(c_Cr_Location, 1:3)
    BaseLine_y_Vec = Solid_El_Tip_BaseLine_y_Vec(ref_elem)%row(c_Cr_Location, 1:3)
    BaseLine_z_Vec = Solid_El_Tip_BaseLine_z_Vec(ref_elem)%row(c_Cr_Location, 1:3)
    c_T_Matrx = Solid_El_Tip_BaseLine_T_Matrix(ref_elem)%row(c_Cr_Location, 1:3, 1:3)

    !----------------------------------------------------------
    ! Get the old crack-tip local system from global OLD arrays.
    !----------------------------------------------------------
    if(.not. allocated(Solid_El_Tip_BaseLine_OLD(ref_elem)%row)) then
        print *, '    ERROR-2026051404 :: Solid_El_Tip_BaseLine_OLD(ref_elem)%row not allocated!'
        call Warning_Message('S', Keywords_Blank)
    endif
    BaseLine_A_OLD   = Solid_El_Tip_BaseLine_OLD(ref_elem)%row(c_Cr_Location, 1, 1:3)
    BaseLine_B_OLD   = Solid_El_Tip_BaseLine_OLD(ref_elem)%row(c_Cr_Location, 2, 1:3)
    BaseLine_Mid_OLD = (BaseLine_A_OLD + BaseLine_B_OLD) / TWO
    if(.not. allocated(Solid_El_Tip_BaseLine_x_Vec_OLD(ref_elem)%row)) then
        print *, '    ERROR-2026051405 :: Solid_El_Tip_BaseLine_x_Vec_OLD(ref_elem)%row not allocated!'
        call Warning_Message('S', Keywords_Blank)
    endif
    BaseLine_x_Vec_OLD = Solid_El_Tip_BaseLine_x_Vec_OLD(ref_elem)%row(c_Cr_Location, 1:3)
    BaseLine_y_Vec_OLD = Solid_El_Tip_BaseLine_y_Vec_OLD(ref_elem)%row(c_Cr_Location, 1:3)
    BaseLine_z_Vec_OLD = Solid_El_Tip_BaseLine_z_Vec_OLD(ref_elem)%row(c_Cr_Location, 1:3)
    if(.not. allocated(Solid_El_Tip_BaseLine_T_Matrix_OLD(ref_elem)%row)) then
        print *, '    ERROR-2026051405 :: Solid_El_Tip_BaseLine_T_Matrix_OLD(ref_elem)%row not allocated!'
        call Warning_Message('S', Keywords_Blank)
    endif
    c_T_Matrx_OLD = Solid_El_Tip_BaseLine_T_Matrix_OLD(ref_elem)%row(c_Cr_Location, 1:3, 1:3)
    
!    BaseLine_A_OLD   = Solid_El_Tip_BaseLine(ref_elem)%row(c_Cr_Location, 1, 1:3)
!    BaseLine_B_OLD   = Solid_El_Tip_BaseLine(ref_elem)%row(c_Cr_Location, 2, 1:3)
!    BaseLine_Mid_OLD = (BaseLine_A_OLD + BaseLine_B_OLD) / TWO
!    BaseLine_x_Vec_OLD = Solid_El_Tip_BaseLine_x_Vec(ref_elem)%row(c_Cr_Location, 1:3)
!    BaseLine_y_Vec_OLD = Solid_El_Tip_BaseLine_y_Vec(ref_elem)%row(c_Cr_Location, 1:3)
!    BaseLine_z_Vec_OLD = Solid_El_Tip_BaseLine_z_Vec(ref_elem)%row(c_Cr_Location, 1:3)
!    c_T_Matrx_OLD = Solid_El_Tip_BaseLine_T_Matrix(ref_elem)%row(c_Cr_Location, 1:3, 1:3)

    ! Material number/type uses current step directly.
    mat_num    = Elem_Mat(ref_elem)
    c_mat_type = Material_Type(mat_num)

    do g1 = 1, num_gp
        do g2 = 1, num_gp
            do g3 = 1, num_gp
                kesi   = gp_loc(g1)
                yita   = gp_loc(g2)
                zeta   = gp_loc(g3)
                weight = gp_w(g1) * gp_w(g2) * gp_w(g3)

call Cal_N_dNdkesi_J_detJ_3D( kesi, yita, zeta, c_X_NODES, c_Y_NODES, c_Z_NODES, detJ, Jmat, Nmat, dNdkesi)

                N8(1:8) = Nmat(1, 1:24:3)

                Global_coor_Gauss(1) = dot_product(N8, c_X_NODES)
                Global_coor_Gauss(2) = dot_product(N8, c_Y_NODES)
                Global_coor_Gauss(3) = dot_product(N8, c_Z_NODES)

                !--------------------------------------
                ! New tip function at Gauss point
                !--------------------------------------
                Temp_Point = Global_coor_Gauss - BaseLine_Mid
                Gauss_Coor_Local_New = matmul(c_T_Matrx, Temp_Point)
                r_Gauss_new     = sqrt(Gauss_Coor_Local_New(1)**2 + Gauss_Coor_Local_New(2)**2)
                theta_Gauss_new = atan2(Gauss_Coor_Local_New(2), Gauss_Coor_Local_New(1))
call Cal_F_dFdx_dFdy_dFdz_3D( r_Gauss_new, theta_Gauss_new, c_T_Matrx, c_mat_type, &
F_Gauss_new, dFdx_tip, dFdy_tip, dFdz_tip)

                !--------------------------------------
                ! Old tip function at Gauss point
                !--------------------------------------
                Temp_Point = Global_coor_Gauss - BaseLine_Mid_OLD
                Gauss_Coor_Local_Old = matmul(c_T_Matrx_OLD, Temp_Point)
                r_Gauss_old     = sqrt(Gauss_Coor_Local_Old(1)**2 + Gauss_Coor_Local_Old(2)**2)
                theta_Gauss_old = atan2(Gauss_Coor_Local_Old(2), Gauss_Coor_Local_Old(1))
call Cal_F_dFdx_dFdy_dFdz_3D( r_Gauss_old, theta_Gauss_old, c_T_Matrx_OLD, c_mat_type, &
F_Gauss_old, dFdx_tip, dFdy_tip, dFdz_tip)

                do i = 1, number_nodes_projected
                    c_node = nodes_to_be_projected(i)
                    loc_i = local_pos_of_proj_nodes(i)
                    if (loc_i == 0) cycle

                    !--------------------------------------
                    ! New tip basis contribution for node i
                    !--------------------------------------
                    Node_Point(1) = c_X_NODES(loc_i)
                    Node_Point(2) = c_Y_NODES(loc_i)
                    Node_Point(3) = c_Z_NODES(loc_i)
                    Temp_Point = Node_Point - BaseLine_Mid
                    Node_Coor_Local_New = matmul(c_T_Matrx, Temp_Point)
                    r_Node_new     = sqrt(Node_Coor_Local_New(1)**2 + Node_Coor_Local_New(2)**2)
                    theta_Node_new = atan2(Node_Coor_Local_New(2), Node_Coor_Local_New(1))
                    call Cal_F(r_Node_new, theta_Node_new, omega, c_mat_type, F_Node_new_i)
                    phi_tip_new_i = F_Gauss_new(1) - F_Node_new_i(1)

                    do j = 1, number_nodes_projected
                        node_j = nodes_to_be_projected(j)
                        loc_j = local_pos_of_proj_nodes(j)
                        if (loc_j == 0) cycle

                        !--------------------------------------
                        ! New tip basis contribution for node j
                        !--------------------------------------
                        Node_Point(1) = c_X_NODES(loc_j)
                        Node_Point(2) = c_Y_NODES(loc_j)
                        Node_Point(3) = c_Z_NODES(loc_j)
                        Temp_Point = Node_Point - BaseLine_Mid
                        Node_Coor_Local_New = matmul(c_T_Matrx, Temp_Point)
                        r_Node_new     = sqrt(Node_Coor_Local_New(1)**2 + Node_Coor_Local_New(2)**2)
                        theta_Node_new = atan2(Node_Coor_Local_New(2), Node_Coor_Local_New(1))
                        call Cal_F(r_Node_new, theta_Node_new, omega, c_mat_type, F_Node_new_j)
                        phi_tip_new_j = F_Gauss_new(1) - F_Node_new_j(1)

                        !--------------------------------------
                        ! Old tip basis contribution for node j
                        !--------------------------------------
                        Temp_Point = Node_Point - BaseLine_Mid_OLD
                        Node_Coor_Local_Old = matmul(c_T_Matrx_OLD, Temp_Point)
                        r_Node_old     = sqrt(Node_Coor_Local_Old(1)**2 + Node_Coor_Local_Old(2)**2)
                        theta_Node_old = atan2(Node_Coor_Local_Old(2), Node_Coor_Local_Old(1))
                        call Cal_F(r_Node_old, theta_Node_old, omega, c_mat_type, F_Node_old_j)
                        phi_tip_old_j = F_Gauss_old(1) - F_Node_old_j(1)

scalar_M = N8(loc_i) * N8(loc_j) * phi_tip_new_i * phi_tip_new_j * detJ * weight

scalar_P = N8(loc_i) * N8(loc_j) * phi_tip_new_i * phi_tip_old_j * detJ * weight

                        M_tt((i-1)*3 + 1, (j-1)*3 + 1) = M_tt((i-1)*3 + 1, (j-1)*3 + 1) + scalar_M
                        M_tt((i-1)*3 + 2, (j-1)*3 + 2) = M_tt((i-1)*3 + 2, (j-1)*3 + 2) + scalar_M
                        M_tt((i-1)*3 + 3, (j-1)*3 + 3) = M_tt((i-1)*3 + 3, (j-1)*3 + 3) + scalar_M

                        P_tt((i-1)*3 + 1, (j-1)*3 + 1) = P_tt((i-1)*3 + 1, (j-1)*3 + 1) + scalar_P
                        P_tt((i-1)*3 + 2, (j-1)*3 + 2) = P_tt((i-1)*3 + 2, (j-1)*3 + 2) + scalar_P
                        P_tt((i-1)*3 + 3, (j-1)*3 + 3) = P_tt((i-1)*3 + 3, (j-1)*3 + 3) + scalar_P
                    end do
                end do

            end do
        end do
    end do
enddo

!**************************
! Get right hand side.
!**************************
rhs_vec = matmul(P_tt, b_old)

!*************************************
! Calculate the determinant (optional).
!*************************************
call Matrix_Det(number_nodes_projected*3, M_tt, det_Mtt)
print *, '    Det(M_tt)=', det_Mtt

!***********************************************
! Calculate the number of conditions (optional).
!***********************************************
call Matrix_Condition_Number_v2(1,number_nodes_projected*3,M_tt,Condition_Num_Norm)
print *, '    Condition number (Norm1) =', Condition_Num_Norm
call Matrix_Condition_Number_v2(2,number_nodes_projected*3,M_tt,Condition_Num_Norm)
print *, '    Condition number (Norm2) =', Condition_Num_Norm

!******************************
! Get the projected result.
!******************************
call Matrix_Solve_LSOE(0,1,6,M_tt,rhs_vec,b_new,number_nodes_projected*3)

!**************************************
! Obtain the function return value.
!**************************************
nodes_new_enriched_value(1:number_nodes_projected*3) = b_new(1:number_nodes_projected*3)

print *,"    -------------------------------------"
print *,"    DOFs of projected tip enriched nodes:"
print *,"    -------------------------------------"
do i = 1, number_nodes_projected
write(*,'("     Node=", I6, ", X=", F12.6, ", Y=", F12.6, ", Z=", F12.6)') nodes_to_be_projected(i), &
nodes_new_enriched_value((i-1)*3 + 1), nodes_new_enriched_value((i-1)*3 + 2), nodes_new_enriched_value((i-1)*3 + 3)
end do

!********************
! Clear the data.
!********************
deallocate(b_old)
deallocate(b_new)
deallocate(rhs_vec)
deallocate(M_tt)
deallocate(P_tt)
deallocate(M_tt_inv)
deallocate(Domain_Elements_List)
deallocate(gp_loc)
deallocate(gp_w)

return
end subroutine Local_L2_projection_for_intra_element