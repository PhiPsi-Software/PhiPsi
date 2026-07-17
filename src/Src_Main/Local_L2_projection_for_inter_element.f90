!-----------------------------------------------------------
! Brief: 3D history-preserving local L2 projection that converts
!        obsolete crack-tip enriched DOFs to Heaviside enriched
!        DOFs when a crack tip crosses an element boundary.
!
! Parameters:
!   Input:  isub                   - sub-step index
!           crack_number           - propagating crack index
!           num_gp                 - Gauss points per direction
!           c_DOFs_Value           - global DOF vector
!           number_nodes_projected - size of the projection set
!           nodes_to_be_projected  - node IDs to be re-enriched
!   Output: nodes_new_enriched_value - new Heaviside DOFs (xyz)
!
! Notes:   Solves M_HH * a_new = P_Ht * b_old over the union of
!          elements owning the projected nodes; only the first tip
!          branch function is used.
!-----------------------------------------------------------

subroutine Local_L2_projection_for_inter_element( isub, crack_number,num_gp, c_DOFs_Value, number_nodes_projected, &
nodes_to_be_projected, nodes_new_enriched_value)
!-----------------------------------------------------------------------
! Subroutine Local_L2_projection_for_inter_element
!
! Purpose:
!   Perform a history-preserving L2 projection of crack-tip enriched DOFs
!   onto Heaviside enriched DOFs during inter-element crack-tip crossing.
!   This routine implements the enrichment-type transition described in
!   the paper, where obsolete crack-tip enriched nodes are converted to 
!   Heaviside enriched nodes before the crack tip enters the adjacent 
!   element.
!
!   The projection solves:  M_HH * a_new = P_Ht * b_old, where
!     - b_old : vector of old crack-tip enriched displacement/velocity/
!               acceleration DOFs (collected from c_DOFs_Value)
!     - a_new : resulting Heaviside enriched DOFs (output)
!   The matrices M_HH and P_Ht are assembled by numerical integration
!   over the support domain (union of elements containing the projected
!   nodes). The same operator is applied to displacement, velocity and
!   acceleration to preserve kinematic consistency.
!
! Assumptions:
!   - Only a single crack-tip enrichment function F (the first branch
!     function, typically sqrt(r)*sin(theta/2)) is used.
!   - The crack is represented in 3D; nodes are enriched with a single
!     tip enrichment DOF set per node.
!   - The current crack-tip position and orientation are stored in global
!     data structures (Solid_El_Tip_BaseLine_*).
!   - Heaviside enrichment is defined using the sign of the signed
!     distance to the crack surface.
!
! Usage example:
!   crack_number   = 1
!   number_nodes_projected   = 4
!   allocate(nodes_to_be_projected(number_nodes_projected))
!   allocate(nodes_new_enriched_value(number_nodes_projected*3))
!   nodes_to_be_projected(1:number_nodes_projected) = [265,625,626,266]
!   nodes_new_enriched_value(1:number_nodes_projected*3) = 0.0_FT
!   num_gp = 14   ! Gauss integration order (2..15 supported)
!   call Local_L2_projection_for_inter_element(isub, crack_number, num_gp, &
!          DISP(1:Total_FD), number_nodes_projected, &
!          nodes_to_be_projected(1:number_nodes_projected), &
!          nodes_new_enriched_value(1:number_nodes_projected*3))
!
! Arguments:
!   isub                       (in)  integer : sub-stepping index (for dynamic analysis)
!   crack_number               (in)  integer : which crack (fracture) is propagating
!   num_gp                     (in)  integer : number of Gauss points per direction
!                                              (supported: 2,3,...,15)
!   c_DOFs_Value               (in)  real(FT) : global displacement (or velocity/acceleration)
!                                              vector; contains both standard and enriched DOFs.
!   number_nodes_projected     (in)  integer : number of nodes to be converted from
!                                              crack-tip enrichment to Heaviside enrichment.
!   nodes_to_be_projected      (in)  integer(:) : list of node IDs to be converted.
!   nodes_new_enriched_value   (out) real(FT) : output array of length
!                                              number_nodes_projected * 3,
!                                              containing the new Heaviside enriched
!                                              DOFs (x, y, z components interleaved).
!
! Global data dependencies (modules):
!   Global_Float_Type, Global_Common, Global_Filename, Global_Model,
!   Global_Elem_Area_Vol, Global_Crack_Common,Global_Crack_3D, 
!   Global_INTERFACE_D3_Get_Signed_Dis_to_Crack_Mesh,
!   Global_INTERFACE_D3_Get_Signed_Dis_to_Crack_Mesh_for_InPlane,
!   Global_Inter_Cal_N_dNdkesi_J_detJ_3D
!
! Notes:
!   - The projection domain is the union of all elements that contain any
!     of the nodes in nodes_to_be_projected.
!   - The function checks that each node in nodes_to_be_projected is
!     currently a crack-tip enriched node (Enriched_Node_Type_3D == 1);
!     otherwise an error is printed and the program may stop.
!   - Numerical integration uses Gauss-Legendre quadrature; the number of
!     Gauss points can be adjusted via num_gp (higher values increase
!     accuracy but also computational cost).
!   - The linear system M_HH * a_new = rhs is solved by the linear solver
!     (Matrix_Solve_LSOE). The determinant and condition number of M_HH
!     are printed for debugging.
!   - The same projection operator is applied to displacement, velocity,
!     and acceleration fields by calling this routine with the respective
!     global DOF vector.
!
! Implementation references:
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
use Global_INTERFACE_D3_Get_Signed_Dis_to_Crack_Mesh
use Global_INTERFACE_D3_Get_Signed_Dis_to_Crack_Mesh_for_InPlane
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
integer :: info_inv
integer :: c_NN(8)
integer :: temp_vector(Solid_El_Max_num_Crs)
real(kind=FT), allocatable :: b_old(:), a_new(:), rhs_vec(:)
real(kind=FT), allocatable :: M_HH(:,:), P_Ht(:,:), M_HH_inv(:,:)
real(kind=FT) :: c_X_NODES(8), c_Y_NODES(8), c_Z_NODES(8)
real(kind=FT), allocatable  :: gp_loc(:), gp_w(:)
real(kind=FT) :: kesi, yita, zeta, weight
real(kind=FT) :: detJ
real(kind=FT) :: Jmat(3,3), Nmat(3,24), dNdkesi(8,3)
real(kind=FT) :: N8(8)
real(kind=FT) :: scalar_M, scalar_P
real(kind=FT) :: H_gp, H_i, H_j
real(kind=FT) :: phi_H_i, phi_H_j, phi_tip_j
real(kind=FT) :: Global_coor_Gauss(3)
real(kind=FT) :: c_Distance, c_Signed_Dis_v2
real(kind=FT) :: c_PER_Node_to_FS(3), c_n_Vector(3)
logical :: c_Yes_Node_PER_in_FS, Yes_Found_Min_Signed_Dis
real(kind=FT) :: dFdx_tip(4)
real(kind=FT) :: dFdy_tip(4)
real(kind=FT) :: dFdz_tip(4)
real(kind=FT) :: F_Gauss(4)
real(kind=FT) :: F_Node(4)
real(kind=FT) :: BaseLine_A(3), BaseLine_B(3), BaseLine_Mid(3)
real(kind=FT) :: BaseLine_x_Vec(3), BaseLine_y_Vec(3), BaseLine_z_Vec(3)
real(kind=FT) :: c_T_Matrx(3,3)
real(kind=FT) :: Node_Point(3), Node_Coor_Local(3), Gauss_Coor_Local(3)
real(kind=FT) :: Temp_Point(3)
real(kind=FT) :: r_Node, theta_Node, r_Gauss, theta_Gauss, omega
real(kind=FT) :: Check_Ball_R
real(kind=FT) :: det_MHH,Condition_Num_Norm
integer i_E,element_number
integer local_pos_of_proj_nodes(number_nodes_projected)
logical :: found
integer Temp_Domain_Elements(1000),Uniqued_Temp_Domain_Elements(1000),i_E_count
integer Uniqued_n
integer number_of_elements
integer,ALLOCATABLE::Domain_Elements_List(:)

print *,"    Perform L2 projection of crack-tip enriched DOFs onto Heaviside enriched DOFs..."

!****************************
! Check the projected nodes.
!****************************
do i = 1, number_nodes_projected
    c_node = nodes_to_be_projected(i)
    if (Enriched_Node_Type_3D(c_node, crack_number) /= 1) then
        print *, '    ERROR-2026051301 :: projected node is not tip enriched node.'
        print *, '    In Local_L2_projection_for_inter_element.'
        print *, '    Node number: ', c_node
        call Warning_Message('S', Keywords_Blank)
    end if
end do

!****************************
! Initialization.
!****************************
nodes_new_enriched_value(1:number_nodes_projected*3) = ZR
allocate(b_old(number_nodes_projected*3))
allocate(a_new(number_nodes_projected*3))
allocate(rhs_vec(number_nodes_projected*3))
allocate(M_HH(number_nodes_projected*3, number_nodes_projected*3))
allocate(P_Ht(number_nodes_projected*3, number_nodes_projected*3))
allocate(M_HH_inv(number_nodes_projected*3, number_nodes_projected*3))
b_old    = ZR
a_new    = ZR
rhs_vec  = ZR
M_HH     = ZR
P_Ht     = ZR
M_HH_inv = ZR

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
!Remove duplicate elements from the array, generate a new array.
call Vector_Unique_Int(1000,i_E_count,Temp_Domain_Elements(1:1000),Uniqued_Temp_Domain_Elements(1:1000),number_of_elements)   
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
    ! 2-point Gauss-Legendre integration on [-1,1]
    gp_loc(1) = -0.5773502691896258D0
    gp_loc(2) =  0.5773502691896258D0
    gp_w(1) = 1.0000000000000000D0
    gp_w(2) = 1.0000000000000000D0
case (3)
    ! 3-point Gauss-Legendre integration on [-1,1]
    gp_loc(1) = -0.7745966692414834D0
    gp_loc(2) =  0.0000000000000000D0
    gp_loc(3) =  0.7745966692414834D0
    gp_w(1) = 0.5555555555555556D0
    gp_w(2) = 0.8888888888888889D0
    gp_w(3) = 0.5555555555555556D0
case (4)
    ! 4-point Gauss-Legendre integration on [-1,1]
    gp_loc(1) = -0.8611363115940526D0
    gp_loc(2) = -0.3399810435848563D0
    gp_loc(3) =  0.3399810435848563D0
    gp_loc(4) =  0.8611363115940526D0
    gp_w(1) = 0.3478548451374538D0
    gp_w(2) = 0.6521451548625461D0
    gp_w(3) = 0.6521451548625461D0
    gp_w(4) = 0.3478548451374538D0
case (5)
    ! 5-point Gauss-Legendre integration on [-1,1]
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
    ! 6-point Gauss-Legendre integration on [-1,1]
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
    ! 7-point Gauss-Legendre integration on [-1,1]
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
    ! 8-point Gauss-Legendre integration on [-1,1]
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
    ! 9-point Gauss-Legendre integration on [-1,1]
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
    ! 10-point Gauss-Legendre integration on [-1,1]
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
    ! 11-point Gauss-Legendre integration on [-1,1]
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
    ! 12-point Gauss-Legendre integration on [-1,1]
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
    ! 13-point Gauss-Legendre integration on [-1,1]
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
    ! 14-point Gauss-Legendre integration on [-1,1]
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
    ! 15-point Gauss-Legendre integration on [-1,1]
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
    print *, '    ERROR-2026051401 :: Unsupported number of Gauss points: ', num_gp
    print *, '                        Supported values: 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15'
    call Warning_Message('S', Keywords_Blank)
end select


!******************************************************    
!Loop over Domain C to perform integral calculation.
!******************************************************
do i_E = 1,number_of_elements
    element_number = Domain_Elements_List(i_E)
    ! Get the current element nodes and coordinates
    c_NN(1:8)      = G_NN(1:8, element_number)
    omega = ZR
    
    ! Get local number of nodes.
    local_pos_of_proj_nodes(1:number_nodes_projected) = 0
    do i = 1, number_nodes_projected
        c_node = nodes_to_be_projected(i)
        do j = 1, 8
            if (G_NN(j, element_number) == c_node) then
                local_pos_of_proj_nodes(i) = j
            endif
        enddo
    enddo
    
    
    ! Get the coodinates of nodes of element.
    c_X_NODES(1:8) = G_X_NODES(1:8, element_number)
    c_Y_NODES(1:8) = G_Y_NODES(1:8, element_number)
    c_Z_NODES(1:8) = G_Z_NODES(1:8, element_number)
        
    !Loop over integration points.
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
                ! Heaviside function at Gauss point
                Check_Ball_R = 3.0D0 * Ave_Elem_L
                if (Key_InPlane_Growth == 0) then
call D3_Get_Signed_Dis_to_Crack_Mesh( Global_coor_Gauss, crack_number, Check_Ball_R, c_Distance, c_Signed_Dis_v2, &
c_Yes_Node_PER_in_FS, c_PER_Node_to_FS, Yes_Found_Min_Signed_Dis, c_n_Vector)
                else if (Key_InPlane_Growth == 1) then
call D3_Get_Signed_Dis_to_Crack_Mesh_for_InPlane_Growth( Global_coor_Gauss, crack_number, Check_Ball_R, c_Distance, &
c_Yes_Node_PER_in_FS, c_PER_Node_to_FS, Yes_Found_Min_Signed_Dis, c_n_Vector)
                end if

                if (.not. Yes_Found_Min_Signed_Dis) cycle

                call Cal_Sign(c_Distance, H_gp)

                do i = 1, number_nodes_projected
                    c_node = nodes_to_be_projected(i)
                    loc_i = local_pos_of_proj_nodes(i)
                    if (loc_i == 0) cycle
                    call Cal_Sign(Dis_Node_to_FS(c_node)%row(crack_number), H_i)
                    phi_H_i = H_gp - H_i
                    do j = 1, number_nodes_projected
                        node_j = nodes_to_be_projected(j)
                        loc_j = local_pos_of_proj_nodes(j)
                        if (loc_j == 0) cycle
                        ! Heaviside basis contribution
                        call Cal_Sign(Dis_Node_to_FS(node_j)%row(crack_number), H_j)
                        phi_H_j = H_gp - H_j
                        scalar_M = N8(loc_i) * N8(loc_j) * phi_H_i * phi_H_j * detJ * weight
                        ! Tip basis contribution for node_j
                        if (Elem_Type_3D(element_number, crack_number) == 1) then
                            ref_elem = element_number
                        else
                            ref_elem = Ele_Num_Tip_Enriched_Node_3D(node_j)%row(crack_number)
                        end if
                        temp_vector = Solid_El_Crs(ref_elem, 1:Solid_El_Max_num_Crs)
call Vector_Location_Int_v2( Solid_El_Max_num_Crs, temp_vector, crack_number, c_Cr_Location)
                        BaseLine_A   = Solid_El_Tip_BaseLine(ref_elem)%row(c_Cr_Location, 1, 1:3)
                        BaseLine_B   = Solid_El_Tip_BaseLine(ref_elem)%row(c_Cr_Location, 2, 1:3)
                        BaseLine_Mid = (BaseLine_A + BaseLine_B) / TWO
                        BaseLine_x_Vec = Solid_El_Tip_BaseLine_x_Vec(ref_elem)%row(c_Cr_Location, 1:3)
                        BaseLine_y_Vec = Solid_El_Tip_BaseLine_y_Vec(ref_elem)%row(c_Cr_Location, 1:3)
                        BaseLine_z_Vec = Solid_El_Tip_BaseLine_z_Vec(ref_elem)%row(c_Cr_Location, 1:3)
                        c_T_Matrx = Solid_El_Tip_BaseLine_T_Matrix(ref_elem)%row(c_Cr_Location, 1:3, 1:3)
                        !c_ThetaX  = Solid_El_Tip_BaseLine_T_theta(ref_elem)%row(c_Cr_Location, 1)
                        !c_ThetaY  = Solid_El_Tip_BaseLine_T_theta(ref_elem)%row(c_Cr_Location, 2)
                        !c_ThetaZ  = Solid_El_Tip_BaseLine_T_theta(ref_elem)%row(c_Cr_Location, 3)
                        mat_num    = Elem_Mat(ref_elem)
                        c_mat_type = Material_Type(mat_num)
                        ! Gauss point local coordinates in tip system
                        Temp_Point = Global_coor_Gauss - BaseLine_Mid
                        Gauss_Coor_Local = matmul(c_T_Matrx, Temp_Point)
                        r_Gauss     = sqrt(Gauss_Coor_Local(1)**2 + Gauss_Coor_Local(2)**2)
                        theta_Gauss = atan2(Gauss_Coor_Local(2), Gauss_Coor_Local(1))
                        ! Node local coordinates in tip system
                        Node_Point(1) = c_X_NODES(loc_j)
                        Node_Point(2) = c_Y_NODES(loc_j)
                        Node_Point(3) = c_Z_NODES(loc_j)
                        Temp_Point = Node_Point - BaseLine_Mid
                        Node_Coor_Local = matmul(c_T_Matrx, Temp_Point)
                        r_Node     = sqrt(Node_Coor_Local(1)**2 + Node_Coor_Local(2)**2)
                        theta_Node = atan2(Node_Coor_Local(2), Node_Coor_Local(1))
call Cal_F_dFdx_dFdy_dFdz_3D( r_Gauss, theta_Gauss, c_T_Matrx, c_mat_type, F_Gauss, dFdx_tip, dFdy_tip, dFdz_tip)
                        call Cal_F(r_Node, theta_Node, omega, c_mat_type, F_Node)
                        phi_tip_j = F_Gauss(1) - F_Node(1)
                        scalar_P = N8(loc_i) * N8(loc_j) * phi_H_i * phi_tip_j * detJ * weight
                        ! Assemble M_HH
                        M_HH((i-1)*3 + 1, (j-1)*3 + 1) = M_HH((i-1)*3 + 1, (j-1)*3 + 1) + scalar_M
                        M_HH((i-1)*3 + 2, (j-1)*3 + 2) = M_HH((i-1)*3 + 2, (j-1)*3 + 2) + scalar_M
                        M_HH((i-1)*3 + 3, (j-1)*3 + 3) = M_HH((i-1)*3 + 3, (j-1)*3 + 3) + scalar_M
                        ! Assemble P_Ht
                        P_Ht((i-1)*3 + 1, (j-1)*3 + 1) = P_Ht((i-1)*3 + 1, (j-1)*3 + 1) + scalar_P
                        P_Ht((i-1)*3 + 2, (j-1)*3 + 2) = P_Ht((i-1)*3 + 2, (j-1)*3 + 2) + scalar_P
                        P_Ht((i-1)*3 + 3, (j-1)*3 + 3) = P_Ht((i-1)*3 + 3, (j-1)*3 + 3) + scalar_P
                    end do
                end do

            end do
        end do
    end do

enddo

!do i = 1, number_nodes_projected*3
!end do
!
!
!do i = 1, number_nodes_projected*3
!end do
!
!
!
!


!**************************
! Get right hand side.
!**************************
rhs_vec = matmul(P_Ht, b_old)


!*************************************
! Calculate the determinant (optional).
!*************************************
call Matrix_Det(number_nodes_projected*3, M_HH, det_MHH)
print *, '    Det(M_HH)=', det_MHH

!***********************************************
! Calculate the number of conditions (optional).
!***********************************************
call Matrix_Condition_Number_v2(1,number_nodes_projected*3,M_HH,Condition_Num_Norm)
print *, '    Condition number (Norm1) =', Condition_Num_Norm
call Matrix_Condition_Number_v2(2,number_nodes_projected*3,M_HH,Condition_Num_Norm)
print *, '    Condition number (Norm2) =', Condition_Num_Norm

!******************************
! Get the projected result.
!******************************
!OPTION 1: Solve the linear system.
call Matrix_Solve_LSOE(0,1,6,M_HH,rhs_vec,a_new,number_nodes_projected*3) 
!OPTION 2: Find the inverse matrix and then obtain the result.
!call Matrix_Inverse(M_HH, M_HH_inv, number_nodes_projected*3)
!a_new = matmul(M_HH_inv, rhs_vec)

!**************************************
! Obtain the function return value.
!**************************************
nodes_new_enriched_value(1:number_nodes_projected*3) = a_new(1:number_nodes_projected*3)
print *,"    -------------------------------------"
print *,"    DOFs of projected Heaviside nodes:   "
print *,"    -------------------------------------"
do i = 1, number_nodes_projected
write(*,'("     Node=", I6, ", X=", F12.6, ", Y=", F12.6, ", Z=", F12.6)') nodes_to_be_projected(i), &
nodes_new_enriched_value((i-1)*3 + 1), nodes_new_enriched_value((i-1)*3 + 2), nodes_new_enriched_value((i-1)*3 + 3)
end do

!********************
! Clear the data.
!********************
deallocate(b_old)
deallocate(a_new)
deallocate(rhs_vec)
deallocate(M_HH)
deallocate(P_Ht)
deallocate(M_HH_inv)
deallocate(Domain_Elements_List)

return
end subroutine Local_L2_projection_for_inter_element