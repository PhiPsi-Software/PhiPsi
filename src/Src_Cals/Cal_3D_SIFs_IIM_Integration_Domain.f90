!-----------------------------------------------------------
! Brief: Build IIM integration domain around a crack-front vertex.
!
! Parameters:
!   Input:  Tip_Point  - coordinates of the crack-front vertex
!           e1, e2, e3 - local basis (e3 along the front)
!           R_in, R_out- inner/outer radii of the domain
!           L_f        - half-length along the front (0 = radial)
!   Output: Num_Elems  - number of elements in the domain
!           Elem_List  - list of element IDs
!           q_vals     - nodal q-function values per element
!
! Notes:   When L_f = 0 a radial spherical domain is used; otherwise
!          a cylindrical (front-aware) domain is constructed.
!-----------------------------------------------------------

subroutine Cal_3D_SIFs_IIM_Integration_Domain(Tip_Point, e1, e2, e3, R_in, R_out, L_f, Num_Elems, Elem_List, q_vals)
!-------------------------------------------------------------------------------
! PURPOSE: Build integration domain around crack front vertex
!
! SHAPE: 
!   - If L_f = 0: Simple radial sphere (for point-wise calculation)
!   - If L_f > 0: Cylindrical domain (full 3D method)
!
! INPUT:
!   Tip_Point(3)  - Crack front vertex coordinates
!   e1, e2, e3    - Local coordinate basis vectors
!   R_in, R_out   - Inner and outer radii
!   L_f           - Half-length along crack front (0 = radial only)
!
! OUTPUT:
!   Num_Elems     - Number of elements in domain
!   Elem_List(:)  - List of element IDs
!   q_vals(:,:)   - q-function values at element nodes
!----------------------------
! Rules for assigning the q-function values
!
! (1) Value in the radial direction
! For any point in space, first calculate its radial distance r to the crack front vertex 
!  (Note: This is the distance perpendicular to the crack front tangent):
!
!
! Special case: L_front = 0
! When L_front = 0 is set, it degenerates to a simple spherical domain (actually the application of
! a 2D radial domain in 3D):
!
! Only considers the Euclidean distance to the vertex
! Does not distinguish between axial and radial directions
! Suitable for cases where the crack front is locally approximated as a straight line
! Practical engineering significance
! R_inner: Typically set to 0.5 times the mesh size, ensuring sufficient sampling points inside
! R_outer: Typically set to 3 times the mesh size, ensuring capture of the main influence range of
! the stress field
! L_front:
!   When set to 0, assumes the stress field does not vary along the front (suitable for long
!   straight cracks)
!   When set to 2 times the mesh size, considers the influence of front curvature (suitable for
!   complex 3D cracks)
!-------------------------------------------------------------------------------
use Global_Float_Type
use Global_Model

implicit none
real(kind=FT), intent(in) :: Tip_Point(3), e1(3), e2(3), e3(3)
real(kind=FT), intent(in) :: R_in, R_out, L_f
integer, intent(out) :: Num_Elems
integer, allocatable, intent(out) :: Elem_List(:)
real(kind=FT), allocatable, intent(out) :: q_vals(:,:)

integer :: i_Elem, i_Node, Node_ID
real(kind=FT) :: Node_Coords(3), Vec_to_Node(3)
real(kind=FT) :: r_radial, s_along_front, dist
real(kind=FT) :: q_node
logical :: Elem_in_Domain
integer, allocatable :: Temp_List(:)
integer :: Count

allocate(Temp_List(num_Elem))
allocate(q_vals(num_Elem, 8))
q_vals = ZR
Count = 0

do i_Elem = 1, num_Elem
    
    Elem_in_Domain = .false.
    
    do i_Node = 1, 8
        Node_ID = G_NN(i_Node, i_Elem)
        Node_Coords(1) = Coor(Node_ID, 1)
        Node_Coords(2) = Coor(Node_ID, 2)
        Node_Coords(3) = Coor(Node_ID, 3)
        
        Vec_to_Node = Node_Coords - Tip_Point
        
        if (L_f > 1.0D-12) then
            s_along_front = dot_product(Vec_to_Node, e3)
            r_radial = norm2(Vec_to_Node - s_along_front * e3)
            if (r_radial <= R_out .AND. abs(s_along_front) <= 2.0D0*L_f) then
                Elem_in_Domain = .true.
                call Cal_3D_SIFs_IIM_q_3D(r_radial, s_along_front, R_in, R_out, L_f, q_node)
                q_vals(i_Elem, i_Node) = q_node
            endif
            
        else
            dist = norm2(Vec_to_Node)
            
            if (dist <= R_out) then
                Elem_in_Domain = .true.
                
                if (dist <= R_in) then
                    q_node = ONE
                else
                    q_node = (R_out - dist) / (R_out - R_in)
                endif
                q_vals(i_Elem, i_Node) = q_node
            endif
        endif
        
    enddo
    
    if(Elem_in_Domain) then
        Count = Count + 1
        Temp_List(Count) = i_Elem
    endif
enddo

Num_Elems = Count
allocate(Elem_List(Num_Elems))
Elem_List(1:Num_Elems) = Temp_List(1:Num_Elems)

deallocate(Temp_List)
    
    
end subroutine Cal_3D_SIFs_IIM_Integration_Domain