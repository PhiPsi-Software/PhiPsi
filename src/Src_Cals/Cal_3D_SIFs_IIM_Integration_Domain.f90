
subroutine Cal_3D_SIFs_IIM_Integration_Domain(Tip_Point, e1, e2, e3, R_in, R_out, L_f, &
                                    Num_Elems, Elem_List, q_vals)
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



