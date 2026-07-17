!-----------------------------------------------------------
! Brief: Extract the outer boundary of a 3D crack surface mesh.
!
! Parameters:
!   Input:  i_C                  - crack index
!           num_of_ele           - triangle count
!           num_of_nodes         - node count
!           c_Crack3D_Meshed_Ele - triangle node indices (num_of_ele x 3)
!
! Notes:   Edges appearing once in the triangle list form the outline;
!          chain-sorted end-to-end and stored in
!          Crack3D_Meshed_Outline(i_C)%row with count in
!          Crack3D_Meshed_Outline_num(i_C).
!-----------------------------------------------------------

subroutine D3_Find_Crack_Mesh_Outline(i_C,num_of_ele,num_of_nodes, c_Crack3D_Meshed_Ele)
!     Find the outer boundary of the crack surface mesh and store it in the global variable Crack3D_Meshed_Outline(Max_Num_Cr,1000,3) for each crack.
!     The number of boundary lines is stored in the global variable Crack3D_Meshed_Outline_num(Max_Num_Cr).
!
!     Firstly written by Fang Shi on 2019-03-28.
!     integer Crack3D_Meshed_Outline(Max_Num_Cr_3D, Max_N_Node_3D, 4)   !3D crack outer boundary after discretization
                                                                  ! data 1 is the first point on the boundary line of the crack front edge
                                                                  ! data 2 is the second point on the leading edge boundary line of the crack
                                                                  ! data 3 corresponds to the discrete fracture element number
                                                                  ! data 4 is used to mark whether the two points of the boundary line are allowed to extend,
                                                                  ! extending in very small steps (2021-08-20)
!     integer Crack3D_Meshed_Outline_num(Max_Num_Cr)                    !Number of 3D crack boundary lines after discretization


! The following variables Phi and Psi were not used
!     real(kind=FT) Crack3D_Meshed_Outline_Vertex(Max_Num_Cr_3D, Max_N_Node_3D, 3)  !Coordinates of the 3D crack outer boundary vertices after discretization
!     integer Crack3D_Meshed_Outline_Vertex_Ele_num(Max_Num_Cr_3D, Max_N_Node_3D)  !Solid element number of the 3D crack outer boundary vertex after discretization

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
integer,intent(in)::i_C,num_of_ele,num_of_nodes
integer,intent(in)::c_Crack3D_Meshed_Ele(num_of_ele,3)
integer :: tem1(3*num_of_ele,2)
integer :: tem2(3*num_of_ele)
integer :: i,j,all_num_Outline,Out_num_Outline
integer,allocatable::All_Outline(:,:)
integer,allocatable::Temp_Outline(:,:)
integer :: location
logical :: Yes_In_1,Yes_In_2
integer temp_vector(3)


!     ---------------------
!     Main program section
!     ---------------------
! The three sides of the element are stored in the temporary variable tem1
tem1(1:num_of_ele,1)   = c_Crack3D_Meshed_Ele(1:num_of_ele,1)
tem1(1:num_of_ele,2)   = c_Crack3D_Meshed_Ele(1:num_of_ele,2)
tem1(num_of_ele+1:2*num_of_ele,1)   = c_Crack3D_Meshed_Ele(1:num_of_ele,2)
tem1(num_of_ele+1:2*num_of_ele,2)   = c_Crack3D_Meshed_Ele(1:num_of_ele,3)
tem1(2*num_of_ele+1:3*num_of_ele,1) = c_Crack3D_Meshed_Ele(1:num_of_ele,3)
tem1(2*num_of_ele+1:3*num_of_ele,2) = c_Crack3D_Meshed_Ele(1:num_of_ele,1)

! Sort
call Matrix_Sort_Int(3*num_of_ele,2,tem1)

! Count the occurrences of each row in the matrix and store them in tem2
call Matrix_Count_Row_Int(3*num_of_ele,2, tem1,tem2,all_num_Outline)

! Extract edges that appear only once
allocate(All_Outline(all_num_Outline,2))
allocate(Temp_Outline(all_num_Outline,2))
j=0
do i=1,3*num_of_ele
    if (tem2(i).eq.1)then
        j=j+1
        All_Outline(j,:) = tem1(i,:)
    end if
end do

! Connected end to end, containing only the outer boundary
call Tool_Sort_by_End_to_End(all_num_Outline,all_num_Outline, All_Outline,Temp_Outline, Out_num_Outline)

Crack3D_Meshed_Outline_num(i_C) = Out_num_Outline

! Find the elements corresponding to the outer boundary of the crack
do i=1,Out_num_Outline
    Crack3D_Meshed_Outline(i_C)%row(i,1:2)= Temp_Outline(i,1:2)
    do j=1,num_of_ele
        temp_vector(1:3) = c_Crack3D_Meshed_Ele(j,1:3)
call Vector_Location_Int(3,temp_vector(1:3), Temp_Outline(i,1),location,Yes_In_1)
call Vector_Location_Int(3,temp_vector(1:3), Temp_Outline(i,2),location,Yes_In_2)
        ! if the two nodes on the current outer boundary belong to a certain element, then that element is
        ! the one we are looking for.
        if (Yes_In_1 .and. Yes_In_2)then
            Crack3D_Meshed_Outline(i_C)%row(i,3)  = j
            exit
        end if
    end do
end do
if (allocated(All_Outline))deallocate(All_Outline)
if (allocated(Temp_Outline))deallocate(Temp_Outline)



return
end subroutine D3_Find_Crack_Mesh_Outline
