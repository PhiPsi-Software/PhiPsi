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

SUBROUTINE Cut_and_Save_Boundary_Cracks_3D(isub,c_DISP)
!================================================================
! First written on 2026-02-08. Updated on 2026-02-10.
! Cut 3D boundary cracks: clip the crack triangular mesh against
! the outer boundary of the 8-node hex FE mesh.
! Only portions inside the model are retained.
!
! Algorithm:
!   1. Extract all hex faces, sort & find boundary faces (shared by 1 elem)
!   2. Split each boundary quad into 2 triangles
!   3. For each crack triangle:
!        - test vertices inside/outside (ray casting)
!        - clip partial triangles at boundary intersections
!   4. Update crack mesh data
!================================================================
!
!New file name: cnox -> cnxc; cnoy -> cnyc; cnoz -> cnzc;
!               cms1 -> cmc1; cms2 -> cmc2; cms3 -> cmc3; 
!               cmap -> ccap. 
!
!Regarding the connectivity of hexahedral faces: hex_face defines
!     the 6 faces of an 8-node hexahedron, each with 4 local node numbers.
!     If your finite element program uses a different node 
!     numbering convention, you only need to modify these 6 lines. 
!     For boundary face extraction (finding duplicate faces), the 
!     orientation of the faces does not affect the result because 
!     node numbers are sorted during comparison.

!Ray-casting method for inside/outside determination: bcrk_pt_in_model
!     casts a ray from the query point in an "irrational number" direction 
!     and counts the number of intersections with all boundary &
!     triangles. An odd number of intersections indicates the point is 
!     inside. The ray direction is chosen as (1.0, 0.1234..., &
!     0.0987...) to minimize the chance of exactly passing through 
!     an edge or vertex of a triangle. For fracture vertices located &
!     exactly on the model's boundary surface, the condition tv > TOL 
!     automatically excludes contacts with t, usually'// &
!     'ensuring the correctness of the inside/outside judgment.

!Clipping logic: For the "1 vertex inside" case, one triangle is 
!     retained. For the "2 vertices inside" case, the clipped trapezoid
!     is split into two triangles. bcrk_seg_isect always starts from 
!     the inside point to find the first intersection with the boundary 
!     (minimum t), ensuring the clipping location is correct.

!Performance: Currently, the inside/outside determination for each 
!     fracture vertex and the intersection calculation for each edge 
!     are performed by iterating through all boundary triangles 
!     (brute-force method). For large models, pre-screening with 
!     bounding boxes or spatial indexing (e.g., octrees, BVH) can be 
!     added for acceleration.

use Global_Float_Type
use Global_Common
use Global_Crack_Common
use Global_Crack_3D
use Global_Crack
use Global_Model
use Global_Filename
use Global_POST
use Global_Elem_Area_Vol
use Global_Ragged_Array_Real_Classs
use Global_Ragged_Array_Int_Classs
use module_Cal_Crack_Point_Aperture_3D

implicit none
integer, intent(in) :: isub
real(kind=FT),intent(in)::c_DISP(Total_FD)

! ---- local scalars ----
integer :: i_C, i_E, i, k, f, ii
integer :: Number_of_3D_crack_nodes, Number_of_3D_crack_Eles
integer :: n1_idx, n2_idx, n3_idx
integer :: n_in
integer :: num_all_faces, num_bfaces, num_btri
integer :: max_nn, max_ne, cnt_n, cnt_e
integer :: idxA, idxB, idxC, ni1, ni2

! ---- small fixed-size arrays ----
real(kind=FT) :: c_X(8), c_Y(8), c_Z(8)
integer :: c_NN(8)
integer :: hex_face(4,6)
real(kind=FT) :: v1(3), v2(3), v3(3)
real(kind=FT) :: pA(3), pB(3), pC(3)
real(kind=FT) :: ip1(3), ip2(3)
logical :: in1, in2, in3, found1, found2

! ---- allocatable work arrays ----
integer, allocatable :: face_srt(:,:)
real(kind=FT), allocatable :: face_xyz(:,:,:)
logical, allocatable :: is_bnd(:)
integer, allocatable :: srt_ix(:)
real(kind=FT), allocatable :: bt(:,:,:)
real(kind=FT), allocatable :: nn_xyz(:,:)
integer, allocatable :: ne_conn(:,:)
type(Ragged_Array_2D),allocatable::Crack3D_Meshed_Node_CUT(:)
type(Ragged_Int_Array_2D),allocatable::Crack3D_Meshed_Ele_CUT(:)
type(Ragged_Array_2D),allocatable::Crack3D_Meshed_Node_Value_CUT(:)
integer Crack3D_Meshed_Node_num_CUT(num_Crack)
integer Crack3D_Meshed_Ele_num_CUT(num_Crack)
character(200) Filename_1,Filename_2,Filename_3
character(5) temp    
integer j,i_Crack_Node
real(kind=FT) ::Crack_Node_Coor(3),c_Aperture
logical on_boundary
logical, allocatable :: Crack_Boundary_Node_Flag(:,:)
integer :: max_nodes_all
real(kind=FT) :: Relative_Disp(3)
logical, allocatable :: is_boundary_node(:)
integer, allocatable :: old_to_new(:)
real(kind=FT) :: ori_n(3)
real(kind=FT), allocatable :: ele_normals(:,:)
real(kind=FT) :: e1(3), e2(3), norm_vec(3), norm_len
real(kind=FT) :: total_area
real(kind=FT) :: cross_prod(3)
real(kind=FT) :: area
real(kind=FT) :: moved_point(3), move_dist
integer :: adjacent_elem_count
real(kind=FT), allocatable :: ele_centers(:,:)
real(kind=FT) :: MAX_MOVE_RATIO
real(kind=FT) :: ref_normal(3), tri_normal(3), dot_ref
logical :: ref_normal_set
integer :: temp_idx
real(kind=FT) :: temp_v(3)
logical, allocatable :: ele_valid(:)
integer n1,n2,n3,ref_idx

!====================================================
!If all cracks are not boundary cracks, then return.
!====================================================
if (all(.not. Boundary_Cracks)) then
    goto 200
end if

print *,'    Cut and save boundary cracks...'
allocate(Crack3D_Meshed_Node_CUT(num_Crack))
allocate(Crack3D_Meshed_Ele_CUT(num_Crack))
allocate(Crack3D_Meshed_Node_Value_CUT(num_Crack))

! ============================================================
! 1. Define local face connectivity for 8-node hexahedron
!
!         5-------8          Face 1 (bottom): 1 4 3 2
!        /|      /|          Face 2 (top)   : 5 6 7 8
!       / |     / |          Face 3 (front) : 1 2 6 5
!      6-------7  |          Face 4 (right) : 2 3 7 6
!      |  1----|--4          Face 5 (back)  : 3 4 8 7
!      | /     | /           Face 6 (left)  : 4 1 5 8
!      |/      |/
!      2-------3
!
!   NOTE: Adjust if your code uses a different node ordering.
! ============================================================
hex_face(:,1) = (/1, 4, 3, 2/)
hex_face(:,2) = (/5, 6, 7, 8/)
hex_face(:,3) = (/1, 2, 6, 5/)
hex_face(:,4) = (/2, 3, 7, 6/)
hex_face(:,5) = (/3, 4, 8, 7/)
hex_face(:,6) = (/4, 1, 5, 8/)

! =============================
! 2. Extract all element faces
! =============================
num_all_faces = Num_Elem * 6
allocate(face_srt(4, num_all_faces))
allocate(face_xyz(3, 4, num_all_faces))
allocate(is_bnd(num_all_faces))
allocate(srt_ix(num_all_faces))

k = 0
do i_E = 1, Num_Elem
    c_NN = G_NN(1:8, i_E)
    c_X  = G_X_NODES(1:8, i_E)
    c_Y  = G_Y_NODES(1:8, i_E)
    c_Z  = G_Z_NODES(1:8, i_E)
    do f = 1, 6
        k = k + 1
        do i = 1, 4
            ii = hex_face(i, f)
            face_srt(i, k) = c_NN(ii)
            face_xyz(1, i, k) = c_X(ii)
            face_xyz(2, i, k) = c_Y(ii)
            face_xyz(3, i, k) = c_Z(ii)
        end do
        ! Sort 4 node IDs for duplicate detection
        call bcrk_sort4(face_srt(:, k))
        srt_ix(k) = k
    end do
end do

! ==============================================================
! 3. Sort faces lexicographically, then find boundary faces
!    (faces that appear only ONCE = not shared by two elements)
! ==============================================================
call bcrk_qsort_idx(face_srt, srt_ix, 1, num_all_faces)

is_bnd = .false.
i = 1
do while (i <= num_all_faces)
    if (i < num_all_faces) then
        if (bcrk_face_eq(face_srt(:, srt_ix(i)), &
                         face_srt(:, srt_ix(i+1)))) then
            ! Internal face (shared by two elements) -> skip both
            i = i + 2
            cycle
        end if
    end if
    ! Boundary face (appears only once)
    is_bnd(srt_ix(i)) = .true.
    i = i + 1
end do

deallocate(srt_ix)

! =============================================
! 4. Build boundary triangle list
!    Each quad face is split into 2 triangles:
!       Triangle A : nodes 1-2-3
!       Triangle B : nodes 1-3-4
! =============================================
num_bfaces = 0
do i = 1, num_all_faces
    if (is_bnd(i)) num_bfaces = num_bfaces + 1
end do
num_btri = num_bfaces * 2
allocate(bt(3, 3, num_btri))

k = 0
do i = 1, num_all_faces
    if (.not. is_bnd(i)) cycle
    ! Triangle A: quad nodes 1-2-3
    k = k + 1
    bt(:, 1, k) = face_xyz(:, 1, i)
    bt(:, 2, k) = face_xyz(:, 2, i)
    bt(:, 3, k) = face_xyz(:, 3, i)
    ! Triangle B: quad nodes 1-3-4
    k = k + 1
    bt(:, 1, k) = face_xyz(:, 1, i)
    bt(:, 2, k) = face_xyz(:, 3, i)
    bt(:, 3, k) = face_xyz(:, 4, i)
end do

deallocate(face_srt, face_xyz, is_bnd)



! ======================================================
! 5. Loop over each boundary crack and perform clipping
! ======================================================
max_nodes_all = 0
do i_C = 1, num_Crack
    if (.not. Boundary_Cracks(i_C)) cycle
    max_nodes_all = max(max_nodes_all, &
                        Crack3D_Meshed_Node_num(i_C) + 2 * Crack3D_Meshed_Ele_num(i_C))
end do
allocate(Crack_Boundary_Node_Flag(max_nodes_all, num_Crack))
Crack_Boundary_Node_Flag = .false.

do i_C = 1, num_Crack
    if (.not. Boundary_Cracks(i_C)) cycle

    Number_of_3D_crack_nodes = Crack3D_Meshed_Node_num(i_C)
    Number_of_3D_crack_Eles  = Crack3D_Meshed_Ele_num(i_C)
    if (Number_of_3D_crack_Eles == 0) cycle

    max_nn = Number_of_3D_crack_nodes + 2 * Number_of_3D_crack_Eles
    max_ne = 2 * Number_of_3D_crack_Eles
    allocate(nn_xyz(max_nn, 3))
    allocate(ne_conn(max_ne, 3))
    allocate(is_boundary_node(max_nn))
    is_boundary_node = .false.

    allocate(old_to_new(Number_of_3D_crack_nodes))
    old_to_new = 0
    
    cnt_n = 0
    cnt_e = 0
    
    !Initialize reference normal
    ref_normal_set = .false.
    ref_normal = (/0.0d0, 0.0d0, 0.0d0/)

    
    do i = 1, Number_of_3D_crack_Eles
        n1_idx = Crack3D_Meshed_Ele(i_C)%row(i, 1)
        n2_idx = Crack3D_Meshed_Ele(i_C)%row(i, 2)
        n3_idx = Crack3D_Meshed_Ele(i_C)%row(i, 3)

        v1(1:3) = Crack3D_Meshed_Node(i_C)%row(n1_idx, 1:3)
        v2(1:3) = Crack3D_Meshed_Node(i_C)%row(n2_idx, 1:3)
        v3(1:3) = Crack3D_Meshed_Node(i_C)%row(n3_idx, 1:3)
        
        !Compute the normal vector of the current triangle.
        e1 = v2 - v1
        e2 = v3 - v1
        tri_normal(1) = e1(2)*e2(3) - e1(3)*e2(2)
        tri_normal(2) = e1(3)*e2(1) - e1(1)*e2(3)
        tri_normal(3) = e1(1)*e2(2) - e1(2)*e2(1)
        norm_len = sqrt(tri_normal(1)**2 + tri_normal(2)**2 + tri_normal(3)**2)
        
        if (norm_len > 1.0d-12) then
            tri_normal = tri_normal / norm_len
            
            if (.not. ref_normal_set) then
                ref_normal = tri_normal
                ref_normal_set = .true.
            else
                dot_ref = tri_normal(1)*ref_normal(1) + &
                          tri_normal(2)*ref_normal(2) + &
                          tri_normal(3)*ref_normal(3)
                
                if (dot_ref < 0.0d0) then
                    temp_idx = n2_idx
                    n2_idx = n3_idx
                    n3_idx = temp_idx
                    
                    temp_v = v2
                    v2 = v3
                    v3 = temp_v
                end if
            end if
        end if
        
        call bcrk_pt_in_model(v1, bt, num_btri, in1)
        call bcrk_pt_in_model(v2, bt, num_btri, in2)
        call bcrk_pt_in_model(v3, bt, num_btri, in3)

        n_in = 0
        if (in1) n_in = n_in + 1
        if (in2) n_in = n_in + 1
        if (in3) n_in = n_in + 1

        if (n_in == 0) cycle

        ! +++++ Case 3: All points inside +++++
        if (n_in == 3) then
            if (old_to_new(n1_idx) == 0) then
                cnt_n = cnt_n + 1
                old_to_new(n1_idx) = cnt_n
                nn_xyz(cnt_n, :) = v1
            end if
            if (old_to_new(n2_idx) == 0) then
                cnt_n = cnt_n + 1
                old_to_new(n2_idx) = cnt_n
                nn_xyz(cnt_n, :) = v2
            end if
            if (old_to_new(n3_idx) == 0) then
                cnt_n = cnt_n + 1
                old_to_new(n3_idx) = cnt_n
                nn_xyz(cnt_n, :) = v3
            end if
            
            cnt_e = cnt_e + 1
            ne_conn(cnt_e, :) = (/old_to_new(n1_idx), &
                                  old_to_new(n2_idx), &
                                  old_to_new(n3_idx)/)
            cycle
        end if

        ! +++++ Case 1: One point inside +++++
        if (n_in == 1) then
            if (in1) then
                idxA = n1_idx; pA = v1
                idxB = n2_idx; pB = v2
                idxC = n3_idx; pC = v3
            else if (in2) then
                idxA = n2_idx; pA = v2
                idxB = n3_idx; pB = v3
                idxC = n1_idx; pC = v1
            else
                idxA = n3_idx; pA = v3
                idxB = n1_idx; pB = v1
                idxC = n2_idx; pC = v2
            end if

            if (old_to_new(idxA) == 0) then
                cnt_n = cnt_n + 1
                old_to_new(idxA) = cnt_n
                nn_xyz(cnt_n, :) = pA
            end if

            call bcrk_seg_isect(pA, pB, bt, num_btri, ip1, found1)
            if (.not. found1) ip1 = 0.5d0 * (pA + pB)

            call bcrk_seg_isect(pA, pC, bt, num_btri, ip2, found2)
            if (.not. found2) ip2 = 0.5d0 * (pA + pC)

            cnt_n = cnt_n + 1;  ni1 = cnt_n;  nn_xyz(cnt_n, :) = ip1
            cnt_n = cnt_n + 1;  ni2 = cnt_n;  nn_xyz(cnt_n, :) = ip2
            is_boundary_node(ni1) = .true.
            is_boundary_node(ni2) = .true.
            
            !Check normal of new triangle.
            e1 = nn_xyz(ni1, :) - nn_xyz(old_to_new(idxA), :)
            e2 = nn_xyz(ni2, :) - nn_xyz(old_to_new(idxA), :)
            tri_normal(1) = e1(2)*e2(3) - e1(3)*e2(2)
            tri_normal(2) = e1(3)*e2(1) - e1(1)*e2(3)
            tri_normal(3) = e1(1)*e2(2) - e1(2)*e2(1)
            norm_len = sqrt(tri_normal(1)**2 + tri_normal(2)**2 + tri_normal(3)**2)
            
            if (norm_len > 1.0d-12) then
                tri_normal = tri_normal / norm_len
                dot_ref = tri_normal(1)*ref_normal(1) + &
                          tri_normal(2)*ref_normal(2) + &
                          tri_normal(3)*ref_normal(3)
                
                cnt_e = cnt_e + 1
                if (dot_ref > 0.0d0) then
                    ! Normal in the correct direction - keep node order.
                    ne_conn(cnt_e, :) = (/old_to_new(idxA), ni1, ni2/)
                else
                    ! Normal is flipped swap the last two nodes to correct orientation.
                    ne_conn(cnt_e, :) = (/old_to_new(idxA), ni2, ni1/)
                end if
            else
                ! Degenerate triangle (zero or near-zero area); still add it as-is.
                cnt_e = cnt_e + 1
                ne_conn(cnt_e, :) = (/old_to_new(idxA), ni1, ni2/)
            end if
            cycle
        end if

        !       Case 2: Two on the inside      
        if (n_in == 2) then
            if (.not. in1) then
                idxA = n1_idx; pA = v1
                idxB = n2_idx; pB = v2
                idxC = n3_idx; pC = v3
            else if (.not. in2) then
                idxA = n2_idx; pA = v2
                idxB = n3_idx; pB = v3
                idxC = n1_idx; pC = v1
            else
                idxA = n3_idx; pA = v3
                idxB = n1_idx; pB = v1
                idxC = n2_idx; pC = v2
            end if

            if (old_to_new(idxB) == 0) then
                cnt_n = cnt_n + 1
                old_to_new(idxB) = cnt_n
                nn_xyz(cnt_n, :) = pB
            end if
            if (old_to_new(idxC) == 0) then
                cnt_n = cnt_n + 1
                old_to_new(idxC) = cnt_n
                nn_xyz(cnt_n, :) = pC
            end if

            call bcrk_seg_isect(pB, pA, bt, num_btri, ip1, found1)
            if (.not. found1) ip1 = 0.5d0 * (pB + pA)

            call bcrk_seg_isect(pC, pA, bt, num_btri, ip2, found2)
            if (.not. found2) ip2 = 0.5d0 * (pC + pA)

            cnt_n = cnt_n + 1;  ni1 = cnt_n;  nn_xyz(cnt_n, :) = ip1
            cnt_n = cnt_n + 1;  ni2 = cnt_n;  nn_xyz(cnt_n, :) = ip2
            is_boundary_node(ni1) = .true.
            is_boundary_node(ni2) = .true.
            
            ! +++++ Triangle 1: B -> ni1 -> ni2 +++++
            e1 = nn_xyz(ni1, :) - nn_xyz(old_to_new(idxB), :)
            e2 = nn_xyz(ni2, :) - nn_xyz(old_to_new(idxB), :)
            tri_normal(1) = e1(2)*e2(3) - e1(3)*e2(2)
            tri_normal(2) = e1(3)*e2(1) - e1(1)*e2(3)
            tri_normal(3) = e1(1)*e2(2) - e1(2)*e2(1)
            norm_len = sqrt(tri_normal(1)**2 + tri_normal(2)**2 + tri_normal(3)**2)
            
            cnt_e = cnt_e + 1
            if (norm_len > 1.0d-12) then
                tri_normal = tri_normal / norm_len
                dot_ref = tri_normal(1)*ref_normal(1) + &
                          tri_normal(2)*ref_normal(2) + &
                          tri_normal(3)*ref_normal(3)
                
                if (dot_ref > 0.0d0) then
                    ne_conn(cnt_e, :) = (/old_to_new(idxB), ni1, ni2/)
                else
                    ne_conn(cnt_e, :) = (/old_to_new(idxB), ni2, ni1/)
                end if
            else
                ne_conn(cnt_e, :) = (/old_to_new(idxB), ni1, ni2/)
            end if
            
            ! +++++ Triangle 2: B -> ni2 -> C +++++
            e1 = nn_xyz(ni2, :) - nn_xyz(old_to_new(idxB), :)
            e2 = nn_xyz(old_to_new(idxC), :) - nn_xyz(old_to_new(idxB), :)
            tri_normal(1) = e1(2)*e2(3) - e1(3)*e2(2)
            tri_normal(2) = e1(3)*e2(1) - e1(1)*e2(3)
            tri_normal(3) = e1(1)*e2(2) - e1(2)*e2(1)
            norm_len = sqrt(tri_normal(1)**2 + tri_normal(2)**2 + tri_normal(3)**2)
            
            cnt_e = cnt_e + 1
            if (norm_len > 1.0d-12) then
                tri_normal = tri_normal / norm_len
                dot_ref = tri_normal(1)*ref_normal(1) + &
                          tri_normal(2)*ref_normal(2) + &
                          tri_normal(3)*ref_normal(3)
                
                if (dot_ref > 0.0d0) then
                    ne_conn(cnt_e, :) = (/old_to_new(idxB), ni2, old_to_new(idxC)/)
                else
                    ne_conn(cnt_e, :) = (/old_to_new(idxB), old_to_new(idxC), ni2/)
                end if
            else
                ne_conn(cnt_e, :) = (/old_to_new(idxB), ni2, old_to_new(idxC)/)
            end if
            
            cycle
        end if
        
    end do
    
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! +++++ Additional step: Remove duplicate intersection nodes.   +++++
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if (cnt_n > 0 .and. cnt_e > 0) then
        call merge_duplicate_nodes(nn_xyz, ne_conn, is_boundary_node, &
                                   cnt_n, cnt_e, 1.0d-8)
    end if
    
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! +++++ Additional step: Re-check and unify normals of all triangles +++++
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if (cnt_e > 0) then
        ! Step 1: Compute normals for all triangles.
        if (allocated(ele_normals)) deallocate(ele_normals)
        if (allocated(ele_valid)) deallocate(ele_valid)
        allocate(ele_normals(cnt_e, 3))
        allocate(ele_valid(cnt_e))
        ele_valid = .false.
        
        do i = 1, cnt_e
            n1 = ne_conn(i, 1)
            n2 = ne_conn(i, 2)
            n3 = ne_conn(i, 3)
            
            v1 = nn_xyz(n1, :)
            v2 = nn_xyz(n2, :)
            v3 = nn_xyz(n3, :)
            
            e1 = v2 - v1
            e2 = v3 - v1
            
            tri_normal(1) = e1(2)*e2(3) - e1(3)*e2(2)
            tri_normal(2) = e1(3)*e2(1) - e1(1)*e2(3)
            tri_normal(3) = e1(1)*e2(2) - e1(2)*e2(1)
            norm_len = sqrt(tri_normal(1)**2 + tri_normal(2)**2 + tri_normal(3)**2)

            if (norm_len > 1.0d-12) then
                ele_normals(i, :) = tri_normal / norm_len
                ele_valid(i) = .true.
            else
                ele_normals(i, :) = 0.0d0
                ele_valid(i) = .false.
            end if
        end do
        ! Step 2: Find the first valid triangle as reference.
        ref_idx = 0
        do i = 1, cnt_e
            if (ele_valid(i)) then
                ref_idx = i
                ref_normal = ele_normals(i, :)
                exit
            end if
        end do
        ! Step 3: Check and correct all other triangles.
        if (ref_idx > 0) then
            do i = 1, cnt_e
                if (.not. ele_valid(i)) cycle
                if (i == ref_idx) cycle
                dot_ref = ele_normals(i, 1)*ref_normal(1) + &
                          ele_normals(i, 2)*ref_normal(2) + &
                          ele_normals(i, 3)*ref_normal(3)
                if (dot_ref < 0.0d0) then
                    ! Normal is opposite; swap node 2 and node 3.
                    temp_idx = ne_conn(i, 2)
                    ne_conn(i, 2) = ne_conn(i, 3)
                    ne_conn(i, 3) = temp_idx
                end if
            end do
        end if
        deallocate(ele_normals, ele_valid)
    end if
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    ! Write data.
    Crack3D_Meshed_Node_num_CUT(i_C) = cnt_n
    Crack3D_Meshed_Ele_num_CUT(i_C)  = cnt_e

    if (allocated(Crack3D_Meshed_Node_CUT(i_C)%row)) &
        deallocate(Crack3D_Meshed_Node_CUT(i_C)%row)
    if (allocated(Crack3D_Meshed_Ele_CUT(i_C)%row)) &
        deallocate(Crack3D_Meshed_Ele_CUT(i_C)%row)
    if (allocated(Crack3D_Meshed_Node_Value_CUT(i_C)%row)) &
        deallocate(Crack3D_Meshed_Node_Value_CUT(i_C)%row)

    allocate(Crack3D_Meshed_Node_CUT(i_C)%row(cnt_n, 3))
    allocate(Crack3D_Meshed_Node_Value_CUT(i_C)%row(cnt_n, 3))
    allocate(Crack3D_Meshed_Ele_CUT(i_C)%row(cnt_e, 3))

    Crack3D_Meshed_Node_CUT(i_C)%row(1:cnt_n, 1:3) = nn_xyz(1:cnt_n, 1:3)
    
    do i = 1, cnt_e
        Crack3D_Meshed_Ele_CUT(i_C)%row(i, 1:3) = ne_conn(i, 1:3)
    end do
    
    ! Save boundary node flags.
    do i = 1, cnt_n
        if (is_boundary_node(i)) then
            Crack3D_Meshed_Node_Value_CUT(i_C)%row(i, 2) = 1.0d0
        else
            Crack3D_Meshed_Node_Value_CUT(i_C)%row(i, 2) = 0.0d0
        end if
    end do
    
    ! Save boundary flags to global array.
    Crack_Boundary_Node_Flag(1:cnt_n, i_C) = is_boundary_node(1:cnt_n)
    
    deallocate(nn_xyz, ne_conn, is_boundary_node, old_to_new)
end do




! ======================================================
! 6. Calculating apertures of new mesh nodes of cracks.
! ======================================================
do i_C = 1, num_Crack
    if (.not. Boundary_Cracks(i_C)) cycle
    
    print *,'    Computing apertures for boundary cutting cracks...'
    
    ! Precompute normal vectors for all elements.
    if(allocated(ele_normals)) deallocate(ele_normals)
    if(allocated(ele_centers)) deallocate(ele_centers)
    allocate(ele_normals(3, Crack3D_Meshed_Ele_num_CUT(i_C)))
    allocate(ele_centers(3, Crack3D_Meshed_Ele_num_CUT(i_C)))
    do i = 1, Crack3D_Meshed_Ele_num_CUT(i_C)
        n1_idx = Crack3D_Meshed_Ele_CUT(i_C)%row(i, 1)
        n2_idx = Crack3D_Meshed_Ele_CUT(i_C)%row(i, 2)
        n3_idx = Crack3D_Meshed_Ele_CUT(i_C)%row(i, 3)
        v1 = Crack3D_Meshed_Node_CUT(i_C)%row(n1_idx, :)
        v2 = Crack3D_Meshed_Node_CUT(i_C)%row(n2_idx, :)
        v3 = Crack3D_Meshed_Node_CUT(i_C)%row(n3_idx, :)
        ! Compute centroid of the triangular fracture face element.
        ele_centers(:, i) = (v1 + v2 + v3) / 3.0d0
        ! Compute two edge vectors.
        e1 = v2 - v1
        e2 = v3 - v1
        ! Compute normal vector via cross product.
        norm_vec(1) = e1(2)*e2(3) - e1(3)*e2(2)
        norm_vec(2) = e1(3)*e2(1) - e1(1)*e2(3)
        norm_vec(3) = e1(1)*e2(2) - e1(2)*e2(1)
        !! Normalize the vector.
        norm_len = sqrt(norm_vec(1)**2 + norm_vec(2)**2 + norm_vec(3)**2)
        if (norm_len > 1.0d-12) then
            ele_normals(:, i) = norm_vec / norm_len
        else
            ele_normals(:, i) = (/0.0d0, 0.0d0, 1.0d0/)
        end if
    end do
    
    !Loop over new crack points.
    do i_Crack_Node = 1, Crack3D_Meshed_Node_num_CUT(i_C)
        Crack_Node_Coor = Crack3D_Meshed_Node_CUT(i_C)%row(i_Crack_Node, 1:3)
        ! Check if it is a boundary node.
        on_boundary = Crack_Boundary_Node_Flag(i_Crack_Node, i_C)
        ! Compute nodal normal as area-weighted average of adjacent element normals.
        ori_n = 0.0d0
        total_area = 0.0d0
        adjacent_elem_count = 0
        moved_point = 0.0d0
        
        do i = 1, Crack3D_Meshed_Ele_num_CUT(i_C)
            ! Check if the current node belongs to this element.
            if (any(Crack3D_Meshed_Ele_CUT(i_C)%row(i, :) == i_Crack_Node)) then
                adjacent_elem_count = adjacent_elem_count + 1
                n1_idx = Crack3D_Meshed_Ele_CUT(i_C)%row(i, 1)
                n2_idx = Crack3D_Meshed_Ele_CUT(i_C)%row(i, 2)
                n3_idx = Crack3D_Meshed_Ele_CUT(i_C)%row(i, 3)
                v1 = Crack3D_Meshed_Node_CUT(i_C)%row(n1_idx, :)
                v2 = Crack3D_Meshed_Node_CUT(i_C)%row(n2_idx, :)
                v3 = Crack3D_Meshed_Node_CUT(i_C)%row(n3_idx, :)
                e1 = v2 - v1
                e2 = v3 - v1
                cross_prod(1) = e1(2)*e2(3) - e1(3)*e2(2)
                cross_prod(2) = e1(3)*e2(1) - e1(1)*e2(3)
                cross_prod(3) = e1(1)*e2(2) - e1(2)*e2(1)
                area = 0.5d0 * sqrt(cross_prod(1)**2 + cross_prod(2)**2 + cross_prod(3)**2)
                ori_n = ori_n + area * ele_normals(:, i)
                total_area = total_area + area
                ! For boundary nodes, accumulate element centers (used later for point relocation).
                if (on_boundary) then
                    moved_point = moved_point + ele_centers(:, i)
                end if
            end if
        end do

        ! Normalize the computed nodal normal vector.
        if (total_area > 1.0d-12) then
            ori_n = ori_n / sqrt(ori_n(1)**2 + ori_n(2)**2 + ori_n(3)**2)
        else
            ori_n = (/0.0d0, 0.0d0, 1.0d0/)  
        end if
        
        !......................................................
        ! On boundary: move point and then calculate aperture.
        !......................................................
        if (on_boundary) then
            ! Move toward the center of the triangles containing the current point
            if (adjacent_elem_count > 0) then
                    ! Calculate the average position of adjacent element centers
                    moved_point = moved_point / real(adjacent_elem_count, kind=FT)
                    ! Limit the movement distance (avoid moving too far)
                    move_dist = sqrt(sum((moved_point - Crack_Node_Coor)**2))
                    MAX_MOVE_RATIO = 0.1D0
                    if (move_dist > MAX_MOVE_RATIO * Ave_Elem_L_Enrich) then
                        ! Limit the movement distance
                        moved_point = Crack_Node_Coor + &
                            (moved_point - Crack_Node_Coor)*(MAX_MOVE_RATIO*Ave_Elem_L_Enrich/move_dist)
                    end if
                    ! Calculate aperture at the moved point
                    call Cal_Crack_Point_Aperture_3D(c_DISP, i_C, moved_point,Relative_Disp(1:3), 0, 0, 0, 0)
                    ! Calculate normal aperture
                    c_Aperture = Relative_Disp(1)*ori_n(1) + Relative_Disp(2)*ori_n(2) + Relative_Disp(3)*ori_n(3)
                    ! Save aperture value
                    Crack3D_Meshed_Node_Value_CUT(i_C)%row(i_Crack_Node, 1) = c_Aperture
                else
                    ! Isolated boundary point: set aperture to 0
                    Crack3D_Meshed_Node_Value_CUT(i_C)%row(i_Crack_Node, 1) = 0.0d0
                    Crack3D_Meshed_Node_Value_CUT(i_C)%row(i_Crack_Node, 3) = 0.0d0
                end if
        ! Not on boundary.
        elseif (.not. on_boundary) then
            call Cal_Crack_Point_Aperture_3D(c_DISP, i_C, Crack_Node_Coor, &
                                             Relative_Disp(1:3), 0, 0, 0, 0)
            !Get aperture.
            c_Aperture = Relative_Disp(1)*ori_n(1)+ Relative_Disp(2)*ori_n(2)+ Relative_Disp(3)*ori_n(3)                             
            ! Save aperture to Crack3D_Meshed_Node_Value_CUT.
            Crack3D_Meshed_Node_Value_CUT(i_C)%row(i_Crack_Node, 1) = c_Aperture
        endif
        
        
    end do
end do

deallocate(bt)

! =================
! 7. Saving files.
! =================
200 continue

print *,'    Saving cnxc, cnyc, cnzc files for boundary cracks cutting...'
write(temp,'(I5)') isub
Filename_1=trim(Full_Pathname)//'.cnxc'//'_'//ADJUSTL(temp)      
Filename_2=trim(Full_Pathname)//'.cnyc'//'_'//ADJUSTL(temp)   
Filename_3=trim(Full_Pathname)//'.cnzc'//'_'//ADJUSTL(temp)        
open(101,file=Filename_1,status='unknown')    
open(102,file=Filename_2,status='unknown') 
open(103,file=Filename_3,status='unknown') 
do i_C=1,num_Crack
    if (.not. Boundary_Cracks(i_C)) then
        write(101, '(50000E20.12)') (Crack3D_Meshed_Node(i_C)%row(j,1),j=1,Crack3D_Meshed_Node_num(i_C))       
        write(102, '(50000E20.12)') (Crack3D_Meshed_Node(i_C)%row(j,2),j=1,Crack3D_Meshed_Node_num(i_C))   
        write(103, '(50000E20.12)') (Crack3D_Meshed_Node(i_C)%row(j,3),j=1,Crack3D_Meshed_Node_num(i_C))        
    else
        write(101, '(50000E20.12)') (Crack3D_Meshed_Node_CUT(i_C)%row(j,1),j=1,Crack3D_Meshed_Node_num_CUT(i_C))       
        write(102, '(50000E20.12)') (Crack3D_Meshed_Node_CUT(i_C)%row(j,2),j=1,Crack3D_Meshed_Node_num_CUT(i_C))   
        write(103, '(50000E20.12)') (Crack3D_Meshed_Node_CUT(i_C)%row(j,3),j=1,Crack3D_Meshed_Node_num_CUT(i_C))     
    endif
end do
close(101)          
close(102)    
close(103)  

print *,'    Saving cmc1, cmc2, cmc3 files for boundary cracks cutting...'
write(temp,'(I5)') isub
Filename_1   =  trim(Full_Pathname)//'.cmc1_'//ADJUSTL(temp)  
Filename_2   =  trim(Full_Pathname)//'.cmc2_'//ADJUSTL(temp)      
Filename_3   =  trim(Full_Pathname)//'.cmc3_'//ADJUSTL(temp)         
open(101,file=Filename_1,status='unknown')    
open(102,file=Filename_2,status='unknown') 
open(103,file=Filename_3,status='unknown')     
do i_C=1,num_Crack
    if (.not. Boundary_Cracks(i_C)) then
        write(101, '(50000I10)') (Crack3D_Meshed_Ele(i_C)%row(j,1),j=1,Crack3D_Meshed_Ele_num(i_C))
        write(102, '(50000I10)') (Crack3D_Meshed_Ele(i_C)%row(j,2),j=1,Crack3D_Meshed_Ele_num(i_C))
        write(103, '(50000I10)') (Crack3D_Meshed_Ele(i_C)%row(j,3),j=1,Crack3D_Meshed_Ele_num(i_C))
    else
        write(101, '(50000I10)') (Crack3D_Meshed_Ele_CUT(i_C)%row(j,1),j=1,Crack3D_Meshed_Ele_num_CUT(i_C))
        write(102, '(50000I10)') (Crack3D_Meshed_Ele_CUT(i_C)%row(j,2),j=1,Crack3D_Meshed_Ele_num_CUT(i_C))
        write(103, '(50000I10)') (Crack3D_Meshed_Ele_CUT(i_C)%row(j,3),j=1,Crack3D_Meshed_Ele_num_CUT(i_C))
    endif
end do
close(101)
close(102)   
close(103)  

print *,'    Saving ccap file for boundary cracks cutting...'
write(temp,'(I5)') isub
Filename_1=trim(Full_Pathname)//'.ccap'//'_'//ADJUSTL(temp)    
open(101,file=Filename_1,status='unknown') 
do i_C=1,num_Crack
    if (.not. Boundary_Cracks(i_C)) then
        if (Key_Del_Neg_Aperture == 1) then
            ! Adjust negative values to 0 before writing
            write(101,'(50000E20.12)')(max(0.0d0, Crack3D_Meshed_Node_Value(i_C)%row(j,1)), &
                                       j=1,Crack3D_Meshed_Node_num(i_C))
        else
            ! Write original values directly
            write(101,'(50000E20.12)')(Crack3D_Meshed_Node_Value(i_C)%row(j,1), &
                                       j=1,Crack3D_Meshed_Node_num(i_C))
        endif
    else
        if (Key_Del_Neg_Aperture == 1) then
            ! Adjust negative values to 0 before writing
            write(101,'(50000E20.12)')(max(0.0d0, Crack3D_Meshed_Node_Value_CUT(i_C)%row(j,1)), &
                                       j=1,Crack3D_Meshed_Node_num_CUT(i_C))
        else
            ! Write original values directly
            write(101,'(50000E20.12)')(Crack3D_Meshed_Node_Value_CUT(i_C)%row(j,1), &
                                       j=1,Crack3D_Meshed_Node_num_CUT(i_C))
        endif
    endif
end do
close(101)

!Deallocate variables.
if (all(.not. Boundary_Cracks)) then
    return
else
    do i_C=1,num_Crack
        if (allocated(Crack3D_Meshed_Ele_CUT(i_C)%row))  deallocate(Crack3D_Meshed_Ele_CUT(i_C)%row)
        if (allocated(Crack3D_Meshed_Node_CUT(i_C)%row)) deallocate(Crack3D_Meshed_Node_CUT(i_C)%row)
        if (allocated(Crack3D_Meshed_Node_Value_CUT(i_C)%row)) deallocate(Crack3D_Meshed_Node_Value_CUT(i_C)%row)
    enddo
endif

    
RETURN

! ===========================================
!                   INTERNAL HELPER ROUTINES
! ===========================================
CONTAINS

!-----------------------------------------------------
! Sort 4 integers in ascending order (selection sort)
!-----------------------------------------------------
subroutine bcrk_sort4(a)
    implicit none
    integer, intent(inout) :: a(4)
    integer :: p, q, tmp
    do p = 1, 3
        do q = p + 1, 4
            if (a(q) < a(p)) then
                tmp = a(p); a(p) = a(q); a(q) = tmp
            end if
        end do
    end do
end subroutine bcrk_sort4

!-------------------------------------------
! Check equality of two sorted 4-node faces
!-------------------------------------------
logical function bcrk_face_eq(f1, f2)
    implicit none
    integer, intent(in) :: f1(4), f2(4)
    bcrk_face_eq = (f1(1)==f2(1) .and. f1(2)==f2(2) .and. &
                    f1(3)==f2(3) .and. f1(4)==f2(4))
end function bcrk_face_eq

!-----------------------------------------------------
! Lexicographic comparison of two sorted 4-node faces
! Returns: -1 (f1 < f2),  0 (equal),  +1 (f1 > f2)
!-----------------------------------------------------
integer function bcrk_face_cmp(f1, f2)
    implicit none
    integer, intent(in) :: f1(4), f2(4)
    integer :: m
    bcrk_face_cmp = 0
    do m = 1, 4
        if (f1(m) < f2(m)) then
            bcrk_face_cmp = -1; return
        else if (f1(m) > f2(m)) then
            bcrk_face_cmp = 1; return
        end if
    end do
end function bcrk_face_cmp

!-----------------------------------------------------------------
! Quicksort an index array by the face node IDs (Hoare partition)
!-----------------------------------------------------------------
recursive subroutine bcrk_qsort_idx(fns, ix, lo, hi)
    implicit none
    integer, intent(in)    :: fns(4, *)
    integer, intent(inout) :: ix(*)
    integer, intent(in)    :: lo, hi
    integer :: i_qs, j_qs, tmp_qs, pivot_f(4)

    if (lo >= hi) return

    pivot_f = fns(:, ix((lo + hi) / 2))
    i_qs = lo;  j_qs = hi

    do while (i_qs <= j_qs)
        do while (bcrk_face_cmp(fns(:, ix(i_qs)), pivot_f) < 0)
            i_qs = i_qs + 1
        end do
        do while (bcrk_face_cmp(fns(:, ix(j_qs)), pivot_f) > 0)
            j_qs = j_qs - 1
        end do
        if (i_qs <= j_qs) then
            tmp_qs = ix(i_qs)
            ix(i_qs) = ix(j_qs)
            ix(j_qs) = tmp_qs
            i_qs = i_qs + 1
            j_qs = j_qs - 1
        end if
    end do

    if (lo  < j_qs) call bcrk_qsort_idx(fns, ix, lo, j_qs)
    if (i_qs < hi)  call bcrk_qsort_idx(fns, ix, i_qs, hi)
end subroutine bcrk_qsort_idx

!----------------------------------------------------------------
! Point-in-model test via ray casting
!
! Cast a ray from pt in a fixed "irrational" direction and count
! intersections with the boundary triangle soup.
! Odd count => inside the closed surface.
!----------------------------------------------------------------
subroutine bcrk_pt_in_model(pt, btri, ntri, inside)
    implicit none
    real(kind=FT), intent(in) :: pt(3)
    integer, intent(in)       :: ntri
    real(kind=FT), intent(in) :: btri(3, 3, ntri)
    logical, intent(out)      :: inside

    real(kind=FT) :: rd(3), tv, uv, vv
    integer :: m_rc, cnt_rc
    logical :: hit_rc
    real(kind=FT), parameter :: TOL = 1.0d-10

    ! Direction chosen to avoid alignment with mesh edges
    rd(1) = 1.0d0
    rd(2) = 0.123456789d0
    rd(3) = 0.0987654321d0

    cnt_rc = 0
    do m_rc = 1, ntri
        call bcrk_ray_tri(pt, rd, &
             btri(:,1,m_rc), btri(:,2,m_rc), btri(:,3,m_rc), &
             tv, uv, vv, hit_rc)
        if (hit_rc .and. tv > TOL) cnt_rc = cnt_rc + 1
    end do

    inside = (mod(cnt_rc, 2) == 1)
end subroutine bcrk_pt_in_model

!------------------------------------------------------------------
! Ray-triangle intersection  (Moller-Trumbore algorithm)
!
! Tests:  ray_o + t * ray_d  against triangle (t0, t1, t2)
! Returns t (ray parameter), u & v (barycentric coords), hit flag.
! Detects both front-face and back-face intersections.
!------------------------------------------------------------------
subroutine bcrk_ray_tri(ray_o, ray_d, t0, t1, t2, t, u, v, hit)
    implicit none
    real(kind=FT), intent(in)  :: ray_o(3), ray_d(3)
    real(kind=FT), intent(in)  :: t0(3), t1(3), t2(3)
    real(kind=FT), intent(out) :: t, u, v
    logical, intent(out)       :: hit

    real(kind=FT) :: e1(3), e2(3), pv(3), sv(3), qv(3)
    real(kind=FT) :: det, inv_det
    real(kind=FT), parameter :: EPS = 1.0d-12

    hit = .false.
    t = 0.0d0;  u = 0.0d0;  v = 0.0d0

    e1 = t1 - t0
    e2 = t2 - t0

    ! pv = cross(ray_d, e2)
    pv(1) = ray_d(2)*e2(3) - ray_d(3)*e2(2)
    pv(2) = ray_d(3)*e2(1) - ray_d(1)*e2(3)
    pv(3) = ray_d(1)*e2(2) - ray_d(2)*e2(1)

    det = e1(1)*pv(1) + e1(2)*pv(2) + e1(3)*pv(3)
    if (abs(det) < EPS) return

    inv_det = 1.0d0 / det
    sv = ray_o - t0

    u = (sv(1)*pv(1) + sv(2)*pv(2) + sv(3)*pv(3)) * inv_det
    if (u < -EPS .or. u > 1.0d0 + EPS) return

    ! qv = cross(sv, e1)
    qv(1) = sv(2)*e1(3) - sv(3)*e1(2)
    qv(2) = sv(3)*e1(1) - sv(1)*e1(3)
    qv(3) = sv(1)*e1(2) - sv(2)*e1(1)

    v = (ray_d(1)*qv(1) + ray_d(2)*qv(2) + ray_d(3)*qv(3)) * inv_det
    if (v < -EPS .or. u + v > 1.0d0 + EPS) return

    t = (e2(1)*qv(1) + e2(2)*qv(2) + e2(3)*qv(3)) * inv_det
    hit = .true.
end subroutine bcrk_ray_tri

!---------------------------------------------------------------
! Line-segment / boundary intersection
!
! For segment P1 -> P2, find the FIRST intersection (smallest t
! in the open interval (0,1)) with all boundary triangles.
! P1 is expected to be inside the model.
!---------------------------------------------------------------
subroutine bcrk_seg_isect(P1, P2, btri, ntri, Pint, found)
    implicit none
    real(kind=FT), intent(in)  :: P1(3), P2(3)
    integer, intent(in)        :: ntri
    real(kind=FT), intent(in)  :: btri(3, 3, ntri)
    real(kind=FT), intent(out) :: Pint(3)
    logical, intent(out)       :: found

    real(kind=FT) :: dir(3), tv, uv, vv, tmin
    integer :: m_si
    logical :: hit_si
    real(kind=FT), parameter :: TOL = 1.0d-10

    dir   = P2 - P1
    found = .false.
    tmin  = 2.0d0

    do m_si = 1, ntri
        call bcrk_ray_tri(P1, dir, &
             btri(:,1,m_si), btri(:,2,m_si), btri(:,3,m_si), &
             tv, uv, vv, hit_si)
        if (hit_si .and. tv > TOL .and. tv < 1.0d0 - TOL) then
            if (tv < tmin) then
                tmin  = tv
                found = .true.
            end if
        end if
    end do

    if (found) then
        Pint = P1 + tmin * dir
    end if
end subroutine bcrk_seg_isect

!-----------------------------------------------------------------
! Merge duplicate nodes (mainly for boundary intersection points)
! tolerance: distance tolerance, nodes closer than this are 
!            considered duplicates.
!-----------------------------------------------------------------
subroutine merge_duplicate_nodes(nodes, elements, bdry_flags, &
                                 num_nodes, num_elems, tolerance)
    implicit none
    real(kind=FT), intent(inout) :: nodes(:,:)
    integer, intent(inout) :: elements(:,:)
    logical, intent(inout) :: bdry_flags(:)
    integer, intent(inout) :: num_nodes, num_elems
    real(kind=FT), intent(in) :: tolerance
    
    integer, allocatable :: node_map(:)
    logical, allocatable :: node_kept(:)
    integer :: i, j, k, new_count
    real(kind=FT) :: dist
    
    allocate(node_map(num_nodes))
    allocate(node_kept(num_nodes))
    
    node_kept = .true.
    node_map = 0
    
    ! Step 1: Mark duplicate nodes (only check boundary nodes for efficiency)
    do i = 1, num_nodes
        if (.not. node_kept(i)) cycle
        node_map(i) = i
        
        ! Only perform duplicate check for boundary nodes
        if (bdry_flags(i)) then
            do j = i + 1, num_nodes
                if (.not. node_kept(j)) cycle
                if (.not. bdry_flags(j)) cycle
                
                ! Calculate distance
                dist = sqrt((nodes(i,1) - nodes(j,1))**2 + &
                           (nodes(i,2) - nodes(j,2))**2 + &
                           (nodes(i,3) - nodes(j,3))**2)
                
                if (dist < tolerance) then
                    ! j is duplicate node, map to i
                    node_kept(j) = .false.
                    node_map(j) = i
                end if
            end do
        end if
    end do
    
    ! Step 2: Renumber kept nodes
    new_count = 0
    do i = 1, num_nodes
        if (node_kept(i)) then
            new_count = new_count + 1
            if (new_count /= i) then
                nodes(new_count, :) = nodes(i, :)
                bdry_flags(new_count) = bdry_flags(i)
            end if
            node_map(i) = new_count
        else
            ! Duplicate node: map to the new number of its master node
            node_map(i) = node_map(node_map(i))
        end if
    end do
    
    ! Step 3: Update element connectivity
    do i = 1, num_elems
        do j = 1, 3
            elements(i, j) = node_map(elements(i, j))
        end do
    end do
    
    ! Step 4: Remove degenerate elements (with duplicate nodes)
    k = 0
    do i = 1, num_elems
        if (elements(i,1) /= elements(i,2) .and. &
            elements(i,2) /= elements(i,3) .and. &
            elements(i,3) /= elements(i,1)) then
            k = k + 1
            if (k /= i) then
                elements(k, :) = elements(i, :)
            end if
        end if
    end do
    num_nodes = new_count
    num_elems = k
    deallocate(node_map, node_kept)
    
end subroutine merge_duplicate_nodes

END SUBROUTINE Cut_and_Save_Boundary_Cracks_3D
