!-----------------------------------------------------------
! Brief: Test intersection between a 3D hydraulic and a natural fracture.
!
! Parameters:
!   Input:  i_HF          - HF crack index
!           i_NF          - NF crack index
!   Output: Yes_Connect   - true if any triangle pair intersects
!           num_Inters    - number of intersection points found
!           Inter_Points  - intersection point coordinates (up to 10)
!
! Notes:   Decomposes NF into triangles and tests each against HF mesh
!          triangles via Tool_Intersections_of_Two_Triangles_3D.
!-----------------------------------------------------------

SUBROUTINE D3_Check_hf_nf_connect(i_HF,i_NF,Yes_Connect,num_Inters,Inter_Points)
! Used to detect the intersection status of hydraulic fractures and natural fractures.
!2023-01-09.

!*****************************
! Read public variable module
!*****************************
use Global_Float_Type
use Global_Crack_Common
use Global_Crack_3D
use Global_Model
use Global_Common

!**********************
! Variable Declaration
!**********************
implicit none
integer,intent(in)::i_HF,i_NF
logical,intent(out)::Yes_Connect
integer,intent(out)::num_Inters
real(kind=FT),intent(out)::Inter_Points(10,3)

integer i_Cr_Node1,i_Cr_Node2,i_Cr_Node3
real(kind=FT) c_Tri_1(3,3),c_Tri_2(3,3)
real(kind=FT) c_Inter_Points(10,3)
logical c_Logical_Inter
integer c_num_Inters
integer i_Crack_Ele
logical c_Logical_Parallel
integer i_Tri                                   

integer num_Vertex


Yes_Connect = .False.
num_Inters  = 0
Inter_Points= ZR

!#################################################
! Determine the intersection points between HF 
! discrete fracture surface triangles and natural 
! fractures.
!#################################################
do i_Crack_Ele =1,Crack3D_Meshed_Ele_num(i_HF)
    i_Cr_Node1 = Crack3D_Meshed_Ele(i_HF)%row(i_Crack_Ele,1)
    i_Cr_Node2 = Crack3D_Meshed_Ele(i_HF)%row(i_Crack_Ele,2)
    i_Cr_Node3 = Crack3D_Meshed_Ele(i_HF)%row(i_Crack_Ele,3)
    c_Tri_1(1,1:3)=Crack3D_Meshed_Node(i_HF)%row(i_Cr_Node1,1:3)
    c_Tri_1(2,1:3)=Crack3D_Meshed_Node(i_HF)%row(i_Cr_Node2,1:3)
    c_Tri_1(3,1:3)=Crack3D_Meshed_Node(i_HF)%row(i_Cr_Node3,1:3)    

    ! Number of edges of natural fractures.
    num_Vertex = Each_NaCr3D_Poi_Num(i_NF)  

    ! Divide the natural cracks into triangles for assessment
    do i_Tri=1, num_Vertex-2
        c_Tri_2(1,1:3) = Na_Crack3D_Coor(i_NF,1,1:3)
        c_Tri_2(2,1:3) = Na_Crack3D_Coor(i_NF,i_Tri+1,1:3)
        c_Tri_2(3,1:3) = Na_Crack3D_Coor(i_NF,i_Tri+2,1:3)
        ! Check whether the discrete fracture surface triangles of the two cracks intersect.
        call Tool_Intersections_of_Two_Triangles_3D(c_Tri_1,c_Tri_2,c_Logical_Inter,c_num_Inters, &
        c_Inter_Points,c_Logical_Parallel)  
        ! If there is an intersection.
        if(c_Logical_Inter .eqv. .True.)then        
            Yes_Connect = .True.
            num_Inters  = c_num_Inters
            Inter_Points(1:num_Inters,1:3) = c_Inter_Points(1:num_Inters,1:3)
            return
        endif
    enddo
enddo



RETURN
END SUBROUTINE D3_Check_hf_nf_connect
